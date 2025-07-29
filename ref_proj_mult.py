""" Stitch to a bunch of stations using a radius around station

Can simply rename the directory to include 'Stitch' to make it work
"""

import shutil
from collections import OrderedDict
import h5py
from shapely.geometry import Point, Polygon

# fails depending on directory
from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from VLM.bzFRInGE.ref_proj import proj_one

from mintpy.utils import readfile, writefile, utils0, utils as ut

from BZ import bbGIS
from BZ.bbLogger import logger
import geopandas as gpd
gdal.UseExceptions()


def plot_crop(da, geom=None):
    """ Plot the cropped dataset and the geometry used to crop it on map"""
    import cartopy.crs as ccrs
    from cartopy.io import img_tiles as cimgt
    from BZ import bbPlot
    basemap = cimgt.GoogleTiles(
        url='https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}'
        )

    proj = ccrs.PlateCarree()
    WESN = [da.lon.min(), da.lon.max(), da.lat.min(), da.lat.max()]
    # axes.set_extent(WESN, crs=proj)
    fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': basemap.crs})
    # gl = axes.gridlines(draw_labels=True)
    # bbPlot.fmt_gridlines(gl, bottom=True)
    axes.add_image(basemap, 10)

    continuous = True
    norm = mpl.colors.TwoSlopeNorm(0, -5, 5) if continuous else \
        mpl.colors.BoundaryNorm(np.arange(-5, 6, 1), 256)
    pparms = {'cmap': 'cmc.vik', 'norm': norm}
    im   = axes.pcolormesh(da.lon, da.lat, da, shading='nearest',
                    transform=proj, **pparms)
    bbPlot.cartopy_cbar(im, ylabel='Vertical Rate (mm/yr)')
    return


def proj_rad(ExpA, mp_exp, stas, neofs=25, rad=25):
    """" Project in a rad (km) radius around each station """
    epsg    = 'epsg:4326'
    lst_rates, lst_uncs, lst_exps = [], [], []
    for i, sta in enumerate(stas):
        Exp = ExpBase(ExpA, mp_exp, sta, neofs)
        path_vup = proj_one(Exp, npix=0)

        da_rate = h5_to_xr(path_vup, 'velocity')
        da_rate.rio.write_crs(epsg, inplace=True)

        da_unc  = h5_to_xr(path_vup, 'velocityStd')
        da_unc.rio.write_crs(epsg, inplace=True)

        # da_res  = h5_to_xr(path_vup, 'residue')*1000
        # da_res.rio.write_crs(epsg, inplace=True)

        ## clip the rates to polygon
        df_gps   = prep_gps(Exp.path_gps, Exp.reg, horiz=True)
        ser_ref  = df_gps[df_gps.sta == Exp.ref_sta]
        gdf_crop = bbGIS.df2gdf(ser_ref).to_crs(4087).buffer(rad*1000).to_crs(4326)

        da_rate_crop = da_rate.rio.clip(gdf_crop.geometry, epsg, drop=False)
        da_unc_crop  = da_unc.rio.clip(gdf_crop.geometry, epsg, drop=False)
        # da_res_crop  = da_res.rio.clip(gdf_crop.geometry, epsg, drop=False)

        # plot_crop(da_rate_crop*1000)

        lst_rates.append(da_rate_crop)
        lst_uncs.append(da_unc_crop)
        lst_exps.append(Exp)
        if i == 0:
            ExpREF = Exp

    rate_stitch = np.nanmean(lst_rates, 0)
    unc_stitch  = np.nanmean(lst_uncs, 0)

    ix = np.logical_and(~np.isnan(lst_rates[0]), ~np.isnan(lst_rates[1]))

    overlap = np.where(ix, 1, np.nan)

    # do a weighted mean in the overlap region
    # technically should propagate variance of GPS in overlap but
    # same for NYC (0.8 for NJHT, NYBK), overlap in water for HR
    if not np.isnan(overlap).all():
        log.warning('There is overlap between the radius, averaging overlap.')
        # unc_ovl = np.stack([lst_uncs[0].data*overlap, lst_uncs[1].data*overlap])
        unc_ovl = np.stack([lst_unc.data*overlap for lst_unc in lst_uncs])
        wgts     = 1/(unc_ovl**2)
        rate_ovl = np.stack([lst_rate.data*overlap for lst_rate in lst_rates])
        rate_ovl = np.average(rate_ovl, 0, wgts)
        unc_ovl  = np.nanmean(unc_ovl, 0)

        rate_stitch = np.where(np.isnan(rate_ovl), rate_stitch, rate_ovl)
        unc_stitch  = np.where(np.isnan(unc_ovl), unc_stitch, unc_ovl)

    da_rate_crop.data = rate_stitch
    if not socket.gethostname().startswith('leffe'):
        plot_crop(da_rate_crop*1000)

    for Expi in [ExpREF]:
        shutil.copy(Expi.path_vup_geo, Expi.path_vup_geo_stitch)
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, np.flipud(rate_stitch), 'velocity')
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, np.flipud(unc_stitch), 'velocityStd')
        # writefile.write_hdf5_block(Expi.path_vup_geo_stitch, res_stitch, 'residue')

        with h5py.File(Expi.path_vup_geo_stitch, 'r+') as h5:
            h5.attrs['ref_stas'] = stas

        logger.info ('Wrote: %s', Expi.path_vup_geo_stitch)

    log.info('All written _Stitch exp h5 velocities will be the same.')
    log.critical('Need to write netcdfs and kmzs.')

## the intention was to rereference on the fly
## but DEM err corr introduces non-neglible transforms ~0.5 mm/yr that makes just
## subtracting new reference point not work
def proj_rad2(ExpA, mp_exp, stas, neofs=25, rad=25):
    """" Project in a rad (km) radius around each station """
    from mintpy.cli import reference_point
    epsg    = 'epsg:4326'
    lst_rates, lst_uncs, lst_exps = [], [], []

    ## get original GPS experiment
    neofs = 0
    ExpB = ExpBase(ExpA, mp_exp, neofs=neofs)
    ## get the LOS velocity of ExpB
    path_velB = op.join(ExpB.path_mp_exp_geo,
                        f'geo_velocity_recon{ExpB.neofs}.h5') if \
                        ExpB.neofs else ExpB.path_vlos_geo
    ## reference to LOY2
    lalo  = DCT_GPS['LOY2'].split()
    # cmd  = f'{path_velB} --lat {lalo[0]} --lon {lalo[1]}'
    cmd  = f'{path_velB} --lat {lalo[0]} --lon {lalo[1]} -o temp_vel.h5'
    reference_point.main(cmd.split())
    arrB = readfile.read('temp_vel.h5')[0]
    # arrB = readfile.read(path_velB)[0]


    ## get a already projected one using full method and see if you can match
    Exp = ExpBase(ExpA, mp_exp, 'LOY2', neofs)
    # path_vel = '/u/leffe-data2/buzzanga/data/VLM/Sentinel1/HR/MintPy_2alks_5rlks_33_15/LOY2_ERA5_SET_PM_Stitch_ex_Fast/geo/geo_velocity_recon20.h5'
    path_vel = Exp.path_vlos_geo
    arr_loy2 = readfile.read(path_vel)[0]
    print (np.allclose(arr_loy2, arrB, equal_nan=True))
    resid = arr_loy2 - arrB
    print (np.nanmax(resid))
    print (np.nanmin(resid))
    breakpoint()



    for sta in stas:
        breakpoint()
        Exp = ExpBase(ExpA, mp_exp, sta, neofs)
        path_vup = proj_one(Exp, npix=0, path_vel=path_vel)

        da_rate = h5_to_xr(path_vup, 'velocity')
        da_rate.rio.write_crs(epsg, inplace=True)

        da_unc  = h5_to_xr(path_vup, 'velocityStd')
        da_unc.rio.write_crs(epsg, inplace=True)

        # da_res  = h5_to_xr(path_vup, 'residue')*1000
        # da_res.rio.write_crs(epsg, inplace=True)

        ## clip the rates to polygon
        df_gps   = prep_gps(Exp.path_gps, Exp.reg, horiz=True)
        ser_ref  = df_gps[df_gps.sta == Exp.ref_sta]
        da_rate_crop = da_rate.rio.clip(gdf_crop.geometry, epsg, drop=False)
        da_unc_crop  = da_unc.rio.clip(gdf_crop.geometry, epsg, drop=False)
        # da_res_crop  = da_res.rio.clip(gdf_crop.geometry, epsg, drop=False)

        # plot_crop(da_rate_crop*1000)

        lst_rates.append(da_rate_crop)
        lst_uncs.append(da_unc_crop)
        lst_exps.append(Exp)
        if sta == 'SPVA':
            ExpSPVA = Exp

    rate_stitch = np.nanmean(lst_rates, 0)
    unc_stitch  = np.nanmean(lst_uncs, 0)

    ix = np.logical_and(~np.isnan(lst_rates[0]), ~np.isnan(lst_rates[1]))

    overlap = np.where(ix, 1, np.nan)

    # do a weighted mean in the overlap region
    # technically should propagate variance of GPS in overlap but
    # same for NYC (0.8 for NJHT, NYBK), overlap in water for HR
    if not np.isnan(overlap).all():
        log.warning('There is overlap between the radius, averaging overlap.')
        # unc_ovl = np.stack([lst_uncs[0].data*overlap, lst_uncs[1].data*overlap])
        unc_ovl = np.stack([lst_unc.data*overlap for lst_unc in lst_uncs])
        wgts     = 1/(unc_ovl**2)
        rate_ovl = np.stack([lst_rate.data*overlap for lst_rate in lst_rates])
        rate_ovl = np.average(rate_ovl, 0, wgts)
        unc_ovl  = np.nanmean(unc_ovl, 0)

        rate_stitch = np.where(np.isnan(rate_ovl), rate_stitch, rate_ovl)
        unc_stitch  = np.where(np.isnan(unc_ovl), unc_stitch, unc_ovl)

    da_rate_crop.data = rate_stitch
    plot_crop(da_rate_crop*1000)

    ## have to do this for all the stas; prob save the exps above
    for Expi in [ExpSPVA]:
        shutil.copy(Expi.path_vup_geo, Expi.path_vup_geo_stitch)
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, np.flipud(rate_stitch), 'velocity')
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, np.flipud(unc_stitch), 'velocityStd')
        # writefile.write_hdf5_block(Expi.path_vup_geo_stitch, res_stitch, 'residue')

        with h5py.File(Expi.path_vup_geo_stitch, 'r+') as h5:
            h5.attrs['ref_stas'] = stas

        logger.info ('Wrote: %s', Expi.path_vup_geo_stitch)


if __name__ == '__main__':
    # proj_rad(HR_SR, 'ERA5_SET_PM_Stitch_ex_Fast', 'LOY2 SPVA VAHP DRV5 LS03'.split(), neofs=20, rad=25)
    # proj_rad2(HR_SR, 'Base_Fast', 'LOY2 SPVA VAHP DRV5 LS03'.split(), neofs=20, rad=25)

    proj_rad(Houston_SR, 'ERA5_SET_PM_Stitch_ex_Fast',
             'NASA UHDT CFJV SESG TXTG CSTE DMFB TXB6 TXLQ'.split(), neofs=10, rad=20)
    # plt.show()
