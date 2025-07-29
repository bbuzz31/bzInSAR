#!/usr/bin/env python3
import shutil
from collections import OrderedDict
import h5py
from shapely.geometry import Point, Polygon
from BZ import bbGIS

# fails depending on directory
from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase

from mintpy.utils import readfile, writefile, utils0, utils as ut

from BZ.bbLogger import logger
import geopandas as gpd
gdal.UseExceptions()


# same as in ref_proj_mult
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

## same as in ref_proj.py
def proj_one(ExpA, path_vel=None, npix=0, na12=False):
    """ Tie to a single GPS station (optionally npix surrounding pixels) """
    ## identical to the other class version
    if path_vel is None:
        path_vel = op.join(ExpA.path_mp_exp_geo,
                           f'geo_velocity_recon{ExpA.neofs}.h5') if \
                            ExpA.neofs else ExpA.path_vlos_geo

    logger.info ('Projection velocity at: %s', path_vel)

    ## remove the station from the exclude when using e.g., stitching
    df_gps   = prep_gps(ExpA.path_gps, ExpA.reg, horiz=True)
    ser_ref  = df_gps[df_gps.sta == ExpA.ref_sta]

    assert not ser_ref.empty, f'No station matching {ExpA.ref_sta} in {ExpA.path_gps}'

    logger.info ('Tieing to: %s', ser_ref.sta.item())

    rate, meta = readfile.read(path_vel, datasetName='velocity')
    unc  = readfile.read(path_vel, datasetName='velocityStd')[0]
    inc  = readfile.read(ExpA.path_geom_mp_geo, datasetName='incidenceAngle')[0]
    az   = readfile.read(ExpA.path_geom_mp_geo, datasetName='azimuthAngle')[0]

    # 73, 470
    y, x     = int(meta['REF_Y']), int(meta['REF_X'])
    if na12:
        rate_ref = utils0.enu2los(ser_ref.e_vel, ser_ref.n_vel, ser_ref.u_vel, inc[y, x], az_angle=az[y,x]).item()
        unc_ref  = utils0.enu2los(ser_ref.e_sig, ser_ref.n_sig, ser_ref.u_sig, inc[y, x], az_angle=az[y,x]).item()
        path_vup_geo = ExpA.path_vup_geo_na12

    ## residue will be copied over if it exists
    # res  = readfile.read(path_vel, datasetName='residue')[0]
    else:
        rate /= np.cos(np.deg2rad(inc))
        unc  /= np.cos(np.deg2rad(inc)) # bumps it up slightly
        rate_ref = ser_ref.u_vel.item()
        unc_ref  = ser_ref.u_sig.item()
        path_vup_geo = ExpA.path_vup_geo

    ### average a few pixels around the GPS station (at GPS, insar=0)
    if npix:
        rate_ref = np.nanmean(rate[y-npix:y+npix, x-npix:x+npix]) + rate_ref
        unc_ref  = np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])  + unc_ref

    logger.info (f'GPS Offset: {rate_ref*1000:.2f} +/- {unc_ref*1000:.2f} mm/yr')
    rate += rate_ref
    unc   = np.sqrt(unc**2 + unc_ref**2)

    # use na12 gps
    if na12:
        rate /= np.cos(np.deg2rad(inc))
        unc  /= np.cos(np.deg2rad(inc)) # bumps it up slightly

    # copy and update layers
    shutil.copy(path_vel, path_vup_geo)
    writefile.write_hdf5_block(path_vup_geo, rate, 'velocity')
    writefile.write_hdf5_block(path_vup_geo, unc, 'velocityStd')
    # writefile.write_hdf5_block(ExpA.path_vup_geo, res, 'residue')

    # add reference station to attributes

    with h5py.File(path_vup_geo, 'r+') as h5:
        h5.attrs['ref_stas'] = ser_ref.sta.item()

    logger.info ('Wrote: %s', path_vup_geo)
    return path_vup_geo


## same as proj_NYC
def proj_HR2(ExpA, ExpB, npix=0, proj=True):
    """ Project to north/south and stitch together

    Optionally skip projection (set to False) for use within another method
    """
    epsg     = 'epsg:4326'
    lst_rates, lst_uncs, lst_ress = [], [], []
    for Expi in [ExpA, ExpB]:
        if proj:
            path_vup = proj_one(Expi, npix=npix)
        else:
            path_vup = Expi.path_vup_geo
        da_rate = h5_to_xr(path_vup, 'velocity')
        da_rate.rio.write_crs(epsg, inplace=True)

        da_unc  = h5_to_xr(path_vup, 'velocityStd')
        da_unc.rio.write_crs(epsg, inplace=True)

        # da_res  = h5_to_xr(path_vup, 'residue')*1000
        # da_res.rio.write_crs(epsg, inplace=True)

        ## clip the rates to polygon
        path_crop    = op.join(Expi.path_crops, f'{CROP_MAP[Expi.ref_sta]}.GeoJSON')
        gdf_crop     = gpd.read_file(path_crop)
        da_rate_crop = da_rate.rio.clip(gdf_crop.geometry, epsg, drop=False)
        da_unc_crop  = da_unc.rio.clip(gdf_crop.geometry, epsg, drop=False)
        # da_res_crop  = da_res.rio.clip(gdf_crop.geometry, epsg, drop=False)

        lst_rates.append(da_rate_crop)
        lst_uncs.append(da_unc_crop)

        logger.info (f'Mean Velocity in {Expi.ref_sta} {CROP_MAP[Expi.ref_sta]} '\
               f'{da_rate_crop.mean().item()*1000:.2f} mm/yr')

        # lst_ress.append(da_res_crop)

    rate_stitch = np.nanmean(lst_rates, 0)
    unc_stitch  = np.nanmean(lst_uncs, 0)

    # overlap = np.sum(lst_rates, 0)
    ix = np.logical_and(~np.isnan(lst_rates[0]), ~np.isnan(lst_rates[1]))
    overlap = np.where(ix, 1, np.nan)

    # do a weighted mean in the overlap region
    # technically should propagate variance of GPS in overlap but
    # same for NYC (0.8 for NJHT, NYBK), overlap in water for HR
    if not np.isnan(overlap).all():
        log.warning('There is overlap between the polygons, averaging overlap.')
        unc_ovl = np.stack([lst_uncs[0].data*overlap, lst_uncs[1].data*overlap])
        wgts     = 1/(unc_ovl**2)
        rate_ovl = np.stack([lst_rates[0].data*overlap, lst_rates[1].data*overlap])
        rate_ovl = np.average(rate_ovl, 0, wgts)
        unc_ovl  = np.nanmean(unc_ovl, 0)

        rate_stitch= np.where(np.isnan(rate_ovl), rate_stitch, rate_ovl)
        unc_stitch= np.where(np.isnan(unc_ovl), unc_stitch, unc_ovl)

    ## need to be able to handle stitched (multiple) and single stations
    ref_stas  = []
    for Expi in [ExpA, ExpB]:
        meta = readfile.read_attribute(Expi.path_vup_geo)
        ref_stas0 = meta.get('ref_stas', '')
        ref_stas1 = [ref_stas0] if isinstance(ref_stas0, str) else ref_stas0
        ref_stas.extend(ref_stas1)

    ## just make both; they'll be the same
    for Expi in [ExpA, ExpB]:
        shutil.copy(Expi.path_vup_geo, Expi.path_vup_geo_stitch)
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, np.flipud(rate_stitch), 'velocity')
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, np.flipud(unc_stitch), 'velocityStd')
        # writefile.write_hdf5_block(Expi.path_vup_geo_stitch, res_stitch, 'residue')
        with h5py.File(Expi.path_vup_geo_stitch, 'r+') as h5:
            h5.attrs['ref_stas'] = ref_stas

        logger.info ('Wrote: %s', Expi.path_vup_geo_stitch)
    return


def proj_HRn(ExpA, stas=['LS03'], npix=0):
    """ Use a single reference and then to other stations with offset using other GPS

    Called by 2020
    """
    df_gps   = prep_gps(ExpA.path_gps, ExpA.reg, horiz=True)

    lst_rates, lst_uncs = [], []

    ## eventually write these to the individual MintPy directory
    path_vel = op.join(ExpA.path_mp_exp_geo, f'geo_velocity_recon{ExpA.neofs}.h5') \
        if ExpA.neofs else ExpA.path_vlos_geo

    rate, meta = readfile.read(path_vel, datasetName='velocity')
    unc  = readfile.read(path_vel, datasetName='velocityStd')[0]
    inc  = readfile.read(ExpA.path_geom_mp_geo, datasetName='incidenceAngle')[0]
    y, x  = int(meta['REF_Y']), int(meta['REF_X'])
    coord = ut.coordinate(meta, lookup_file=ExpA.path_geom_mp_geo)

    ## convert it to LOS
    rate /= np.cos(np.deg2rad(inc))
    unc  /= np.cos(np.deg2rad(inc)) # bumps it up slightly

    # now iterate over the ref and extra station to get an offset
    offsets, weights = [], []
    stas = [ExpA.ref_sta, stas] if isinstance(stas, str) else \
        np.insert(stas, 0, ExpA.ref_sta).tolist()
    for i, sta in enumerate(stas):
        ser_sta  = df_gps[df_gps.sta == sta]
        rate_ref = ser_sta.u_vel.item()
        unc_ref  = ser_sta.u_sig.item()
        # unc_ref will be 0 at the real reference station, so then it's propagated later

        ## get the other GPS coordinate
        if i > 0:
            y, x = coord.geo2radar(ser_sta.lat, ser_sta.lon)[0:2]

        ## use the difference between the GPS rate and the pixels as offset
        # this is correct
        if npix:
            offset = rate_ref - np.nanmean(rate[y-npix:y+npix, x-npix:x+npix])
            weight = unc_ref**2 + np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])**2
        else:
            offset = rate_ref - rate[y, x]
            weight = unc_ref**2 + unc[y,x]**2

        ## use the abs difference between GPS and pixels as weight
        offsets.append(offset)
        weights.append(np.abs(weight))

    ## ~~~~~ PERFORM WEIGHTED MEAN OF OFFSETS USING GNSS ONLY ~~~~~ ##
    wgts       = 1/np.array(weights) # variance
    offset     = np.average(offsets, weights=wgts) # equivalent to polyfit
    offset_unc = np.average(np.sqrt(weights), weights=wgts)

    logger.info (f'Fit Offset: {offset*1000:.2f} +/- {offset_unc*1000:.2f} mm/yr')

    rate += offset
    unc   = np.sqrt(unc**2 + offset_unc**2)

    shutil.copy(path_vel, ExpA.path_vup_geo)
    writefile.write_hdf5_block(ExpA.path_vup_geo, rate, 'velocity')
    writefile.write_hdf5_block(ExpA.path_vup_geo, unc, 'velocityStd')

    with h5py.File(ExpA.path_vup_geo, 'r+') as h5:
        h5.attrs['ref_stas'] = stas

        print (stas)

    logger.info ('Wrote: %s', ExpA.path_vup_geo)

    ## pretty sure this is here by mistake
    # epsg    = 'epsg:4326'
    # da_rate = h5_to_xr(ExpA.path_vup_geo, 'velocity')
    # da_rate.rio.write_crs(epsg, inplace=True)

    # da_unc  = h5_to_xr(ExpA.path_vup_geo, 'velocityStd')
    # da_unc.rio.write_crs(epsg, inplace=True)
    # logger.info ('Wrote: %s', ExpA.path_vup_geo)

    return ExpA.path_vup_geo


def proj_HR2020(ExpN, ExpS, stas=['SPVA','LS03'], npix=2):
    """ Project using the method of Buzzanga et al., 2020 """
    ## project the north
    # writes ExpN.path_vup_geo
    path_vup_n = proj_one(ExpN, npix=npix)

    ## project the south  using the weighted offset method
    # writes ExpS.path_vup_geo
    path_vup_s = proj_HRn(ExpS, stas, npix)

    # writes ExpN/S.path_vup_geo_stitch
    proj_HR2(ExpN, ExpS, proj=False)
    return




if __name__ == '__main__':
    # Exp0       = ExpBase(HR_SR, 'Base', 'SPVA', neofs=20)
    # proj_one(Exp0, npix=0)

    # Exp0 = ExpBase(HR_SR, 'ERA5_SET_PM_Stitch1_ex_Fast', 'VAHP', neofs=25)
    # Exp1 = ExpBase(HR_SR, 'ERA5_SET_PM_Stitch1_ex_Fast', 'LOY2', neofs=25)
    # proj_HR2(Exp0, Exp1)

    # proj_rad(HR_SR, 'ERA5_SET_PM_Stitch_ex_Fast', 'LOY2 SPVA VAHP DRV5 LS03'.split(), neofs=20, rad=25)
    # proj_rad2(HR_SR, 'ERA5_SET_PM_Stitch_ex_Fast', 'LOY2 SPVA VAHP DRV5 LS03'.split(), neofs=20, rad=25)
    proj_rad2(HR_SR, 'Base_Fast', 'LOY2 SPVA VAHP DRV5 LS03'.split(), neofs=20, rad=25)
    # plt.show()

