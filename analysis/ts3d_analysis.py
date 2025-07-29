""" Some basic functions for analyzing 3D time series """
from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase

from BZ import bbPlot

import h5py
import xarray as xr
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cmcrameri, cmocean

from subprocess import call
import shutil

BASEMAP = cimgt.GoogleTiles(url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
PROJ    = {'projection': ccrs.PlateCarree()}
# OCEAN = cfeature.NaturalEarthFeature('physical', 'ocean', \
#             scale='10m', edgecolor='black', facecolor='black')
# LAND  = cfeature.NaturalEarthFeature('physical', 'land', \
#             scale='10m', edgecolor='#4d4d4d', facecolor='#4d4d4d')

DCT_WELLS = {'well0':
    'https://cida.usgs.gov/ngwmn/provider/USGS/site/405101073343401/'
    }


def calc_eofs(da, neofs=5):
    """" Quick EOF calculation, no fancy options like in CESM_EOFs """
    from eofs.xarray import Eof
    coslat = np.cos(np.deg2rad(da.lat.data))

    # numerical issues cause negatives instead of 0s which raise annoying warnings
    coslat = np.where(coslat < 0, 0, coslat)
    wgts   = np.sqrt(coslat)[..., np.newaxis]

    solver = Eof(da, weights=wgts, center=True)
    eofs  = solver.eofs(neofs=neofs, eofscaling=2) # in original units
    pcs   = solver.pcs(npcs=neofs, pcscaling=1)    # scaled to unit variance
    per_var = solver.varianceFraction()*100
    return eofs, pcs, per_var


def plot_rate(Exp, vmin=-3, vmax=3):
    """ Just plot the vertical velocity """
    da   = xr.open_dataset(Exp.path_rate_msk_nc)['Band1']*1000

    norm      = mpl.colors.TwoSlopeNorm(0, vmin, vmax)
    fig, axes = plt.subplots(figsize=(12,9), subplot_kw=PROJ)
    im   = axes.pcolormesh(da.lon, da.lat, da, shading='auto',
                norm=norm, transform=ccrs.PlateCarree(), cmap='cmc.roma_r')
    axes.add_image(BASEMAP, 11)

    gl = axes.gridlines(draw_labels=True)
    _fmt_gridlines(gl, bottom=True)
    cbar = bbPlot.cartopy_cbar(im, xlabel='mm/yr')
    axes.set_title(f'Vup')


def plot_eofs_pcts(da, neofs=5, vmin=-5, vmax=5):
    """ Plot rate and a bunch of eofs spatial patterns and their pcts """
    eofs, pcs, per_var = calc_eofs(da, neofs)
    eofs = [eofs] if not 'xarray' in str(type(eofs)) else eofs

    print (f'First {neofs} modes contain {per_var[:neofs].sum().item():.2f}% of the variance\n')

    norm = mpl.colors.TwoSlopeNorm(0, vmin, vmax)
    figs = []
    for eof in eofs:
        eof = eof.where(eof != 0, np.nan)
        fig, axes = plt.subplots(figsize=(12,9), subplot_kw=PROJ)
        im   = axes.pcolormesh(eof.lon, eof.lat, eof, shading='auto',
                    norm=norm, transform=ccrs.PlateCarree(), cmap='cmo.tarn')
        axes.add_image(BASEMAP, 11)
        gl = axes.gridlines(draw_labels=True)
        _fmt_gridlines(gl, bottom=True)
        cbar = bbPlot.cartopy_cbar(im, xlabel=da.units)
        axes.set_title(f'EOF Mode: {eof.mode.item()+1}')
        figs.append(fig)

    fig, axes = plt.subplots(figsize=(12,9), nrows=neofs, sharex=True, sharey=True)
    for ax, pc in zip(axes, pcs.T):
        ax.plot(da.time, pc, color='k')
        ax.set_ylabel(f'PCTS: {pc.mode.item()+1} (none)')

    # fig.set_label(Exp.reg

    return figs


def recon_Vup(Exp, da, neofs):
    """ Reconstruct the timeseries using a subset of the EOFs

    Write to timeseries file
    Calculate velocity with MintPy
    Reference and project
    """
    from eofs.xarray import Eof
    coslat = np.cos(np.deg2rad(da.lat.data))

    # numerical issues cause negatives instead of 0s which raise annoying warnings
    coslat = np.where(coslat < 0, 0, coslat)
    wgts   = np.sqrt(coslat)[..., np.newaxis]

    solver = Eof(da, weights=wgts, center=True)
    recon  = solver.reconstructedField(neofs)

    # stick it into a new hdf5

    dst    = op.join(op.dirname(Exp.path_ts_geo), f'geo_timeseries_recon{neofs}.h5')
    shutil.copy(Exp.path_ts_geo, dst)

    with h5py.File(dst, 'r+') as h5:
        del h5['timeseries-ann']
        data    = h5['timeseries']
        data[:] = np.fliplr(recon.data/1000)
        h5.attrs['FILE_PATH'] = dst


    print ('Wrote:', dst)
    src = dst
    dst = op.join(op.dirname(Exp.path_ts_geo), f'geo_vlos_recon{neofs}.h5')

    cmd  = f'timeseries2velocity.py {src} --ref-lalo {DCT_GPS[Exp.ref_sta]} '
    cmd += f'-o {dst}'

    call(cmd.split())

    per_var = solver.varianceFraction()*100
    print (f'First {neofs} modes contain {per_var[:neofs].sum().item():.2f}% of the variance\n')

    from contrib.ref_proj import RefProj

    lbl = f'{Exp.lbl}_recon{neofs}'
    Obj = RefProj(Exp.path_mp_exp_geo, Exp.reg, Exp.path_mask_mp_nc, corr=lbl)

    if Exp.reg == 'HR':
        Obj.proj_HR(radius=150)
    #     Obj.proj_one(ref_sta)
    #     Obj.proj_two()

    else:
        Obj.proj_one(Exp.ref_sta)


    return


## not implemented
def plot_epoch(ts, dates, epoch=-1, cmap='cmc.roma', vmin=-30, vmax=30):
    """ Plot one epoch of ts as a test """
    ts1 = ts[epoch]
    dt1 = dates[epoch]

    fig, axes = plt.subplots(figsize=(10,10))
    axes.imshow(ts1, interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
    return


## think this is replace by Jupyter
def load_ts_nc(Exp, dataset='timeseries-ann', mask=True, overwrite=False):
    dst     = f'{op.splitext(Exp.path_ts_geo)[0]}.nc'
    if op.exists(dst) and not overwrite:
        da  =  xr.open_dataset(dst)[dataset]
        print ('Got:', dst)
        return da

    from mintpy.utils import readfile
    with h5py.File(Exp.path_ts_geo, 'r') as h5:
        arr_ts  = h5['timeseries'][:]
        arr_tsa = h5['timeseries-ann'][:]
        arr_dt  = [dt.decode('utf-8') for dt in h5['date'][:]]

    dates = pd.to_datetime(arr_dt)
    arr0  = arr_ts*1000
    arr1  = arr_tsa*1000

    if mask:
        with h5py.File(Exp.path_mask_vup, 'r') as h5:
            mask = h5['mask'][:]
            arr0 *= np.where(np.isclose(mask, 0), np.nan, 1)
            arr1 *= np.where(np.isclose(mask, 0), np.nan, 1)

    attrs = readfile.read_attribute(Exp.path_ts_geo)
    if attrs['ORBIT_DIRECTION'] == 'ASCENDING':
        arr0 = np.fliplr(arr0)
        arr1 = np.fliplr(arr1)

    ## for matching the geom
    ds     = xr.open_dataset(Exp.path_rate_nc)
    das    = []
    for arr, name in zip([arr0, arr1], 'timeseries timeseries-ann'.split()):
        das.append(xr.DataArray(arr, name=name, dims=['time', 'lat', 'lon'],
                        coords={'time': dates, 'lat': ds.lat, 'lon': ds.lon}
                        ).assign_attrs(units='mm'))

    ds = xr.merge(das)
    ds.to_netcdf(dst)
    print ('Wrote timeseries to netcdf:', dst)

    return ds[dataset]






if __name__ == '__main__':
    Exp   = load_exp(Charleston_SR, 'Base', 'SCHA')
    da_ts = load_ts_nc(Exp)
    # plot_rate(Exp)
    # plot_eofs_pcts(da_ts)
    recon_Vup(Exp, da_ts, 1)

    plt.show()
