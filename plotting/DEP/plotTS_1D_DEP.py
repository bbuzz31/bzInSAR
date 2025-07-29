"""
Plot a single vertical timeseries
Should be a class
"""

import h5py
import xarray as xr
from VLM.bzFRInGE import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from BZ import bbTS, bbPlot, bbGIS
from mintpy.utils import readfile, utils as ut
from mintpy.objects import timeseries


log = logging.getLogger('BZ')
log.setLevel('DEBUG')

# PATH_RES  = op.join(op.expanduser('~'), 'Desktop', 'VLM')

### misc tools
def plot_ts_poly(Expi, loc='LaGuardia'):
    """ Plot the vertical displacement timeseries averaged using a polygon

    Current polys are LaGuardia and Williams; Drawn in GEarth then
    ogr2ogr dst.GeoJSON src.kml

    This has neglible effect on the standard and damps the rates too much
    """
    import geopandas as gpd
    arr_iang, meta  = readfile.read(Expi.path_geom_mp_geo, datasetName='incidenceAngle')

    da_iang = h5_to_xr(Expi.path_geom_mp_geo, 'incidenceAngle')

    if 'stitch' in Expi.mp_exp.lower():
        log.warning('Using stitched, recon time series')
        path_ts_geo = Expi.path_ts_geo_nc_stitch

        ds_rate = xr.open_dataset(Expi.path_rate_nc_stitch)
        ds_unc  = xr.open_dataset(Expi.path_std_nc_stitch)
        da_ts   = xr.open_dataset(path_ts_geo)['timeseries']

    else:
        raise Exception('Not implemented')
        # path_ts_geo = Expi.path_ts_geo_nc if 'recon' in Expi.path_ts_geo \
        #     else Expi.path_ts_geo_ann_nc

        ds_rate  = xr.open_dataset(Expi.path_rate_nc)
        ds_unc   = xr.open_dataset(Expi.path_std_nc)
        da_ts    = xr.open_dataset(path_ts_geo)['timeseries']

    da_rate = ds_rate['Band1']
    da_unc  = ds_unc['Band1']

    ## convert whole timeseries to Vup (mm/yr)
    da_ts /= np.cos(np.deg2rad(da_iang))
    dates = pd.to_datetime(da_ts['time'])
    decyr = bbTS.date2dec(dates)

    ## this makes the agreement worse; maybe something with shape
    if 'PM' in Expi.mp_exp.upper():
        arr_pm, meta1 = readfile.read(op.join(Expi.path_mp_exp_geo, 'ITRF14.h5'))
        assert 'REF_X' in meta1.keys(), 'Need to reference the plate motion'
        arr_pm0= (arr_pm / np.cos(np.deg2rad(arr_iang))).reshape(-1)
        inters = np.zeros_like(arr_pm0)
        # arr_pm = np.polyval([arr_pm.reshape(-1), np.zeros_like(arr_pm.reshape(-1))], decyr)
        arr_pm1 = (arr_pm0[:, np.newaxis] * decyr) + inters[:, np.newaxis]
        arr_pm  = arr_pm1.transpose(1, 0).reshape(da_ts.shape)
        arr_pm -= arr_pm[0] # reference to first date
        da_ts.data -= np.fliplr(arr_pm)


    epsg=4326
    da_rate.rio.write_crs(epsg, inplace=True)
    da_unc.rio.write_crs(epsg, inplace=True)
    da_ts.rio.write_crs(epsg, inplace=True)

    path_crop    = op.join(Expi.path_crops, f'{loc}.GeoJSON')
    gdf_crop     = gpd.read_file(path_crop)
    da_rate_crop = da_rate.rio.clip(gdf_crop.geometry, epsg, drop=True)
    da_unc_crop = da_unc.rio.clip(gdf_crop.geometry, epsg, drop=True)
    da_ts_crop = da_ts.rio.clip(gdf_crop.geometry, epsg, drop=True)

    rate1 = da_rate_crop.mean().item()*1000
    unc1 = da_unc_crop.mean().item()*1000
    ts   = da_ts_crop.mean('lat lon'.split())*1000


    # now get the stable point
    # lalo     = dct_pts[f'{loc}_Stable'][0]
    # coord    = ut.coordinate(meta, lookup_file=Exp.path_geom_mp_geo)
    # y, x     = coord.geo2radar(lalo[0], lalo[1])[0:2]
    # ts_stable = da_ts.sel(lat=lalo[0], lon=lalo[1], method='nearest')
    # ts -= ts_stable

    ts_los  = ts - ts[0]

    df_gps  = prep_gps(Expi.path_gps, Expi.reg).set_index('sta')
    log.warning('Using reference sta: %s', Expi.ref_sta)
    gps_vel = df_gps.loc[Expi.ref_sta, 'u_vel'] * 1000
    gps_sig = df_gps.loc[Expi.ref_sta, 'u_sig'] * 1000

    # make a 'timeseries' from the GPS velocity
    ts_gps  = np.polyval([gps_vel, 0], decyr)
    ts_gps -= ts_gps[0]

    ## shift vertical timeseries by GPS timeseries
    ts_vup     = ts_los + ts_gps
    coeffs_vup, cov_vup = np.polyfit(decyr, ts_vup, 1, cov='unscaled')

    check_rate(rate1, coeffs_vup[0], npix=da_rate_crop.size)

    trend_vup  = np.polyval(coeffs_vup, decyr)

    ## the polyfit error is too low
    # sem2       = np.sqrt(cov_vup[0,0]+gps_sig**2)
    sem2 = np.sqrt(unc1**2+gps_sig**2)

    ## plot the Vup
    fig, axes = plt.subplots(figsize=(10, 4))
    axes.scatter(dates, ts_vup, c='k', s=15, label=f'$\\sigma$={np.std(ts_vup):.2f}')
    axes.plot(dates, trend_vup, color='r', linestyle='--', label=f'{coeffs_vup[0]:.2f}$\\pm${sem2:.2f} (mm/yr)')
    axes.fill_between(dates, trend_vup-sem2, trend_vup+sem2, color='r', alpha=0.3)
    # axes.set_ylabel('Vertical Displacement (mm)')
    axes.grid(color='gray', linestyle = '--', linewidth=0.1)
    axes.legend()
    axes.set_title(loc)
    fig.set_label(f'{loc}_ts_poly')
    return


def plot_ts_vup1_dd_split(Exp1, split='20200101'):
    """ Split the time series because it looks like there's a jump"""
    # fig, axes = plot_ts_vup1_dd(Exp1, 'Williams', 'Williams_Stable')
    fig, axes = plot_ts_vup1_dd(Exp1, 'Woodside', 'Woodside_Stable')
    plt.close('all')

    dates0 = pd.to_datetime(axes.get_lines()[0].get_xdata())
    ts0    = axes.collections[0].get_offsets()[:, 1].data

    ix1, ix2 = dates0<=split, dates0>split
    # ix1 *= dates0>='20170101'


    fig, axes = plt.subplots(nrows=2, figsize=(10,4))
    for i, ix in enumerate([ix1, ix2]):
        dates, ts = dates0[ix], ts0[ix]
        decyr     = bbTS.date2dec(dates)
        unc_dd    = bbTS.bootstrap_1D(ts, decyr, n_boots=1500)[1]
        # hack to add the GPS uncertainty
        unc_dd    = np.sqrt(unc_dd**2 + 0.76**2)
        coeffs    = np.polyfit(decyr, ts, 1)
        trend     = np.polyval(coeffs, decyr)
        axes[i].scatter(dates, ts, c='k', s=15)
        axes[i].plot(dates, trend, color='r', linestyle='--',
                     label=f'{coeffs[0]:.1f}$\\pm${unc_dd:.1f} (mm/yr)')
        axes[i].legend()

    axes[0].set_title(f'Split Date: {str(pd.to_datetime(split).date())}')


## ----------- Wrappers

def concat_ts_csv():
    """ One csv of the indiviudal ts and their stable/target points for emi """
    df0  = pd.read_csv(f'{PATH_RES}/Woodside_dd_raw.csv', index_col=0)
    # df0['loc'] = 'woodside'
    df1  = pd.read_csv(f'{PATH_RES}/Woodside_dd_corr.csv', index_col=0)
    # df1['loc'] = 'woodside'
    df2  = pd.read_csv(f'{PATH_RES}/Williams_dd_raw.csv', index_col=0)
    # df2['loc'] = 'williams'
    df3  = pd.read_csv(f'{PATH_RES}/Williams_dd_corr.csv', index_col=0)
    # df3['loc'] = 'williams'
    df_m1 = pd.concat([df0['date decyr'.split()], df2.filter(like='Williams'),
                       df0.filter(like='Woodside'), df0['corr']], axis=1)
    df_m2 = pd.concat([df1['date decyr'.split()], df3.filter(like='Williams'),
                       df1.filter(like='Woodside'), df1['corr']], axis=1)
    df_m1.columns = df_m1.columns.str.lower()
    df_m2.columns = df_m2.columns.str.lower()
    df_m1.to_csv(f'{PATH_RES}/ts_raw.csv')
    df_m2.to_csv(f'{PATH_RES}/ts_corr.csv')
    print (f'Done concatenating and writing to {PATH_RES}')
    return


def run_NYC():
    """ Wrapper for testing """
    # loc0, loc1 = 'LaGuardia LaGuardia_Stable'.split()
    loc0, loc1 = 'Williams Williams_Stable'.split()
    # loc0, loc1 = 'Ashe Ashe_Stable'.split()
    # loc0, loc1 = 'Woodside Woodside_Stable'.split()

    Exp0     = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20)
    plot_ts_vup1(Exp0, 'LaGuardia')
    plot_ts_vup1(Exp0, 'Williams')
    # plot_ts_vup1(Exp0, 'Ashe_Stable')
    plot_ts_vup1_dd(Exp0, loc0, loc1, npix=3, show_ts=False)
    plot_woodside_dd(Exp0, npix=0)
    # plot_ts_vup1_dd_split(Exp0, '20200101')

    # concat_ts_csv()

    # for npix in range(10):
        # plot_ts_vup1_dd(Exp0, 'LaGuardia', 'LaGuardia_Stable', npix=npix, show_ts=False)

    # find_stable(Exp0, dct_pts['LaGuardia'][0], rad=3000)
    # find_stable(Exp0, dct_pts['Williams'][0], rad=1000)
    # find_stable(Exp0, dct_pts['Ashe'][0], rad=3000)
    # find_stable(Exp0, dct_pts['Woodside'][0], rad=1000)
    # check_stable(Exp0)


if __name__ == '__main__':
    # run_NYC()

    Exp0  = ExpBase(DC_SR, 'ERA5_SET_PM_ex_Fast', 'USN7', neofs=20)

    # plot_ts_vup1(Exp0, 'Navy_Yard')
    # plot_ts_vup1(Exp0, 'Anacostia')
    # plot_ts_vup1(Exp0, 'MD_Up')
    # plot_ts_vup1(Exp0, 'MD_Down')
    plot_ts_vup1(Exp0, 'Anacostia_Stable')
    # find_stable(Exp0, dct_pts['Anacostia'][0])

    # bbPlot.savefigs(PATH_RES, True, True)


    plt.show()
