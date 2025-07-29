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

# PATH_RES  = op.join(op.expanduser('~'), 'Desktop', 'VLM')

# plt.switch_backend('Qt5Agg')

dct_stable = {
    ## coped from Stitch_ex_Fast
    'Base_ex_Fast' : {
        'LaGuardia_Stable': [(40.782755, -73.853617), 'NYC'],
        'Williams_Stable': [(40.712456, -73.919430), 'NYC'],
        'Ashe_Stable': [(40.75025, -73.825000), 'NYC'],
        'Woodside_Stable': [(40.742738, -73.901400), 'NYC']
        },

    'ERA5_SET_PM_Stitch_Fast' : {
           'LaGuardia_Stable': [(40.747077, -73.895486), 'NYC'],
           'Williams_Stable': [(40.714973, -73.914119), 'NYC'],
           'Ashe_Stable': [(40.746921, -73.825273), 'NYC'],
           'Woodside_Stable': [(40.747474, -73.907443), 'NYC']
           },

    'ERA5_SET_PM_Stitch_Fast2017' : {
           'LaGuardia_Stable': [(40.747077, -73.895486), 'NYC'],
           'Williams_Stable': [(40.714973, -73.914119), 'NYC'],
           'Ashe_Stable': [(40.773866, -73.843051), 'NYC'],
           'Woodside_Stable': [(40.747438, -73.907653), 'NYC']
           },

    'ERA5_SET_PM_Stitch_ex_Fast' : {
           'LaGuardia_Stable': [(40.782755, -73.853617), 'NYC'],
           'Williams_Stable': [(40.712456, -73.919430), 'NYC'],
           'Ashe_Stable': [(40.75025, -73.825000), 'NYC'],
           'Woodside_Stable': [(40.742738, -73.901400), 'NYC']
           },

    'ERA5_SET_PM_Stitch_ex_Fast2017' : {
           'LaGuardia_Stable': [(40.773866, -73.843329), 'NYC'],
        #    'LaGuardia_Stable': [(40.792477, -73.849440), 'NYC'], #option2
           'Williams_Stable': [(40.714653, -73.915824), 'NYC'],
           'Ashe_Stable': [(40.773866, -73.843329), 'NYC'],
           'Woodside_Stable': [(40.742191, -73.898324), 'NYC']
           },

    'ERA5_SET_PM_Stitch_ex90_Fast' : {
           'LaGuardia_Stable': [(40.782755, -73.853607), 'NYC'],
           'Williams_Stable': [(40.714973, -73.914119), 'NYC'],
           'Ashe_Stable': [(40.746921, -73.82527), 'NYC'],
           'Woodside_Stable': [(40.741646, -73.900496), 'NYC']
           },

    'ERA5_SET_PM_Stitch_ex90_Fast2017' : {
           'LaGuardia_Stable': [(40.77387, -73.84305), 'NYC'],
           'Williams_Stable': [(40.714973, -73.914119), 'NYC'],
           'Ashe_Stable': [(40.728588, -73.849440), 'NYC'],
           'Woodside_Stable': [(40.741348, -73.8998151), 'NYC']
           },
    'ERA5_SET_PM_ex_Fast' : {},

}


dct_pts = {
           'LaGuardia':  [(40.77575, -73.866642), 'NYC'],
           'Williams': [(40.7137, -73.921894), 'NYC'],
        #    'Williams_Stable': [(40.714391, -73.915750), 'NYC'],
           'Ashe': [(40.749956, -73.847037), 'NYC'],
           'Woodside': [(40.747195, -73.901), 'NYC'], # 435, 972
           'Fabuwood': [(40.73942, -74.12917), 'NYC'],

           'Pink Frog': [(40.717786, -73.954012), 'NYC'],
           'Williamsburg': [(40.707459, -73.956934), 'NYC'],
           'East Williamsburg 1': [(40.706125, -73.935909), 'NYC'],
           'Crown Heights': [(40.661058, -73.945637), 'NYC'],
           'Flatlands': [(40.622270, -73.940868), 'NYC'],
           'Kew Garden': [(40.718839, -73.819076), 'NYC'],
           'Kew Garden1': [(40.711215, -73.813457), 'NYC'],
           'Soho':  [( 40.724199, -74.002366), 'NYC'],
           'Jersey City': [(40.722614, -74.054235), 'NYC'],
           'Bayonne':  [(40.661856, -74.108103), 'NYC'],
           'BayonneS': [(40.663171,  -74.107424), 'NYC'],
           'Westside': [(32.791415, -79.950970), 'Charleston'],
           'NYBK': [(40.7034315045, -73.9789651762), 'NYC'],

           }


def check_rate(rate_true, rate_ts, npix, error=False):
    """ Check rates calculate by hand (rate_ts) match those in velocity file """
    if npix == 0:
        print (f'"Real velocity": {rate_true:.2f}')
        print (f'New velocity: {rate_ts:.2f}')
        print (f'Difference between new and "real" velocity: '\
            f'{np.abs(rate_ts - rate_true):.2f} mm/yr')

        if error:
            assert np.abs(rate_true - rate_ts) < 1e-3, \
            'Velocity estimated from time series does not match velocity from uvel'

    else:
        print (f'"Real velocity": {rate_true:.2f}')
        print (f'New velocity: {rate_ts:.2f}')
        print (f'Difference between new and "real" velocity: '\
               f'{np.abs(rate_ts - rate_true):.2f} mm/yr')
    return


def buffer_point(lat, lon, rad):
    """ Roughly buffer a lalo point given a radius of meters

    Return S, N, W, E bounding box in WGS84
    """
    from shapely.geometry import Point
    import geopandas as gpd
    gser = gpd.GeoSeries(Point(lon, lat), crs='EPSG:4326')
    gdf  = gser.to_crs('EPSG:4087').buffer(rad).to_crs('EPSG:4326')
    W, S, E, N = gdf.squeeze().bounds
    return S, N, W, E


def plot_ts_vup1(Exp, loc, lalo=None, npix=0, ref_sta=None, mask=True, save=False):
    """ Plot the vertical displacement timeseries average (meters) around one loc """
    if 'recon' in Exp.path_ts_geo:
        log.info('Using reconstructed time series...')
        path_ts_geo  = Exp.path_ts_geo
    else:
        log.info('Using time series with annual cycle removed...')
        path_ts_geo  = Exp.path_ts_geo_ann

    path_vup_geo = Exp.path_vup_geo
    if 'stitch' in Exp.mp_exp.lower():
        if Exp.reg != 'NYC':
            raise Exception ('Stitching not working for multiple stations (ok for NYC)')
        log.info('Using stitched time series and velocity...')
        path_ts_geo  = Exp.path_ts_geo_stitch
        path_vup_geo = Exp.path_vup_geo_stitch


    print ('Reading velocity and uncertainty from file:', path_vup_geo)
    vel_true  = readfile.read(path_vup_geo)[0]*1000
    unc, meta = readfile.read(path_vup_geo, datasetName='velocityStd')
    arr_iang  = readfile.read(Exp.path_geom_mp_geo, datasetName='incidenceAngle')[0]
    unc *= 1000

    try:
        tsObj  = timeseries(path_ts_geo)
        arr_ts = tsObj.read()
    except:
        log.warning('Could not get requested stitched timeseries...')
        tsObj  = timeseries(Exp.path_ts_geo)
        arr_ts = tsObj.read()
    dates  = pd.to_datetime(tsObj.get_date_list())
    decyr  = bbTS.date2dec(dates)

    if lalo is None:
        lalo = dct_stable[Exp.mp_exp0][loc][0] if 'stable' in loc.lower() else dct_pts[loc][0]

    if 'PM' in Exp.mp_exp.upper():
        arr_pm, meta1 = readfile.read(op.join(Exp.path_mp_exp_geo, 'ITRF14.h5'))
        arr_pm = arr_pm.astype(np.float32)
        assert 'REF_X' in meta1.keys(), 'Need to reference the plate motion'
        arr_pm0= (1000 * arr_pm / np.cos(np.deg2rad(arr_iang)).astype(np.float32)).reshape(-1)
        inters = np.zeros_like(arr_pm0, dtype=np.float32)
        # arr_pm = np.polyval([arr_pm.reshape(-1), np.zeros_like(arr_pm.reshape(-1))], decyr)
        arr_pm1 = (arr_pm0[:, np.newaxis] * decyr) + inters[:, np.newaxis]
        arr_pm  = arr_pm1.transpose(1, 0).reshape(arr_ts.shape)
        del arr_pm1, inters

    else:
        arr_pm = 0

    ## convert whole timeseries to Vup, subtract potential plate motion (mm/yr)
    arr  = 1000 * arr_ts / np.cos(np.deg2rad(arr_iang))
    arr -= arr_pm

    if mask:
        with h5py.File(Exp.path_mask_vup, 'r') as h5:
            mask = np.where(np.isclose(h5['waterMask'][:], 0), np.nan, 1)
            arr *= mask
            # arr1 *= np.where(np.isclose(mask, 0), np.nan, 1)
    else:
        mask = 1

    # average around the reference point (bbox or radius)
    # ref_lalo = [float(meta['REF_LAT']), float(meta['REF_LON'])]
    ref_y, ref_x  = int(meta['REF_Y']), int(meta['REF_X'])

    ## the pixel at the actual reference GPS is 0; can average around it a bit
    # ref_ts  = arr[:, ref_y, ref_x]
    # ref_ts  = np.nanmean(arr[:, ref_y-1:ref_y+1, ref_x-1:ref_x+1], axis=(1, 2))
    # ref_ts1  = avg_ts(arr, ref_lalo, rad, Exp.path_rate_nc)

    ## optionally average a bit around the target point
    coord    = ut.coordinate(meta, lookup_file=Exp.path_geom_mp_geo)
    y, x     = coord.geo2radar(lalo[0], lalo[1])[0:2]
    # loc_ts   = avg_ts(arr, lalo, rad, Exp.path_rate_nc)
    if npix == 0:
        loc_ts   = arr[:, y, x] # no real diff if surrounding points similar
        vel_true = vel_true[y, x]
        unc      = unc[y, x]
        assert not np.isnan(mask[y,x]), f'{loc} Pixel y={y} x={x} is masked!'
    else:
        loc_ts   = np.nanmean(arr[:, y-npix:y+npix, x-npix:x+npix], axis=(1, 2))
        vel_true = np.nanmean(vel_true[y-npix:y+npix, x-npix:x+npix])
        unc      = np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])

    # DONT remove the timeseries at the GPS reference from the target point
        # remember timeseries is not shifted to GPS
        # but is referenced, so theoretically ref_ts should be 0
        # this mostly adds noise so I'M NOT DOING IT

    # ts_los   = loc_ts + ref_ts
    # ts_los  = ts_los - loc_ts[0]
    ts_los  = loc_ts - loc_ts[0]

    df_gps  = prep_gps(Exp.path_gps, Exp.reg).set_index('sta')
    ref_sta = Exp.ref_sta if ref_sta is None else ref_sta
    gps_vel = df_gps.loc[ref_sta, 'u_vel'] * 1000
    gps_sig = df_gps.loc[ref_sta, 'u_sig'] * 1000

    # make a 'timeseries' from the GPS velocity
    ts_gps  = np.polyval([gps_vel, 0], decyr)
    ts_gps -= ts_gps[0]

    ## shift vertical timeseries by GPS timeseries
    ts_vup     = ts_los + ts_gps
    coeffs_vup, cov_vup = np.polyfit(decyr, ts_vup, 1, cov='unscaled')
    # preds, coeffs_vup = bbTS.bzSeason.fit_lsq(decyr, ts_vup) # seasonal

    rss = np.polyfit(decyr, ts_vup, 1, full=True)[1] # residual sum of squares
    rmse = np.sqrt(rss.item() / len(ts_vup)) # mm

    log.info (f'Loc {loc} y/x: %s/%s', y, x)
    check_rate(vel_true, coeffs_vup[0], npix)

    trend_vup  = np.polyval(coeffs_vup, decyr)

    # trend_vup  = preds[:, 0] ## show seasonal

    ## the polyfit is too low
    # sem2       = np.sqrt(cov_vup[0,0]+gps_sig**2)
    # sem2 = np.sqrt(unc**2+gps_sig**2) # this double counting GPS!
    sem2 = unc


    ## plot the Vup
    fig, axes = plt.subplots(figsize=(10, 4))

    col = 'r' if coeffs_vup[0] > 0 else 'b'
    axes.scatter(dates, ts_vup, c='k', s=15)#, label=f'$\sigma$={np.std(ts_vup):.2f}')
    # axes.scatter(dates, ts_vup, alpha=0, label=f'$RMSE$={rmse:.2f}')

    axes.plot(dates, trend_vup, color=col, linestyle='--', label=f'{coeffs_vup[0]:.2f}$\pm${sem2:.2f} (mm/yr)')
    axes.fill_between(dates, trend_vup-sem2, trend_vup+sem2, color=col, alpha=0.3)
    axes.set_ylabel('Vertical Displacement (mm)')
    axes.grid(color='gray', linestyle = '--', linewidth=0.1)
    axes.legend()

    if loc in dct_pts:
        axes.set_title(f'{loc}, {dct_pts[loc][1]}')
    elif loc in dct_stable[Exp.mp_exp0]:
        axes.set_title(f'{loc}, {dct_stable[Exp.mp_exp0][loc][1]}')
    else:
        axes.set_title(f'{loc}')

    fig.set_label(f'{loc}_ts')
    if save:
        path_figs = op.join(PATH_RES, f'{Exp.reg}_2024')
        bbPlot.savefigs(path_figs, True, True)

    del arr, arr_pm
    return fig, axes


def plot_ts_vup1_dd(Exp, loc_target, loc_stable, npix=0, mask=True, axes=None, show_ts=False):
    """ Plot the vertical displacement timeseries average (npix) around one loc

    Doubled differenced; plot 'loc_target' after removing 'loc_stable'
    Note that plate motion may give subtle differences in stable trend

    loc_target/stable can be a single name (gets lalo from dictionary) or a tuple with name / lalo
    """
    from mintpy.utils import readfile, utils as ut
    ## add the houston points to the dictionary
    if 'recon' in Exp.path_ts_geo:
        log.info('Using reconstructed time series...')
        path_ts_geo  = Exp.path_ts_geo
    else:
        log.info('Using time series with annual cycle removed...')
        path_ts_geo  = Exp.path_ts_geo_ann

    path_vup_geo = Exp.path_vup_geo
    if 'stitch' in Exp.mp_exp.lower():
        if Exp.reg != 'NYC':
            raise Exception ('Stitching not working for multiple stations (ok for NYC)')
        log.info('Using stitched time series and velocity...')
        path_ts_geo  = Exp.path_ts_geo_stitch
        path_vup_geo = Exp.path_vup_geo_stitch

    print (f'Reading velocity and uncertainty from file for {loc_target}:', path_vup_geo)

    vup_true  = readfile.read(Exp.path_vup_geo)[0]*1000
    unc, meta = readfile.read(Exp.path_vup_geo, datasetName='velocityStd')
    arr_iang  = readfile.read(Exp.path_geom_mp_geo, datasetName='incidenceAngle')[0]
    unc *= 1000

    ## can pass tuple to locs with name and lalo
    try:
        lalo_target = dct_pts[loc_target][0]
        lalo_stable = dct_stable[Exp.mp_exp0][loc_stable][0]
    except:
        loc_target, lalo_target = loc_target
        loc_stable, lalo_stable = loc_stable

    ## do this now to save a lot of time
    mask   = readfile.read(Exp.path_mask_vup)[0].astype(float) if mask else 1
    coord  = ut.coordinate(meta, lookup_file=Exp.path_geom_mp_geo)
    y, x   = coord.geo2radar(lalo_target[0], lalo_target[1])[0:2]

    if isinstance(mask, np.ndarray):
        assert mask[y,x], f'{loc_target} Pixel y={y} x={x} is masked!'

    try:
        tsObj  = timeseries(path_ts_geo)
        arr_ts = tsObj.read()
    except:
        log.warning('Could not get requested stitched timeseries...')
        tsObj  = timeseries(Exp.path_ts_geo)
        arr_ts = tsObj.read()
    dates  = pd.to_datetime(tsObj.get_date_list())
    decyr  = bbTS.date2dec(dates)

    ## convert whole timeseries to Vup
    arr   = 1000 * arr_ts / np.cos(np.deg2rad(arr_iang)) * mask

    if 'PM' in Exp.mp_exp.upper():
        arr_pm, meta1 = readfile.read(op.join(Exp.path_mp_exp_geo, 'ITRF14.h5'))
        arr_pm = arr_pm.astype(np.float32)
        assert 'REF_X' in meta1.keys(), 'Need to reference the plate motion'
        arr_pm0= (1000 * arr_pm / np.cos(np.deg2rad(arr_iang))).reshape(-1)
        inters = np.zeros_like(arr_pm0, dtype=np.float32)
        # arr_pm = np.polyval([arr_pm.reshape(-1), np.zeros_like(arr_pm.reshape(-1))], decyr)
        arr_pm1 = (arr_pm0[:, np.newaxis] * decyr) + inters[:, np.newaxis]
        arr_pm  = arr_pm1.transpose(1, 0).reshape(arr_ts.shape)
        arr_pm -= arr_pm[0] # reference to first date
        arr    -= arr_pm


    # ## average a bit around the reference point
    # attrs    = readfile.read_attribute(Exp.path_ts_geo)
    # ref_ts   = avg_ts(arr, ref_lalo, rad, Exp.path_rate_nc)

    ## average a bit around the target point
    log.info (f'{loc_target} y/x: %s/%s', y, x)
    if npix == 0:
        ts_vup_target   = arr[:, y, x] # no real diff if surrounding points similar
        vup_true_target = vup_true[y, x]
        unc_target      = unc[y, x]
    else:
        ts_vup_target   = np.nanmean(arr[:, y-npix:y+npix, x-npix:x+npix], axis=(1, 2))
        vup_true_target = np.nanmean(vup_true[y-npix:y+npix, x-npix:x+npix])
        unc_target      = np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])
        assert not np.isnan(vup_true_target), 'Pixel is masked despite averaging {npix} pixels!'

    df_gps    = prep_gps(Exp.path_gps, Exp.reg, units='mm').set_index('sta')
    ref_vel   = df_gps.loc[Exp.ref_sta, 'u_vel']
    log.debug('This may be the incorrect ref station depending on loc of ts point')

    rate_target = np.polyfit(decyr, ts_vup_target, 1)[0]
    check_rate(vup_true_target, rate_target+ref_vel, npix)

    # use the stable point exactly
    y, x     = coord.geo2radar(lalo_stable[0], lalo_stable[1])[0:2]
    log.info (f'{loc_stable} y/x: %s/%s', y, x)

    ts_vup_stable   = arr[:, y, x] # no real diff if surrounding points similar
    vup_true_stable = vup_true[y, x]
    unc_stable      = unc[y, x]
    if isinstance(mask, np.ndarray):
        assert mask[y,x], f'{loc_stable} Pixel y={y} x={x} is masked!'


    # else:
    #     ts_vup_stable   = np.nanmean(arr[:, y-npix:y+npix, x-npix:x+npix], axis=(1, 2))
    #     vup_true_stable = np.nanmean(vup_true[y-npix:y+npix, x-npix:x+npix])
    #     unc_stable      = np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])

    rate_stable = np.polyfit(decyr, ts_vup_stable, 1)[0]
    check_rate(vup_true_stable, rate_stable+ref_vel, npix)

    ## double difference
    ts_dd    = ts_vup_target - ts_vup_stable
    # ser = pd.Series(ts_dd, index=dates, color='')
    # import seaborn as sns
    #         ser = pd.series(ts_dd, index=dates)
    # coeffs_dd = np.polyfit(decyr, ts_dd, 1)
    # trend_dd  = np.polyval(coeffs_dd, decyr)
    # sns.lmplot(data=df, x='decyr', y='ts', scatter=True,
            #    scatter_kws={'color': 'k'}, line_kws={'color': 'red'})

    ## show the two timeseries
    if show_ts:
        fig, axes1 = plt.subplots(figsize=(10, 4))
        axes1.scatter(dates, ts_vup_target, label=f'{loc_target}', color='darkblue', s=15, alpha=0.25)
        axes1.scatter(dates, ts_vup_stable, label=f'{loc_stable}', color='darkgreen', s=15, alpha=0.25)
        axes1.scatter(dates, ts_dd, label=f'DDiff', color='black', s=15, alpha=0.75)
        # axes1.plot(dates, trend_dd, linestyle='--', label=f'DDiff: {coeffs_dd[0]:.2f} mm/yr', color='black')
        axes1.set_title('Vertical')
        axes1.legend()
        axes1.set_title(f'Target/Stable Rate: {rate_target:.2f} / {rate_stable:.2f}')
        print ('Remember, numbers in title arent shifted by GPS')

    ## remove the timeseries at loc0 from loc and the GPS reference

    ## get GPS velocity and project it to a timeseries with 0 intercept
    ## not used (each GPS is same for ts and ts stable point)
    ## will break for two reference stations


    ## use the reference point offset ; may be different for stitched
    # gps_vel = df_gps.loc[Exp.ref_sta, 'u_vel']

    if 'stitch' in Exp.mp_exp.lower():
        ## get the nearest gps to the target point
        import geopandas as gpd
        try:
            ref_sta = df_gps[df_gps.index.isin(Exp.ref_sta)]
        except:
            ref_sta = df_gps[df_gps.index == Exp.ref_sta]
        gdf_target = gpd.GeoDataFrame(crs=4326,
            geometry=gpd.points_from_xy([lalo_target[1]], [lalo_target[0]]))
        gdf = bbGIS.df2gdf(ref_sta, 'all').to_crs(4087)
        gdf_near = gdf_target.to_crs(4087).sjoin_nearest(gdf,
                                        distance_col='dist').to_crs(4326)
        gps_vel = ref_sta[ref_sta.index.isin(gdf_near.index_right)].u_vel.item()
        gps_sig = ref_sta[ref_sta.index.isin(gdf_near.index_right)].u_sig.item()

    ## mean uncertainty of the two stat
    elif not isinstance(Exp.ref_sta, str) and len(Exp.ref_sta) > 1:
        log.warning('Mean GPS velocity/uncertainty of two stations')
        ref_sta = df_gps[df_gps.index.isin(Exp.ref_sta)]
        gps_vel = np.sqrt(np.mean(ref_sta.u_sig ** 2))
        gps_sig = np.mean(ref_sta.u_vel)

    ## these two are same, just if ref sta is list/str
    elif not isinstance(Exp.ref_sta, str) and len(Exp.ref_sta) == 1:
        ref_sta = df_gps[df_gps.index.isin(Exp.ref_sta)]
        gps_vel = ref_sta.u_vel.item()
        gps_sig = ref_sta.u_sig.item()

    else:
        ref_sta = df_gps[df_gps.index == Exp.ref_sta]
        gps_vel = ref_sta.u_vel.item()
        gps_sig = ref_sta.u_sig.item()

    # gps_sig = df_gps.loc[Exp.ref_sta, 'u_sig'] * 1000

    # ts_gps  = np.polyval([gps_vel, 0], decyr)
    # ts_gps -= ts_gps[0]

    ## this matches the Vup rate so we're good; careful with annual, plate motion
    # ts1 = arr[:, 435, 972]
    # ts_vup1 = ts1 + ts_gps
    # coeffs_vup, cov_vup = np.polyfit(decyr, ts_vup1, 1, cov='unscaled')
    # print (coeffs_vup[0])

    ## DONT shift vertical timeseries by GPS timeseries
    # theoretically each pixel shifted to GPS above but then it cancels

    # ts_vup_ref = ts_dd + ts_gps
    log.info (f'Vup Rate Target: {rate_target+gps_vel:.2f}, '\
            f'Vup Rate Stable: {rate_stable+gps_vel:.2f}')

    # for Emi, maybe something else?
    df  = pd.DataFrame({f'{loc_target}_dd': ts_dd, f'{loc_target}_stable': ts_vup_stable,
                       f'{loc_target}_target': ts_vup_target, 'decyr': decyr, 'date':dates})
    corr  = True if not 'Base' in Exp.mp_exp0 else False
    corrl = '_corr' if corr else '_raw'
    df['corr'] = corr
    dst = f'{PATH_RES}/{loc_target}_dd{corrl}.csv'
    df.to_csv(dst)
    log.info(f'Wrote: %s', dst)

    rate_dd, unc_dd = bbTS.bootstrap_1D(ts_dd, decyr, n_boots=1500)
    rss = np.polyfit(decyr, ts_dd, 1, full=True)[1]
    rmse = np.sqrt(rss/len(ts_dd)).item()

    # sem2       = np.sqrt(cov_vup[0,0]+gps_sig**2)
    sem2       = np.sqrt(unc_dd**2+gps_sig**2)
    coeffs_vup, cov_vup = np.polyfit(decyr, ts_dd, 1, cov='unscaled')
    trend_vup  = np.polyval(coeffs_vup, decyr)
    # check_ts_unc(dates, ts_dd) # other experiments

    # using seasonal
    # preds, coeffs_vup = bbTS.bzSeason.fit_lsq(decyr, ts_dd)
    # coeffs_vup = coeffs_vup[-2:]
    # trend_vup = preds[:, 0]
    # rate_dd  = coeffs_vup[0]


    ## plot the Vup
    if axes is None:
        fig, axes = plt.subplots(figsize=(10, 4))
        axes.set_title(f'{loc_target}-{loc_stable}')#-GPS')
    else:
        fig = axes.get_figure()

    col = 'r' if rate_dd > 0 else 'b'
    axes.scatter(dates, ts_dd, c='k', s=15)
    axes.plot(dates, trend_vup, color=col, linestyle='--',
              label=f'{rate_dd:.1f}$\pm${sem2:.1f} (mm/yr)', alpha=0)
    # axes.fill_between(dates, trend_vup-sem2, trend_vup+sem2, color=col, alpha=0.3)
    axes.set_ylabel('Vertical Displacement (mm)', fontsize=13)
    lloc = 'lower right' if loc_target == 'Williams' else 'best'
    axes.legend(prop ={"size": 13}, loc=lloc)
    axes.tick_params(axis='x', labelsize=13)
    fig.set_label(f'{loc_target}_dd_{npix}pix')

    log.info(f'{loc_target} (npix={npix})')

    log.info ('Original Rate, Timeseries Std:, %.1f, %.1f',
                vup_true_target, np.nanstd(ts_vup_target))

    log.info ('Double Diff Rate, Timeseries Std, RMSE: %.1f, %.1f, %.1f',
                coeffs_vup[0], np.nanstd(ts_dd), rmse)

    # axes.set_title(f'Timeseries Std; RMSE; npix: {np.nanstd(ts_dd):.1f}; {rmse:.1f}; {npix}')

    # axes.set_ylim([-40, 10])

    axes.grid(color='gray', linestyle = '--', linewidth=0.1)
    fig.set_label(f'tsdd_{loc_target}_{npix}')
    del arr, arr_pm
    return fig, axes


def plot_woodside_dd(Exp,  npix=0, mask=True, axes=None, split='20200101'):
    """ Plot the vertical displacement timeseries average (npix)

    at woodside with two trends
    """
    from mintpy.utils import readfile, utils as ut
    loc_target = 'Woodside'
    loc_stable = 'Woodside_Stable'

    vup_true  = readfile.read(Exp.path_vup_geo)[0]*1000
    unc, meta = readfile.read(Exp.path_vup_geo, datasetName='velocityStd')
    arr_iang  = readfile.read(Exp.path_geom_mp_geo, datasetName='incidenceAngle')[0]
    unc *= 1000
    df_gps    = prep_gps(Exp.path_gps, Exp.reg, units='mm').set_index('sta')
    ref_vel   = df_gps.loc[Exp.ref_sta, 'u_vel']
    log.debug('This may be the incorrect ref station depending on loc of ts point')

    # path_ts_geo = Exp.path_ts_geo if 'recon' in Exp.path_ts_geo else Exp.path_ts_geo_ann
    # log.info('Not using seasonally corrected')
    path_ts_geo = Exp.path_ts_geo if 'recon' in Exp.path_ts_geo else Exp.path_ts_geo
    if 'stitch' in Exp.mp_exp.lower():
        log.warning('Using stitched, recon time series')
        path_ts_geo = Exp.path_ts_geo_stitch


    with h5py.File(path_ts_geo, 'r') as h5:
        arr_ts = h5['timeseries'][:]
        arr_dt = [dt.decode('utf-8') for dt in h5['date'][:]]

    dates = pd.to_datetime(arr_dt)
    decyr = bbTS.date2dec(dates)

    lalo_target = dct_pts[loc_target][0]
    lalo_stable = dct_stable[Exp.mp_exp0][loc_stable][0]

    ## convert whole timeseries to Vup
    arr   = 1000 * arr_ts / np.cos(np.deg2rad(arr_iang))

    if 'PM' in Exp.mp_exp.upper():
        arr_pm, meta1 = readfile.read(op.join(Exp.path_mp_exp_geo, 'ITRF14.h5'))
        assert 'REF_X' in meta1.keys(), 'Need to reference the plate motion'
        arr_pm0= (1000 * arr_pm / np.cos(np.deg2rad(arr_iang))).reshape(-1)
        inters = np.zeros_like(arr_pm0)
        # arr_pm = np.polyval([arr_pm.reshape(-1), np.zeros_like(arr_pm.reshape(-1))], decyr)
        arr_pm1 = (arr_pm0[:, np.newaxis] * decyr) + inters[:, np.newaxis]
        arr_pm  = arr_pm1.transpose(1, 0).reshape(arr_ts.shape)
        arr_pm -= arr_pm[0] # reference to first date
        arr    -= arr_pm

    if mask:
        mask = readfile.read(Exp.path_mask_vup)[0]
        arr *= np.where(np.isclose(mask, 0), np.nan, 1)

    ## average a bit around the target point
    coord    = ut.coordinate(meta, lookup_file=Exp.path_geom_mp_geo)
    y, x     = coord.geo2radar(lalo_target[0], lalo_target[1])[0:2]
    log.info ('Loc y/x: %s/%s', y, x)
    if npix == 0:
        ts_vup_target   = arr[:, y, x] # no real diff if surrounding points similar
        vup_true_target = vup_true[y, x]
        unc_target      = unc[y, x]
    else:
        ts_vup_target   = np.nanmean(arr[:, y-npix:y+npix, x-npix:x+npix], axis=(1, 2))
        vup_true_target = np.nanmean(vup_true[y-npix:y+npix, x-npix:x+npix])
        unc_target      = np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])

    rate_target = np.polyfit(decyr, ts_vup_target, 1)[0]

    # Check rates calculate by hand match those in velocity file
    check_rate(vup_true_target, rate_target+ref_vel, npix)

    ## get coords of (stable) point
    y, x     = coord.geo2radar(lalo_stable[0], lalo_stable[1])[0:2]

    # use the stable point exactly
    ts_vup_stable   = arr[:, y, x] # no real diff if surrounding points similar
    vup_true_stable = vup_true[y, x]
    unc_stable      = unc[y, x]

    rate_stable = np.polyfit(decyr, ts_vup_stable, 1)[0]
    check_rate(vup_true_stable, rate_stable+ref_vel, npix)

    ## double difference
    ts_dd    = ts_vup_target - ts_vup_stable

    ## get gps for propagating the uncertainty
    if 'stitch' in Exp.mp_exp.lower():
        ## get the nearest gps to the target point
        import geopandas as gpd
        try:
            ref_sta = df_gps[df_gps.index.isin(Exp.ref_sta)]
        except:
            ref_sta = df_gps[df_gps.index == Exp.ref_sta]
        gdf_target = gpd.GeoDataFrame(crs=4326,
            geometry=gpd.points_from_xy([lalo_target[1]], [lalo_target[0]]))
        gdf = bbGIS.df2gdf(ref_sta, 'all').to_crs(4087)
        gdf_near = gdf_target.to_crs(4087).sjoin_nearest(gdf,
                                        distance_col='dist').to_crs(4326)
        gps_vel = ref_sta[ref_sta.index.isin(gdf_near.index_right)].u_vel.item()
        gps_sig = ref_sta[ref_sta.index.isin(gdf_near.index_right)].u_sig.item()

    ## mean uncertainty of the two stat
    elif not isinstance(Exp.ref_sta, str) and len(Exp.ref_sta) > 1:
        ref_sta = df_gps[df_gps.index.isin(Exp.ref_sta)]
        gps_vel = np.mean(ref_sta.u_vel)
        gps_sig = np.sqrt(np.mean(ref_sta.u_sig ** 2))

    elif not isinstance(Exp.ref_sta, str) and len(Exp.ref_sta) == 1:
        ref_sta = df_gps[df_gps.index.isin(Exp.ref_sta)]
        gps_vel = ref_sta.u_vel.item()
        gps_sig = ref_sta.u_sig.item()

    else:
        ref_sta = df_gps[df_gps.index == Exp.ref_sta]
        gps_vel = ref_sta.u_vel.item()
        gps_sig = ref_sta.u_sig.item()

    #  note pm may have been removed, and this is relative to GPS
    log.info (f'Vup Rate Target: {rate_target+gps_vel:.2f}, '\
              f'Vup Rate Stable: {rate_stable+gps_vel:.2f}')

    # for Emi, maybe something else?
    df  = pd.DataFrame({f'{loc_target}_dd': ts_dd, f'{loc_target}_stable': ts_vup_stable,
                       f'{loc_target}_target': ts_vup_target, 'decyr': decyr, 'date':dates})
    corr  = True if not 'Base' in Exp.mp_exp0 else False
    corrl = '_corr' if corr else '_raw'
    df['corr'] = corr
    dst = f'{PATH_RES}/{loc_target}_dd{corrl}.csv'
    df.to_csv(dst)
    log.info(f'Wrote: %s', dst)

    rate_dd, unc_dd = bbTS.bootstrap_1D(ts_dd, decyr, n_boots=1500)
    rss = np.polyfit(decyr, ts_dd, 1, full=True)[1]
    rmse = np.sqrt(rss/len(ts_dd)).item()

    # sem2       = np.sqrt(cov_vup[0,0]+gps_sig**2)
    sem2       = np.sqrt(unc_dd**2+gps_sig**2)
    coeffs_vup, cov_vup = np.polyfit(decyr, ts_dd, 1, cov='unscaled')
    trend_vup  = np.polyval(coeffs_vup, decyr)


    ## plot the Vup
    if axes is None:
        fig, axes = plt.subplots(figsize=(10, 4))
        axes.set_title(f'{loc_target}-{loc_stable}')#-GPS')
    else:
        fig = axes.get_figure()

    ts0, dates0 = ts_dd, dates
    ix1, ix2 = dates0<=split, dates0>split
    # ix1 *= dates0>='20170101'
    colors = 'crimson gray'.split()
    for i, ix in enumerate([ix1, ix2]):
        dates, ts = dates0[ix], ts0[ix]
        decyr     = bbTS.date2dec(dates)
        unc_dd    = bbTS.bootstrap_1D(ts, decyr, n_boots=1500)[1]
        # hack to add the GPS uncertainty
        sem2      = np.sqrt(unc_dd**2 + 0.76**2)
        coeffs    = np.polyfit(decyr, ts, 1)
        rate      = coeffs[0]
        if np.abs(rate) < 1e-2:
            rate = np.abs(rate)
        trend     = np.polyval(coeffs, decyr)
        axes.scatter(dates, ts, c='k', s=15)
        l, = axes.plot(dates, trend, color=colors[i], linestyle='--',
                       label=f'{rate:.1f}$\pm${unc_dd:.1f} (mm/yr)')
        axes.fill_between(dates, trend-sem2, trend+sem2, color=l.get_color(), alpha=0.3)

        # l1, = axes.fill(np.nan, np.nan, l.get_color(), alpha=0.3)

        axes.legend(prop ={"size": 13}, loc='lower right')
        # axes.legend([(l1, l)], )

    # axes.scatter(dates, ts_dd, c='k', s=15)
    # axes.plot(dates, trend_vup, color='r', linestyle='--', label=f'{rate_dd:.1f}$\pm${sem2:.1f} (mm/yr)')
    # axes.fill_between(dates, trend_vup-sem2, trend_vup+sem2, color='r', alpha=0.3)
    # axes.set_ylabel('Vertical Displacement (mm)', fontsize=13)
    # axes.legend(prop ={"size": 12})
    # axes.tick_params(axis='x', labelsize=13)
    # fig.set_label(f'{loc_target}_dd')

    log.info(f'{loc_target} (npix={npix})')
    # log.info ('Double Diff Rate, Timeseries Std, RMSE: %.1f, %.1f, %.1f',
    #             coeffs_vup[0], np.nanstd(ts_dd), rmse)

    axes.set_title(f'Timeseries Std; RMSE; npix: {np.nanstd(ts_dd):.1f}; {rmse:.1f}; {npix}')
    axes.set_ylabel('Vertical Displacement (mm)', fontsize=13)

    axes.grid(color='gray', linestyle = '--', linewidth=0.1)
    fig.set_label(f'tsdd_{loc_target}_{npix}_split')
    return fig, axes


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
    axes.scatter(dates, ts_vup, c='k', s=15, label=f'$\sigma$={np.std(ts_vup):.2f}')
    axes.plot(dates, trend_vup, color='r', linestyle='--', label=f'{coeffs_vup[0]:.2f}$\pm${sem2:.2f} (mm/yr)')
    axes.fill_between(dates, trend_vup-sem2, trend_vup+sem2, color='r', alpha=0.3)
    # axes.set_ylabel('Vertical Displacement (mm)')
    axes.grid(color='gray', linestyle = '--', linewidth=0.1)
    axes.legend()
    axes.set_title(loc)
    fig.set_label(f'{loc}_ts_poly')
    return


def check_ts_unc(dates, ts):
    """ Alternate calculations for computing double difference time series rate/unc"""
    from statsmodels.formula.api import ols
    import statsmodels.api as smv
    decyr = bbTS.date2dec(dates)
    df = pd.DataFrame({'decyr': decyr, 'date': dates, 'ts': ts})
    mod = ols('ts ~ decyr', data=df).fit()
    print (mod.summary())

    coeffs = np.polyfit(decyr, ts, 1)
    trend  = np.polyval(coeffs, decyr)
    residuals = (ts - trend)
    RMSE   = np.sqrt(residuals**2)
    print (f'{RMSE:.2f} mm/yr')

    breakpoint()


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
                     label=f'{coeffs[0]:.1f}$\pm${unc_dd:.1f} (mm/yr)')
        axes[i].legend()

    axes[0].set_title(f'Split Date: {str(pd.to_datetime(split).date())}')


def find_stable(Exp, lalo, rad=50, mask=True, alpha=0.8, loc=None, plot=True):
    """ Find a stable point around another one. Loc just for naming

    Radius in meters around target point
    """
    import geopandas as gpd
    from shapely.geometry import Point, Polygon
    from mintpy.utils import readfile
    # get the insar rates
    epsg    = 'epsg:4326'
    da_rate = h5_to_xr(Exp.path_vup_geo, 'velocity')
    da_rate.rio.write_crs(epsg, inplace=True)

    if mask:
        mask = np.flipud(readfile.read(Exp.path_mask_vup)[0])

        da_rate *= mask
        da_rate  = da_rate.where(~np.isclose(da_rate, 0), np.nan)

    lat, lon = lalo
    # look around the target point
    gser = gpd.GeoSeries(gpd.points_from_xy([lon], [lat]), crs=epsg)
    poly = gser.to_crs(4087).buffer(rad, cap_style=1).to_crs(epsg).geometry
    da_pt = da_rate.rio.clip(poly, epsg, drop=True)

    # first cut of mostly high values; can cut again later
    thresh = 1.0 / 1000 # keep only values +/- thresh m
    da_pt  = da_pt.where(np.abs(da_pt) < thresh, np.nan)
    da_pt *= 1000 # convert to mm
    # da_pt = da_pt.where(np.abs(da_pt) < 0.001, np.nan)
    if da_pt.isnull().all():
        print ('No valid points! Try increasing search radius or values')
        return

    df = da_pt.to_dataframe().reset_index().abs().sort_values('velocity').dropna()
    df.lon *= -1
    gdf = bbGIS.df2gdf(df, 'all')

    # get distance from target point to potential
    gdf_pt0 = gpd.GeoDataFrame(geometry=gpd.points_from_xy([lon], [lat], crs=4326))
    gdf_m  = gpd.sjoin_nearest(gdf.to_crs(4087), gdf_pt0.to_crs(4087),
                                distance_col='distance')
    center = gdf_m.dissolve().centroid.to_crs(4326)
    gdf_m.to_crs(4326, inplace=True)


    ## get time series noise ---
    arr_ts, meta = readfile.read(Exp.path_ts_geo)
    arr_ts *= np.flipud(mask)
    coord = ut.coordinate(meta, lookup_file=Exp.path_geom_mp_geo)

    ts_std = []
    for ix, row in gdf.iterrows():
        y, x = coord.geo2radar(row.lalo.y, row.lalo.x)[0:2]
        ts_std.append(np.nanstd(arr_ts[:, y, x]*1000))

    gdf_m['ts_std'] = ts_std
    if loc is not None:
        gdf_m.index = [loc] * gdf_m.shape[0]

    if plot:
        import cartopy.crs as ccrs
        from cartopy.io import img_tiles as cimgt
        alpha=1

        basemap   = cimgt.GoogleTiles(url='https://server.arcgisonline.com/'\
                        'arcgis/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}.jpg')
        proj    = ccrs.PlateCarree()
        S,N,W,E    = bbGIS.get_extent_nc(da_pt)
        norm = mpl.colors.TwoSlopeNorm(0)#, -thresh, thresh)

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': basemap.crs})
        axes.add_image(basemap, 10)
        axes.set_extent([W-0.2,E+0.2,S-0.2,N+0.2], crs=proj)

        axes.scatter(lon, lat, alpha=alpha, transform=proj, color='k', marker='s', s=1)
        axes.scatter(-center.x, center.y, alpha=alpha, transform=proj, color='red', marker='o', s=5)
        im   = axes.pcolormesh(da_pt.lon, da_pt.lat, da_pt, shading='nearest', alpha=alpha,
                        transform=proj, cmap='cmc.roma_r', norm=norm)
        bbPlot.cartopy_cbar(im, 'V$_{up}$ (mm/yr)')
        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)

        # dont bother getting ix; netcdf is flipped; mintpy coord is right
        # sub = da_rate.sel(lat=40.714391, lon=-73.915750, method='nearest')
        # y, x  = np.where(da_rate == sub)
        # 1446, 1179

        # da_rate.data[y, x]

        breakpoint()

    return gdf_m


## ----------- Wrappers
def check_stable(Exp1):
    assert Exp1.reg == 'NYC', print ('Only NYC supported')

    for pt in 'LaGuardia Williams Woodside Ashe'.split():
        plot_ts_vup1(Exp1, f'{pt}_Stable')
        # calculate the distance
        print (f'Distance to stable pt {pt}', np.round(
                        bbGIS.distance_pt2pt(dct_pts[pt][0],
                        dct_stable[Exp1.mp_exp0][f'{pt}_Stable'][0],
                        units='km'), 2), 'km')


    # plt.close('all')


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

    Exp0  = ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=10)

    df_wells = pd.read_csv(op.join(Exp0.path_wd, 'Houston_wells.txt'),
                delim_whitespace=True, header=None, index_col=0, names='lat lon'.split())

    lat, lon = df_wells.iloc[0].to_numpy()
    loc = df_wells.iloc[0]

    lst_gdfs = []
    for sta, row in df_wells.iterrows():
        gdf_m = find_stable(Exp0, [lat, lon], rad=10000, loc=sta.replace('R6-TX-', ''), plot=False)
        lst_gdfs.append(gdf_m)

    gdf_ms = pd.concat(lst_gdfs).drop(columns='crs index_right'.split())
    gdf_ms.index.rename('name', inplace=True)

    dst = op.join(Exp0.path_wd, 'stable_oil.GeoJSON')
    gdf_ms.to_file(dst)
    print (f'Wrote: {dst}')


    # check_stable(Exp0)



    # bbPlot.savefigs(PATH_RES, True, True)


    plt.show()
