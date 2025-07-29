from VLM import *
from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from VLM.bzFRInGE.plotting.plot_TS1 import PlotTS1
from VLM.bzFRInGE.plotting.plot_TS1 import DCT_PTS
from BZ import bbPlot, bbTS

""" Use with GW_Analysis Notebook """

## would rather not use this; use my file instead
def main0(site_no, TSObj, npix=5, service='dv', freq='MS'):
    """ service = iv, dv (instaneous vs daily avg); gw not supported yet """
    region = TSObj.reg
    # get the well info
    ser_sta = nwis.get_record(sites=site, service='site').squeeze()
    lalo = (ser_sta.geometry.y, ser_sta.geometry.x)
    DCT_PTS['well'] = [lalo, region]

    # get the well timeseries
        # Data-value qualification codes included in this output:
            # A  Approved for publication -- Processing and review completed.
            # P  Provisional data subject to revision.

    k1 = '_Mean' if service == 'dv' else ''

    # df1 = nwis.get_record(sites=site, service='iv', start='1950-01-01')
    df1 = nwis.get_record(sites=site, service=service, start='2013-01-01')
    df1.drop('site_no', axis=1, inplace=True)
    n_prov = (df1[f'72019{k1}_cd'] == 'P').sum()
    print (f"{n_prov} ({100*n_prov/df1.shape[0]:.2f}%) provisional subject to revision data.")
    df1 = df1[f'72019{k1}']

    ser_well0 = df1.squeeze().rename(site) * -0.3048
    ser_well0.index = ser_well0.index.tz_localize(None)

    ts_vup, coeffs_vup, unc = TSObj.calc_ts_vup(loc='well')
    # ts_vup, coeffs_vup, unc = TSObj.calc_ts_vup_dd_atm(loc='well', radius=5000, corr_thresh=0.75, show_ts=False)
    ser_vup0 = pd.Series(ts_vup, index=TSObj.dt)

    # align
    ser_well1, ser_vup1 = bbTS.align_st_en(ser_well0, ser_vup0)[0]
    ser_well2 = ser_well1.resample(freq).mean().dropna().squeeze()
    ser_vup2 = ser_vup1.resample(freq).mean().dropna()
    ser_well2, ser_vup2= bbTS.align_st_en(ser_well2, ser_vup2)[0]

    # Get shared datetime index
    shared_dt = ser_well2.index.intersection(ser_vup2.index)

    # Subset both series to those dates
    ser_well2 = ser_well2.loc[shared_dt]
    ser_vup2 = ser_vup2.loc[shared_dt]

    dct_lags = {'MS': 12*2, 'W': 52*2}
    print (f'Correlation at frequency={freq}: {np.corrcoef(ser_well2, ser_vup2)[0,1]:.2f}')
    corr, lag = bbTS.max_xcorr(ser_well2, ser_vup2, max_lags=dct_lags[freq], sign_corr=1.0)
    print (f'Max Correlation: {corr:.2f} at {lag} {freq}')
    log.warning('Max Corr may not be right yet')

    # ser_vup_p, ser_well_p = ser_vup0, ser_well0 # raw
    # ser_vup_p, ser_well_p = ser_vup1, ser_well1 # shared interval
    ser_vup_p, ser_well_p = ser_vup2, ser_well2 # shared epochs

    return ser_vup_p, ser_well_p


# almost same as above but use my preformatted file
def main1(da_well, TSObj, freq='MS', ret_shared=False):
    """
    Align and compare a groundwater well time series with a vertical uplift (vup) time series
    (e.g., from InSAR or geodetic data) at the well location. The function aligns, resamples,
    and computes the correlation between the two series, returning the aligned data and correlation.

    Parameters
    ----------
    da_well : xarray.DataArray
        Groundwater well time series for a single site. Must have dimension 'time' and attributes:
        - 'site_no': site identifier
        - 'lat': latitude of the well
        - 'lon': longitude of the well
    TSObj : object
        Object providing vertical uplift time series and metadata. Must have:
        - reg: region name
        - dt: datetime index for vup time series
        - calc_ts_vup(loc='well'): method returning (ts_vup, coeffs_vup, unc)
    freq : str, default 'MS'
        Resampling frequency for time series alignment (e.g., 'MS' for month start, 'W' for weekly).
    ret_shared : bool, default False
        If True, return the aligned time series over the shared interval (before resampling).
        If both ret_raw and ret_shared are False, return the resampled and aligned time series
        over shared epochs (default behavior).

    Returns
    -------
    ser_vup_p : pandas.Series
        Vertical uplift time series for the well location, aligned to the well series.
    ser_well_p : pandas.Series
        Groundwater well time series, aligned to the vup series.
    corr : float
        Pearson correlation coefficient between the aligned and resampled well and vup series.

    Notes
    -----
    - The function logs the correlation and maximum cross-correlation (with lag) for the site.
    - The function expects the well and vup time series to be reasonably overlapping in time.
    - The function updates the global DCT_PTS['well'] with the well's lat/lon and region.
    - Only one of ret_raw or ret_shared should be True; if both are False, returns the default
      (resampled, aligned) output.
    """
    region = TSObj.reg
    # get the well info

    site_no = da_well.site_no.item()
    ser_well0 = da_well.dropna('time').to_series().rename(site_no)

    lalo = da_well['lat'].item(), da_well['lon'].item()
    DCT_PTS['well'] = [lalo, region]


    ts_vup, coeffs_vup, unc = TSObj.calc_ts_vup(loc='well')
    # ts_vup, coeffs_vup, unc = TSObj.calc_ts_vup_dd_atm(loc='well', radius=5000, corr_thresh=0.75, show_ts=False)
    ser_vup0 = pd.Series(ts_vup, index=TSObj.dt)

    # align
    ser_well1, ser_vup1 = bbTS.align_st_en(ser_well0, ser_vup0)[0]
    ser_well2 = ser_well1.resample(freq).mean().dropna().squeeze()
    ser_vup2 = ser_vup1.resample(freq).mean().dropna()
    ser_well2, ser_vup2= bbTS.align_st_en(ser_well2, ser_vup2)[0]

    # Get shared datetime index
    shared_dt = ser_well2.index.intersection(ser_vup2.index)

    # Subset both series to those dates
    ser_well2 = ser_well2.loc[shared_dt]
    ser_vup2 = ser_vup2.loc[shared_dt]

    dct_lags = {'MS': 12*2, 'W': 52*2}
    corr = np.corrcoef(ser_well2, ser_vup2)[0,1]
    print ('')
    log.critical(f'**Site: {site_no}**')
    log.info(f'Correlation at frequency={freq}: {corr:.2f}')
    corrm, lag = bbTS.max_xcorr(ser_well2, ser_vup2, max_lags=dct_lags[freq], sign_corr=1.0)
    log.info(f'Max Correlation: {corrm:.2f} at {lag} {freq}')
    # log.warning('Max Corr may not be right yet')
    print ('')

    if ret_shared:
        ser_vup_p, ser_well_p = ser_vup1, ser_well1 # shared interval
    else:
        ser_vup_p, ser_well_p = ser_vup2, ser_well2 # shared epochs
    # log.info('Returning shared data')

    return ser_vup_p, ser_well_p, corr


if __name__ == '__main__':
    freq = 'MS'
    npix = 5
    ExpBest = ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=5)
    TSObj = PlotTS1(ExpBest.dct_exp, ExpBest.mp_exp0, ExpBest.ref_sta, neofs=ExpBest.neofs, npix=npix)

    ds_gw = xr.open_dataset(Path(os.getenv('dataroot')) / 'GW' / f'{TSObj.reg}_GW.nc')
    da_field = ds_gw['field'].dropna('time', how='all').dropna('site_no', how='all')
    da_daily = ds_gw['daily_mean'].dropna('time', how='all').dropna('site_no', how='all')
    dct_corrs1 = {}
    for da_welli in da_field.T:
        corr = main1(da_welli, TSObj, freq=freq)[-1]
        dct_corrs1[da_welli.site_no.item()] = corr

    ser1 = pd.Series(dct_corrs1)
    dct_corrs2 = {}
    for da_welli in da_daily.T:
        corr = main1(da_welli, TSObj, freq=freq)[-1]
        dct_corrs2[da_welli.site_no.item()] = corr

    ser2 = pd.Series(dct_corrs2)
    df_corrs = pd.concat([ser1, ser2], axis=1)
    df_corrs.columns = 'field daily_mean'.split()
    breakpoint()



