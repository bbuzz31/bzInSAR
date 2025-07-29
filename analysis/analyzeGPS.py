from BZ import bzBase, bbGIS, bbTS, bbPlot, bbSci
from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase

from VLM.bzGNSS import calc_MIDAS

from mintpy.utils import readfile

## first avoid jumps; last column is if its known step change in UNR
DCT_JUMPS = {
            # NYC
            'NJHT': [datetime(2020, 10, 7), 'en', False],
            'NYBR': [datetime(2020, 10, 8), 'en', True],
            'NYPR': [datetime(2015, 6, 1), 'st', False],
            'NJI2': [datetime(2013, 2, 1), 'st', True],
            # 'NJI2': [datetime(2022, 1, 1), 'en', False],
            'NYBP': [datetime(2010, 11, 1), 'st', True],

            # NNJ
            'KP14': [datetime(2019, 6, 1), 'en', False],
            'NJDY': [datetime(2012, 12, 11), 'st', True],
            'NJMT': [datetime(2012, 12, 6), 'st', True], # probably unnecessary
            'NJTP': [datetime(2017, 8, 1), 'st', True],
            'SHK5': [datetime(2011, 10, 6), 'st', True],
            'SHK6': [datetime(2007, 12, 20), 'st', True],
            'UCNJ': [datetime(2021, 3, 1), 'en', False], # maybe another


            ## HR
            'LOY2': [datetime(2018, 12, 1), 'st', False],
            # 'LOYZ': [datetime(2017, 9, 1), 'st', True], ## this is broken in midas.f
            # 'LS03': [datetime(2018, 9, 1), 'st', True], ## also broken in midas.f
            'VAHP': [datetime(2012, 3, 1), 'st', False],
            'DRV5': [datetime(2008, 3, 1), 'st', True],
            'DRV6': [datetime(2008, 3, 1), 'st', True],
            'WLP2': [datetime(2015, 1, 1), 'st', False],

            ## Miami
            # 'AOML' # ends in 2004
            # 'BOC1': [datetime(2020, 1, 29), 'st', True], # earthquake far
            # 'FCDE': [datetime(2020, 1, 29), 'st', True],
            'FCI1': [datetime(2023, 3, 1), 'en', False], # outliers, prob not needed
            # 'FCI2': [datetime(2023, 3, 1), 'en', False], # outliers, prob not needed
            'FLC5': [datetime(2017, 9, 9), 'st', True], # earthquake
            'FLC6': [datetime(2017, 9, 9), 'st', True], # earthquake
            'FLD6': [datetime(2020, 1, 29), 'st', True], # earthquake
            'FLF1': [datetime(2020, 1, 29), 'st', True], # earthquake
            'FLMB': [datetime(2020, 1, 29), 'st', True], # earthquake
            'FLN2': [datetime(2020, 1, 29), 'st', True], # earthquake
            'FTLD': [datetime(2020, 1, 29), 'st', True], # earthquake
            'HOM1': [datetime(2020, 1, 29), 'st', True], # earthquake
            'LAUD': [datetime(2020, 1, 29), 'st', True], # earthquake
            'N300': [datetime(2020, 1, 29), 'st', True], # earthquake
            'RMND': [datetime(2020, 1, 29), 'st', True], # earthquake

            ## keep these two
            # 'PBCH': [datetime(2014, 2, 1), 'st', True],
            'PBCH': [datetime(2020, 1, 29), 'st', True], # earthquake
            'ZMA1': [datetime(2007, 7, 1), 'st', True],

            ## Kiritimati
            'KRTM': [datetime(2012, 2, 18), 'st', True], # logfile
            # 'KRTM': [datetime(2012, 2, 16), 'en', True], # logfile

            ## Houston
            'NASA': [datetime(2015, 3, 1), 'st', False], # weird periodicity
            'UHDT': [datetime(2015, 3, 1), 'st', False], # weird periodicity

            # ## these two are on top of each other but very diff
            'UHL1': [datetime(2017, 10, 1), 'st', True], # earthquake
            'UHWL': [datetime(2017, 10, 1), 'st', True], # earthquake

            # # these are for dropping values outside [n] sigma
            'UHC1': [datetime(2017, 10, 1), 'st', True, 2],
            'UHC2': [datetime(2017, 10, 1), 'st', True, 2],
            'UHC3': [datetime(2017, 10, 1), 'st', True, 2],

            'ADKS': [datetime(2014, 1, 1), 'st', True], # antenna/elevation

            ### Kennedy
            'CCV5' : [datetime(2008, 3, 12), 'st', True], # receiver change
            'COKO' : [datetime(2006, 4, 15), 'st', True], # receiver change
            'TTVL' : [datetime(2017, 9, 10), 'st', True], # earthquake

            ## DC
            # 'USN8' : [datetime(2020, 8, 22), 'en', False], # Test for half records
            # 'USN8' : [datetime(2020, 8, 22), 'st', False], # Test for half records

            }

lst_outliers = 'UHC1 UHC2 UHC3'.split() # filter these by 5stds

dct_stitch   = {'Houston': 'UHDT CFJV SESG TXTG NASA CSTE DMFB TXLQ'.split()}


def get_sten(path_vup_geo):
    """ Try to get the timeperiod of insar data"""
    if op.exists(path_vup_geo):
        meta   = readfile.read_attribute(path_vup_geo)
        st, en = [meta[dd] for dd in 'START_DATE END_DATE'.split()]

    else:
        log.info ('Couldnt get InSAR data, defaulting to full S1 period (for shading)')
        today  = datetime.today()
        st, en = '20150101', f'{today.year}{str(today.month).zfill(2)}15'

    return [datetime.strptime(dt, '%Y%m%d') for dt in [st, en]]


def dl_ts(sta, ref_frame='IGS14'):
    """ Download and format a single timeseries from UNR from:

    e.g., http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/NJI2.tenv3
    """
    url0 = 'http://geodesy.unr.edu'
    base = f'{url0}/gps_timeseries/tenv3/{ref_frame}'
    url = f'{base}/{sta}.tenv3'
    try:
        df = pd.read_csv(url, sep=r'\s+')
    except Exception as e:
        raise type(e)(f'Could not download from: {url} due to: {e}') from e

    df= df.iloc[:, [0, 1, 2, 6, 8, 10, 12]+list(np.arange(14, 23, 1))]
    # these names should be changes as they're positions not vels
    df.columns = ['site', 'dates', 'decyr', 'reflon', 'vel_e', 'vel_n', 'vel_u', 'sig_e',
                    'sig_n', 'sig_u', 'corr_en', 'corr_eu', 'corr_nu', 'lat', 'lon', 'hgt']
    df = df.copy() # prevent warning
    df.dates = pd.to_datetime(df.dates, format='%y%b%d')
    return df


def dl_tss(df_gps, stas=[], ref_frame='IGS14'):
    """ Download time series at one or more stations on the fly """
    df_gps['u_vel_lsq'] = df_gps['u_vel_orig'].copy()
    df_gps['u_sig_lsq'] = df_gps['u_sig_orig'].copy()
    stas = stas if any (stas) else df_gps.sta
    stas = [stas] if isinstance(stas, str) else stas

    ## rewrite to read with pandas directly
    lst_dfs = []
    for sta in stas:
        if '/' in sta:
            log.warning('Skipping: %s', sta)
            continue

        df = dl_ts(sta, ref_frame)
        lst_dfs.append(df)

        # update the rates with a simple lsq fit as well
        # preds, coeffs = bbTS.bzSeason.fit_lsq(df.decyr, df.vel_u)
        # rate = coeffs[-2]

        p, V  = np.polyfit(df.decyr, df.vel_u, 1, cov=True)
        rate_ls, unc_ls = p[0], np.sqrt(V[0,0])

        ix   = df_gps[df_gps.sta == sta].index
        df_gps.loc[ix, 'u_vel_lsq'] = rate_ls
        df_gps.loc[ix, 'u_sig_lsq'] = unc_ls

    df_m = pd.concat(lst_dfs) if lst_dfs else None
    return df_m


def dl_coop(sta=None, ref_frame='itrf14'):
    ref_frames = 'nad83 itrf14'.split()
    assert ref_frame.lower() in ref_frames, f'Choose reference frame from: {ref_frames}'
    ref_frame = 'itrf2014' if ref_frame == 'itrf' else 'nad83_2011'
    url = f'https://noaa-cors-pds.s3.amazonaws.com/coord/coord_14/{ref_frame}_geo.comp.txt'
    df  = pd.read_csv(url, sep=r'\s+', skiprows=7, header=None)
    keep = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16]
    df   = df.iloc[:, keep]
    lats = df.loc[:, 2:5].apply(lambda x:bbGIS.dms2dd(*x), axis=1)
    lons = df.loc[:, 6:9].apply(lambda x:bbGIS.dms2dd(*x), axis=1)
    df[2] = df[2].astype(float)
    df[3] = df[3].astype(float)
    df.loc[:, 2] = lats
    df.loc[:, 3] = lons
    df.drop(np.arange(4, 10, 1), axis=1, inplace=True)

    cols = 'site lat lon height n_vel e_vel u_vel status'.split()
    df.columns = cols

    cols  = [col for col in df.columns if 'vel' in col or 'sig' in col]
    df    = df.copy()
    df[cols] /= 1000 # convert to m/y
    if sta:
        df = df[df.site == sta]
    return df


def plot_gps_rates(df_gps, d='Up', reg='?'):
        """ Scatter plot of rates across the stations """
        import geopandas as gpd
        d = d.title()
        if d == 'Up':
            dcol = 'u'
        elif d == 'East':
            dcol = 'e'
        elif d == 'North':
            dcol = 'n'

        fig, axes = plt.subplots(figsize=(10, 4))
        axes.errorbar(df_gps['sta'], df_gps[f'{dcol}_vel']*1000,
                      yerr=df_gps[f'{dcol}_sig']*1000, color='k', fmt='o')
        axes.grid(color='k', linestyle='--', alpha=0.1)
        axes.set_ylabel(f'V$_{{{d}}}$ Rate  (mm/yr)')

        mu = np.polyfit(np.ones(df_gps.shape[0]),
                        df_gps[f'{dcol}_vel'], 0,
                        w=1/df_gps[f'{dcol}_sig']).item() * 1000
        sig = df_gps['u_sig'].mean().item() * 1000

        axes.axhline(mu, linestyle='-', color='r', label=f'Weighted$_\\mu$: {mu:.2f}$\\pm$-{sig:.2f}')
        axes.axhline(mu-sig, linestyle='--', color='r', alpha=0.6)
        axes.axhline(mu+sig, linestyle='--', color='r', alpha=0.6)
        axes.legend()
        fig.set_label(f'{reg}_GPS_{d}')
        return fig, axes


def plot_df_ts_overlap(df_ts, cd='u', min_overlap=True):
    """ Overlay the timeseries of overlapping stations.
    If min_overlap compare only smallest overlapping """
    assert cd in 'e n u'.split(), f'Incorrect direction: {cd}'
    dct_lbls = {'u': 'Up', 'e': 'East', 'n': 'North'}
    for i, (site0, dfi) in enumerate(df_ts.groupby('site')):
        # get all the dates for this station
        df_tsi = df_ts[((df_ts.decyr>=dfi.decyr.min()) & (df_ts.decyr<=dfi.decyr.max()))]
        # get stations with more than 2 year overlap; do this here so min overlap works
        sites_keep = [sitei for sitei, dfi in df_tsi.groupby('site') if \
                      (dfi['dates'].max() - dfi['dates'].min()).days > (2 * 365.25)]
        df_tsi = df_tsi[df_tsi['site'].isin(sites_keep)]
        st0 = df_tsi.groupby('site').dates.min().max()
        en0 = df_tsi.groupby('site').dates.max().min()

        cmap = plt.cm.get_cmap('Set1')
        fig, axes = plt.subplots(figsize=(10, 6))

        for j, (site, dfj) in enumerate(df_tsi.groupby('site')):
            serj = dfj[f'vel_{cd} dates'.split()].set_index('dates').squeeze()
            if min_overlap:
                serj = serj[((serj.index >= st0) & (serj.index<=en0))]

            # filter and convert to mm
            dat = bbSci.nan_outliers(serj, 3)
            serj = pd.Series(dat, serj.index).dropna() * 1000
            # seri = seri[seri.abs() < seri.abs().quantile(q)].squeeze().dropna()*1000

            serj -= serj.mean()
            axes.scatter(serj.index, serj, color=cmap(j), s=5, alpha=0.8)
            ## calculate new least squares fit
            decyr = bbTS.date2dec(serj.index)
            coeffs, V = np.polyfit(decyr, serj, 1, cov=True)
            sem = np.sqrt(V[0,0])
            lsq_fit = np.polyval(coeffs, decyr)

            l, = axes.plot(serj.index, lsq_fit, c=cmap(j), linestyle='-',
                    label=f'{site}: {coeffs[0]:.2f}')#$\\pm${sem:.2f}')
            # break

        axes.legend()
        axes.set_ylabel(f'{dct_lbls[cd]} Disp. (mm)', fontsize=12)
        axes.legend()
        axes.grid(color='k', alpha=0.1, linestyle='--')
        axes.set_title(site0)
        fig.set_label(f'GPS_ts_overlap_{site0}')
    return


def plot_df_ts_dates(df_ts0, df_rates=None, cd='u', st=None, en=None, min_yrs=2, show_trend=True):
    """ Overlay the timeseries during certain dates with a min_yrs overlap """
    assert cd in 'e n u'.split(), f'Incorrect direction: {cd}'
    st = pd.to_datetime(st) if isinstance(st, str) else st
    en = pd.to_datetime(en) if isinstance(en, str) else en
    stt, enn = st.strftime('%b, %Y'), en.strftime('%b, %Y')

    dct_lbls = {'u': 'Up', 'e': 'East', 'n': 'North'}
    df_ts1 = df_ts0[((df_ts0.dates>=st) & (df_ts0.dates<=en))]
    sites_keep = [sitei for sitei, dfi in df_ts1.groupby('site') if \
                  (dfi['dates'].max() - dfi['dates'].min()).days > (min_yrs * 365.25)]
    df_ts1 = df_ts1[df_ts1['site'].isin(sites_keep)]
    print (f'Keeping {len(sites_keep)} of {len(df_ts0.site.unique())} GPS stations that overlap {stt} to {enn}.')

    cmap = plt.cm.get_cmap('Set1')
    fig, axes = plt.subplots(figsize=(10, 5))

    for j, (site, dfj) in enumerate(df_ts1.groupby('site')):
        serj = dfj[f'vel_{cd} dates'.split()].set_index('dates').squeeze()

        # filter and convert to mm
        dat = bbSci.nan_outliers(serj, 3)
        serj = pd.Series(dat, serj.index).dropna() * 1000
        # seri = seri[seri.abs() < seri.abs().quantile(q)].squeeze().dropna()*1000

        serj -= serj.mean()
        parms = {'alpha': 0.1} if show_trend else {'label': site}
        axes.scatter(serj.index, serj, color=cmap(j), s=5, **parms)
        ## calculate new least squares fit
        decyr = bbTS.date2dec(serj.index)
        coeffs, V = np.polyfit(decyr, serj, 1, cov=True)
        sem = np.sqrt(V[0,0])
        # hack to get the full midas uncertaitny for now
        if df_rates is not None:
            sem = df_rates[df_rates['sta'] == site][f'{cd}_sig'].item()*1000

        lsq_fit = np.polyval(coeffs, decyr)

        if show_trend:
            l, = axes.plot(serj.index, lsq_fit, c=cmap(j), linestyle='-',
                label=f'{site}: {coeffs[0]:.2f}$\\pm${sem:.2f}')
        # break

    axes.legend()
    axes.set_ylabel(f'{dct_lbls[cd]} Disp. (mm)', fontsize=12)
    axes.set_title(f'{stt} to {enn}', fontsize=15)
    axes.legend()
    axes.grid(color='k', alpha=0.1, linestyle='--')
    fig.set_label(f'GPS_ts_overlap_{st.date()}_{en.date()}_{dct_lbls[cd]}')
    return


class GPS_TS_Analysis(bzBase):
    def __init__(self, Exp, path_gps=None):
        """ Recompute the MIDAS rates and update 'path_gps'

        - Exp: an initialized Experiment object
        - path_gps: optional (default use Exp.path_gps, which is midas cln)

        Optionally plot timeseries of GPS and compare with InSAR
        Optionally set a lst_force [sta, rate, unc] to overwrite value. (Kiribati TG-ALT)
        """
        super().__init__()
        self.Exp      = Exp
        self.path_gps = Exp.path_gps if path_gps is None else path_gps
        self.df_gps = prep_gps(self.path_gps, Exp.reg, 'm', True) # keep horiz
        self.ins_st, self.ins_en = get_sten(self.Exp.path_vup_geo)

        # in case already updated
        if not 'u_vel_orig' in self.df_gps.columns:
            self.df_gps['u_vel_orig'] = self.df_gps['u_vel'].to_numpy()
            self.df_gps['u_sig_orig'] = self.df_gps['u_sig'].to_numpy()

            self.df_gps['ts_start_orig'] = self.df_gps['ts_start']
            self.df_gps['ts_end_orig']   = self.df_gps['ts_end']

        self.stas   = self.df_gps.sta


    def update_MIDAS(self, tthresh=0.4, start=None, end_insar=True,
                     show=True, write=True, lst_force=False):
        """ Overwrite the u_vel / u_sig column with a new MIDAS calculation

        tthresh is the amount (percent of days) of overlap with insar
        start = force it to start after a certain date
        end_insar = stop the GPS rate computation at same time as InSAR
        lst_force can manually override with a value (somehow)
        """
        # log.setLevel(10)
        df_ts = dl_tss(self.df_gps, self.stas).reset_index()
        if show:
            log.info('Plotting updated MIDAS. Compare with original on website '\
                     'and/or run the GPS_Rate_Analysis.plot_trends_new')
            # self.plot_ts(df_ts, dct_stitch[self.Exp.reg])
            self.plot_ts(df_ts)

        ## Use this to inspect the time series and find the jumps
        # self._find_jump(df_ts, 'USN7')

        ## new midas trends for the edited time series
        stas = self.df_gps.sta.tolist()
        for sta in stas:
            df_ts_sub = df_ts[df_ts.site == sta]
            lst = DCT_JUMPS.get(sta, [])
            if lst:
                ## cut by the jump
                if lst[1] == 'en':
                    st, en = df_ts_sub.dates.iloc[0], lst[0]
                    bad_ix = df_ts_sub[df_ts_sub.dates >= en].index
                    n_sig  = None
                    df_ts.drop(bad_ix, inplace=True)

                elif lst[1] == 'st':
                    st, en = lst[0], df_ts_sub.dates.iloc[-1]
                    bad_ix = df_ts_sub[df_ts_sub.dates <= st].index
                    n_sig  = None
                    df_ts.drop(bad_ix, inplace=True)

                # get number of sigma to cut outliers by
                if len(lst) == 2:
                    st, en = df_ts_sub.dates.iloc[0], df_ts_sub.dates.iloc[-1]
                    n_sig  = lst[1]

                elif len(lst) == 4:
                    n_sig  = lst[3]


            if start is not None:
                start = pd.to_datetime(start) if isinstance(start, str) else start
                df_ts_sub1 = df_ts_sub[df_ts_sub['dates'] >= start]
                if df_ts_sub1.empty:
                    log.warning(f'{sta} has no data after {start}, skipping')
                else:
                    df_ts_sub = df_ts_sub1
                st = df_ts_sub.dates.iloc[0]
                en = df_ts_sub.dates.iloc[-1]
                n_sig = None

            if end_insar:
                df_ts_sub = df_ts_sub[df_ts_sub['dates'] <= self.ins_en]
                st = df_ts_sub.dates.iloc[0]
                en = df_ts_sub.dates.iloc[-1]
                n_sig = None

            else:
                continue

            if df_ts_sub['decyr'].max() - df_ts_sub['decyr'].min() < 1.02:
                log.warning(f'{sta} timeseries is too short for MIDAS, skipping')
                continue
            st = df_ts_sub.dates.iloc[0]
            en = df_ts_sub.dates.iloc[-1]
            n_sig = None

            ## recompute midas
            rate, unc = calc_MIDAS.midas_fortran_one(sta, st, en, n_sig, plot=show)

            ## compute LSQ on new timeseries
            p, V  = np.polyfit(df_ts_sub['decyr'], df_ts_sub['vel_u'], 1, cov=True)
            rate_ls, unc_ls = p[0], np.sqrt(V[0,0])

            ix = self.df_gps[self.df_gps.sta == sta].index
            self.df_gps.loc[ix, 'u_vel'] = rate
            self.df_gps.loc[ix, 'u_sig'] = unc

            self.df_gps.loc[ix, 'u_vel_lsq'] = rate_ls
            self.df_gps.loc[ix, 'u_sig_lsq'] = unc_ls

            self.df_gps.loc[ix, 'ts_start'] = st.date()
            self.df_gps.loc[ix, 'ts_end']   = en.date()

        ## take only the stations that have some overlap with InSAR
        if tthresh>0:
            stas_good   = self.temporal_overlap(df_ts, thresh=tthresh)
            self.df_gps = self.df_gps[self.df_gps['sta'].isin(stas_good)]
        else:
            stas_good = self.df_gps.sta.tolist()

        ## also add on the NOAA stations for a double check
        df_coop = dl_coop(ref_frame='ITRF14')
        df_coop = df_coop[df_coop.site.isin(stas_good)].rename(
            columns={'u_vel': 'u_vel_coop', 'e_vel':'e_vel_coop', 'n_vel': 'n_vel_coop'})

        df_coop1 = pd.concat([df_coop.site, df_coop.filter(like='vel')], axis=1)
        self.df_gps = pd.merge(self.df_gps, df_coop1, left_on='sta',
                               right_on='site', how='outer')

        # for showing the number comparisons
        df_uvel = pd.concat([self.df_gps['sta'], self.df_gps.filter(like='u_vel')*1000], axis=1)

        if lst_force:
            log.critical('Manually overwriting gps u_vel/sig')
            f_sta, f_ref, f_unc = lst_force
            df_gps = self.df_gps.set_index('sta')
            df_gps.loc[f_sta, 'u_vel'] = f_ref
            df_gps.loc[f_sta, 'u_sig'] = f_unc
            self.df_gps = df_gps.reset_index()

        if show:
            print (df_uvel)

        ## updated MIDAS with new midas rates (u_vel1) and hand calculated (u_vel2)
        if write:
            self.df_gps.to_csv(self.path_gps, index=False)
            log.info('Wrote: %s', self.path_gps)
        return self.df_gps


    def plot_ts(self, df_ts, vel='up', stas=[], max_nrows=7, ylims=[-50, 50]):
        """ Plot all (or named subset) of the GPS time series """
        vel = vel.title()
        assert vel in 'Up East North'.split(), f'incorrect vel: {vel}'
        col      = 'site' if 'site' in df_ts.columns else 'sta'
        stas     = df_ts[col].unique() if not stas else stas
        lst_sta  = np.arange(len(stas))
        num_sets = (len(stas) - 1) // max_nrows + 1

        # convert to mm
        cols  = [col for col in df_ts.columns if 'vel' in col or 'sig' in col]
        df_ts = df_ts.copy()
        df_ts[cols] *= 1000

        if vel == 'Up':
            dcol = 'vel_u'
        elif vel == 'East':
            dcol = 'vel_e'
        elif vel == 'North':
            dcol = 'vel_n'

        for num_set in range(num_sets):
            try:
                stas_i = stas[num_set*max_nrows:(num_set*max_nrows)+max_nrows]
            except:
                stas_i = stas[num_set*max_nrows:-1]

            fig, axes = plt.subplots(figsize=(10, 6), nrows=min(max_nrows, len(stas_i)),
                                     sharex=True)
            axes = axes if isinstance(axes, Iterable) else [axes]
            axes[0].set_title('Original MIDAS')

            i = 0
            for sta, df in df_ts.groupby(col):
                if not sta in stas_i:
                    continue

                ## truncate to sort of recent ( after 2010)
                # df = df[df.dates>='20100101']

                dat = df[dcol] - df[dcol].mean()

                axes[i].scatter(df.dates, dat, c='k', alpha=0.5, s=2)
                # plot the original midas trend
                gps_midas = self.df_gps[self.df_gps.sta==sta]
                if vel == 'Up':
                    velp, velstdp = gps_midas.u_vel_orig.item()*1000, gps_midas.u_sig_orig.item()*1000
                else:
                    dcol1 = '_'.join(dcol.split('_')[::-1])
                    velp, velstdp = gps_midas[dcol1].item()*1000, gps_midas[dcol1].item()*1000
                midas_trend = np.polyval([velp, dat.iloc[0]], df.decyr)
                midas_trend -= midas_trend.mean()

                l, = axes[i].plot(df.dates, midas_trend, c='r', linestyle='--',
                        label=f'{velp:.2f}$\\pm${velstdp:.2f} (mm/yr)')
                axes[i].fill_between(df.dates, midas_trend - velstdp,
                                midas_trend+velstdp, color=l.get_color(), alpha=0.3)

                axes[i].legend()
                ax2 = axes[i].twinx()
                ax2.yaxis.set_ticks([])
                if len(stas) < 5:
                    axes[i].set_title(sta, fontsize=13)
                else:
                    ax2.set_ylabel(sta, rotation=270, labelpad=10, fontsize=14)

                # plot the insar temporal period
                axes[i].axvspan(self.ins_st, self.ins_en, alpha=0.3, color='darkgray', label='InSAR Period')

                axes[i].grid(color='gray', linestyle = '--', linewidth=0.1)

                axes[i].set_ylim(ylims)
                i+=1

            mdpt = int(len(axes)/2)
            axes[mdpt].set_ylabel(f'V$_{{{vel}}}$ (mm)', fontsize=14)
            fig.subplots_adjust(hspace=0.2)
            fig.set_label(f'{self.Exp.reg}_GPS_TS_raw{num_set}_{vel}')
        log.info('Showing Original MIDAS rates')


    def temporal_overlap(self, df_ts, thresh=0.4):
        """ Compute the percent of days that GPS time series occurs
        within days spanned by InSAR

        Takes all time series from dl_tss and returns sta if its good
        """
        # number of days covered by InSAR
        insar_days = (self.ins_en - self.ins_st).days
        log.debug ('InSAR spans: %s to %s', self.ins_st.date(), self.ins_en.date())

        stas_good = []
        for sta in df_ts.site.unique():
            df_ts1  = df_ts[df_ts.site == sta]
            ts_st, ts_en = df_ts1.dates.iloc[0], df_ts1.dates.iloc[-1]
            log.info('%s spans: %s to %s', sta, ts_st.date(), ts_en.date())

            # get overlapping gps with insar
            df_ts_ins = df_ts1[((df_ts1.dates>=self.ins_st) & (df_ts1.dates<=self.ins_en))]

            # number of days with GPS measurements
            n_days_gps = df_ts_ins.set_index('dates').drop('site',
                                axis=1).resample('D').mean().dropna().shape[0]


            per_cov  = (n_days_gps / insar_days)*100
            gps_good = per_cov >= (thresh*100)
            log.debug (f'{sta} covers {per_cov:.1f}% of days (keep threshold = {thresh*100})')
            stas_good.append(sta) if gps_good else ''

        sta_bad = [sta for sta in df_ts.site.unique() if not sta in stas_good]

        log.warning (f'{", ".join(sta_bad)} do not meet threshold of {thresh*100}%.')

        return stas_good


    def _find_jump(self, df_ts, sta):
        """ Manually find where the jump location is """
        df = df_ts[df_ts.site == sta]
        ser = df.set_index('dates')['vel_u'].rename(sta)
        # ser.diff().plot()
        # moving trends
        # ser = make_all_ste_en(ser.index)

        # calc_moving_trends_1d(ser, nyrs=2)
        # calc_trend_evolve(ser, 2)
        # jmp = datetime(2020, 10, 5)
        # ser1 = ser[((ser.index>='20201001') & (ser.index<='20201015'))]
        fig, axes = plt.subplots(figsize=(10,4))
        axes.plot(ser.index, ser)
        plt.show()

        breakpoint()


    def plot_nearby_stations(self, buff=100, sta=None, lalo=None, d='Up'):
        import geopandas as gpd
        from shapely.geometry import Point
        assert not (sta is None and lalo is None), 'Use either gps sta or lalo'
        assert  sta is not None and not lalo is not None, 'Use eiher gps sta or lalo'
        gdf_gps = bbGIS.df2gdf(self.df_gps, 'all').to_crs(4087)
        if sta is not None:
            assert sta in self.df_gps.sta.tolist(), f'Incorrect station name: {sta}'
            gser_sta = gdf_gps[gdf_gps['sta'] == sta ]
            lalo = gser_sta.geometry.y.item(), gser_sta.geometry.x.item()
        else:
            sta = f'{lalo[0]:.3f},{lalo[1]:.3f}'


        # get distance of all stations
        gser_pt = gpd.GeoSeries(Point(lalo[::-1]), crs=gdf_gps.crs)
        gdf_gps[f'dist_{sta}'] = gdf_gps.distance(Point(lalo[::-1]))

        gdf_gps_close = gdf_gps[gdf_gps[f'dist_{sta}'] <= buff]

        d = d.title()
        if d == 'Up':
            dcol = 'u'
        elif d == 'East':
            dcol = 'e'
        elif d == 'North':
            dcol = 'n'

        fig, axes = plt.subplots(figsize=(10, 4))
        axes.errorbar(gdf_gps_close['sta'], gdf_gps_close[f'{dcol}_vel']*1000,
                      yerr=gdf_gps_close[f'{dcol}_sig']*1000, color='k', fmt='o')
        axes.grid(color='k', linestyle='--', alpha=0.1)
        axes.set_ylabel(f'V$_{{{d}}}$ Rate  (mm/yr)')
        axes.set_title(f'{gdf_gps_close.shape[0]} stations within {buff} m of {sta}')

        mu = np.polyfit(np.ones(gdf_gps_close.shape[0]),
                        gdf_gps_close[f'{dcol}_vel'], 0,
                        w=1/gdf_gps_close[f'{dcol}_sig']).item() * 1000
        sig = gdf_gps_close['u_sig'].mean().item() * 1000

        axes.axhline(mu, linestyle='-', color='r', label=f'Weighted$_\\mu$: {mu:.2f}$\\pm$-{sig:.2f}')
        axes.axhline(mu-sig, linestyle='--', color='r', alpha=0.6)
        axes.axhline(mu+sig, linestyle='--', color='r', alpha=0.6)
        axes.legend()
        fig.set_label(f'{self.Exp.reg}_GPS_{sta}_{buff}_{d}')
        return fig, axes, gdf_gps_close


    def plot_map(self):
        """ Wrapper around plotVLM"""
        from BZ.VLM.bzFRInGE.plotting.plotVLM import PlotVLM
        PO = PlotVLM(self.Exp.dct_exp)
        PO.plot_GPS()
        PO.plot_interp_GPS()
        return


class GPS_Rate_Analysis(bzBase):
    """ Analyze the rates just updated with GPS_TS_Analysis """
    def __init__(self, Exp, path_gps=None):
        super().__init__()
        self.Exp      = Exp
        self.path_gps = self.Exp.path_gps if path_gps is None else path_gps
        self.df_gps = prep_gps(self.path_gps, self.Exp.reg, units='mm', exclude=True)
        self.ins_st, self.ins_en = get_sten(self.Exp.path_vup_geo)
        # self.plot_trends_new()
        # self.plot_insar_resid()


    def plot_insar_resid(self, npix=0):
        df_gps = self.df_gps.copy()
        # use the new midas rates
        df_gps.u_vel = df_gps.u_vel/1000
        df_gps.u_sig = df_gps.u_sig/1000

        ref_stas = self.Exp.get_ref_stas()
        # manually hack for HR at the moment (show comparisons for offset method)
        if 'Stitch2020' in self.Exp.mp_exp and self.Exp.reg == 'HR':
            ref_stas = ['VAHP']

        if npix > 0:
            ref_stas = []

        df_resid = gps_rmse_npix(self.Exp, df_gps=df_gps, npix=npix, exclude=ref_stas)
        ser_resid = df_resid['resid']
        ser_resid *=1000

        ## plot the RMSE and actual resid for the two comparisons
        fig, axes = plt.subplots(figsize=(10,4))
        mu = ser_resid.abs().mean()
        axes.scatter(ser_resid.index, ser_resid, color='k')
        axes.plot(ser_resid.index, np.tile(mu, len(ser_resid.index)), color='k',
                  label=f'|RMSE|: {mu:.1f} (mm/yr)', linestyle='--')

        axes.grid(color='gray', linestyle = '--', linewidth=0.1)
        axes.legend(ncol=1)
        axes.set_ylabel('Residual$_{GPS-Ins}$ (mm/yr)')
        axes.set_xlabel('GNSS Station')
        fig.set_label(f'{self.Exp.mp_exp0}_gps_v_insar_rmse')


    def plot_gps_vs_insar(self, npix=0, unc_lsq=False, use_nc=False):
        """ Plot actual GPS vs InSAR

        Can use GPS uncertainties calculated via lsq to compare with InSAR
        Optionally use an different datarray (uses diff function)
        """
        df_gps = self.df_gps.copy()
        # use the new midas rates
        unc_col = 'u_sig' if not unc_lsq else 'u_sig_lsq'
        df_gps.u_vel = df_gps['u_vel']/1000
        df_gps.u_sig = df_gps[unc_col]/1000

        ref_stas = self.Exp.get_ref_stas()
        # manually hack for HR at the moment (show comparisons for offset method)
        if 'Stitch2020' in self.Exp.mp_exp and self.Exp.reg == 'HR':
            ref_stas = ref_stas[:2]

        ## maybe could include reference stations if using extra pixels
        # if npix > 0:
        #     ref_stas = []
        nn = df_gps.shape[0] - len(ref_stas)

        if not use_nc:
            df1 = gps_rmse_npix(
                self.Exp, df_gps=df_gps, npix=npix, exclude=ref_stas)
        else:
            df1 = gps_rmse_npix_nc(self.Exp, df_gps=df_gps, npix=npix, exclude=ref_stas)

        ser_resid = df1['resid']

        log.info(f'Considering {ser_resid.shape[0]} stations of {nn} non-reference stations')

        ser_resid *= 1000
        df1 *= 1000

        ## compare using least squares GPS instead of MIDAS
        # too small
        # if self.Exp.reg == 'NYC' and unc_lsq:
        #     ref_sta_unc = self.df_gps.set_index('sta').loc[ref_stas, 'u_sig'].mean()
        #     df1['ins_unc_only'] = np.sqrt(df1['ins_unc']**2 - ref_sta_unc**2)
        #     ref_sta_unc_lsq     = self.df_gps.set_index('sta').loc[ref_stas, 'u_sig_lsq'].mean()
        #     df1['ins_unc'] = np.sqrt(df1['ins_unc_only']**2 + ref_sta_unc_lsq**2)

        ## plot the RMSE and actual resid for the two comparisons
        fig, axes = plt.subplots(figsize=(10,4))
        mu  = ser_resid.abs().mean().item()

        x = np.arange(len(df1.index))
        s = axes.errorbar(x-0.05, df1['gps_vel'], yerr=df1['gps_unc'], color='k', label='GNSS', fmt='o')

        # len(categories)) + scatter_offset
        axes.plot(df1.index, df1['ins_vel'], color='k', label=f'RMSE: {mu:.1f} (mm/yr)', alpha=0)
        axes.errorbar(x+0.05, df1['ins_vel'], yerr=df1['ins_unc'], color='red', label='InSAR', fmt='o')

        # axes.plot(ser_resid.index, np.tile(mu, len(ser_resid.index)), color='k',
        #           label=f'|RMSE|: {mu:.1f} (mm/yr)', linestyle='--')

        axes.grid(color='gray', linestyle='--', alpha=0.1)
        loc = 'lower right' if self.Exp == 'NYC' else 'best'
        axes.legend(ncol=1, loc='lower right')
        axes.set_ylabel('Vertical Rate (mm/yr)')
        axes.set_xlabel('GNSS Station')
        if len(x) > 20:
            axes.set_xticklabels(axes.get_xticklabels(), rotation=45)
        # axes.axhline(0, color='gray', linestyle='--', alpha=0.3)
        dct_ylims = {'NYC': [-4, 1], 'HR': [-8, 1]}
        axes.set_ylim(dct_ylims.get(self.Exp.reg, (None, None))) # NYC
        fig.set_label(f'{self.Exp.mp_exp0}_gps_v_insar')
        return fig, axes


    def plot_gps_vs_insar_qq(self, npix=0, fs=15, llim=-10, ulim=10, axes=None, use_nc=False):
        """ Plot actual GPS vs InSAR on x / y axis """
        ref_stas = self.Exp.get_ref_stas()
        # manually hack for HR at the moment (show comparisons for offset method)
        if 'Stitch2020' in self.Exp.mp_exp and self.Exp.reg == 'HR':
            ref_stas = ref_stas[:2]

        ## maybe could include reference stations if using extra pixels
        # if npix > 0:
        #     ref_stas = []
        nn = self.df_gps.shape[0] - len(ref_stas)

        if not use_nc:
            df_cmp = gps_rmse_npix(
                self.Exp, npix=npix, exclude=ref_stas)
        else:
            df_cmp = gps_rmse_npix_nc(self.Exp, npix=npix, exclude=ref_stas)
        df_cmp  *= 1000
        ser_resid = df_cmp['resid']

        log.info(f'Considering {ser_resid.shape[0]} stations of {nn} non-reference stations')
        rmse = ser_resid.abs().sum() / (len(ser_resid)-1)
        r2   = calc_r2(df_cmp['ins_vel'], df_cmp['gps_vel'])
        if axes is None:
            fig, axes = plt.subplots(figsize=(10,6))
        else:
            fig = axes.get_figure()

        # axes.errorbar(df_cmp['gps_vel'], df_cmp['ins_vel'], yerr=df_cmp['ins_unc'], xerr=df_cmp['gps_unc'])
        axes.scatter(df_cmp['ins_vel'], df_cmp['gps_vel'], color='k')#

        axes.set_xlim([llim, ulim])
        axes.set_ylim([llim, ulim])
        tickn = 1 if ulim - llim < 11 else 2
        axes.set_xticks(np.arange(llim, ulim+tickn, tickn))
        axes.set_yticks(np.arange(llim, ulim+tickn, tickn))

        line = plt.Line2D([0, 1], [0, 1], color='red', linestyle='--', transform=axes.transAxes)
        axes.add_line(line)

        # add stats
        at = axes.transAxes
        axes.text(0.04, 0.9, f'RMSE: {rmse:.2f} (mm/yr)', fontsize=fs, transform=at)
        axes.text(0.04, 0.83, f'R$^2$: {r2:.2f} (mm/yr)', fontsize=fs, transform=at)

        axes.grid(color='k', linestyle='--', alpha=0.1)
        axes.set_ylabel('GPS (mm/yr)', fontsize=fs)
        axes.set_xlabel('InSAR (mm/yr)', fontsize=fs)
        fig.set_label('GPS_InSAR_scatter')
        axes.tick_params(axis='both', which='major', labelsize=fs-2)
        return fig, axes


    def plot_trends_new(self, max_rows=7):
        """ Plot the GPS stations trend before & after update """
        import matplotlib.dates as mdates
        df_gps = self.df_gps[self.df_gps.u_vel != self.df_gps.u_vel_orig]
        if df_gps.empty:
            print ('New trends == old trends, nothing to show.')
            return
        self.log.info(f'Got: {df_gps.shape[0]} different stations')

        stas = df_gps.sta.tolist()
        num_sets = (len(df_gps) - 1) // max_rows + 1

        for num_set in range(num_sets):
            try:
                stas_i = stas[num_set*max_rows:(num_set*max_rows)+max_rows]
            except:
                stas_i = stas[num_set*max_rows:-1]

            fig, axes = plt.subplots(figsize=(10, 6), ncols=2,sharex=True,
                                     sharey='row', nrows=min(max_rows, len(stas_i)))

            axes = axes if isinstance(axes, Iterable) else np.array([axes])
            axes = axes[np.newaxis, ...] if axes.ndim == 1 else axes # in case 1 station

            # print (df_gps)
            for i, sta in enumerate(stas_i):
                df = df_gps.set_index('sta').loc[sta]
                for j, kind in enumerate(['_orig', '']):
                    dates = pd.date_range(df[f'ts_start{kind}'], df[f'ts_end{kind}'])
                    if sta == 'NJI2' and j == 0: # cuz it jumps
                        dates = pd.date_range('20020101', df[f'ts_end{kind}'])

                    df_ts = dl_ts(sta)
                    df_ts = df_ts[df_ts.dates>=dates[0]]
                    df_ts = df_ts[df_ts.dates<=dates[-1]]

                    # get and cut outliers if exists for this station in steps
                    sta_info = DCT_JUMPS.get(sta, [])

                    if not sta_info:
                        df_ts = df_ts[df_ts.dates<=self.ins_en]
                        sta_info = [self.ins_en] # put the date in for plotting cyan line
                        color = 'cyan'

                    elif len(sta_info) in [2, 4]:
                        n_sig, vc = sta_info[len(sta_info)-1], 'vel_u'
                        mu, sig = df_ts[vc].mean(), df_ts[vc].std()
                        mask    = (df_ts[vc] < mu + n_sig*sig) & (df_ts[vc] > mu - n_sig*sig)
                        df_ts   = df_ts[mask]


                    ## cut out outliers for visualization
                    # df_ts = df_ts[df_ts.vel_u.abs() <= df_ts.vel_u.abs().quantile(0.95)]
                    ts_vup = 1000 * (df_ts.vel_u  - df_ts.vel_u.mean())

                    axes[i, j].scatter(df_ts.dates, ts_vup, color='k', s=2, alpha=0.7)

                    decyr = bbTS.date2dec(dates)
                    vel, sig = df[f'u_vel{kind}'], df[f'u_sig{kind}']
                    trend  = np.polyval([vel, ts_vup.iloc[0]], decyr)
                    trend -= trend.mean()
                    l, = axes[i, j].plot(dates, trend, c='r',
                                        label=f'{vel:.1f} $\\pm$ {sig:.1f} mm/yr')
                    axes[i, j].fill_between(dates, trend-sig, trend+sig,
                                            color=l.get_color(), alpha=0.4)
                    axes[i, j].legend(loc='upper right')

                    # plot the jumps
                    # green for start jump, cyan for end jump
                    if len(sta_info) in [3, 4]:
                        color = 'cyan' if sta_info[2] else 'chartreuse'

                    axes[i,j].axvline(sta_info[0], color=color, linestyle='--')

                    if j == 1: # only shade the right
                        axes[i,j].axvspan(self.ins_st, self.ins_en, alpha=0.3, color='darkgray')#, label='InSAR Period')
                        for spine in axes[i,j].spines.values():
                            spine.set_edgecolor('darkblue')

                    if self.Exp.reg in 'NYC NNJ'.split():
                        axes[i, j].set_ylim([-25, 25])
                        st1 = datetime(2012, 1, 1)
                        en1 = datetime(2023, 6, 30)
                        axes[i, j].set_xlim([st1, en1])

                    elif self.Exp.reg == 'HR':
                        axes[i, j].set_ylim([-50, 50]) # HR


                # axes[i,0].set_ylabel(f'{sta} Vup (mm)')
                axes[i,0].set_ylabel(f'{sta}')

            axes[0, 0].set_title('Original MIDAS')
            axes[0, 1].set_title('New MIDAS')
            [ax.grid(color='gray', linestyle = '--', linewidth=0.1) for ax in axes.ravel()]

            fig.subplots_adjust(hspace=0.3, wspace=0.02)
            fig.set_label(f'{self.Exp.reg}_GPS_TS_compare{i}')
        return


    def plot_trends_new_scatter(self, llim=-10, ulim=10):
        """ Plot the different rates vs one another, using unc of shorter (new) """
        df_gps = self.df_gps[self.df_gps.u_vel != self.df_gps.u_vel_orig]
        if df_gps.empty:
            print ('New trends == old trends, nothing to show.')
            return
        self.log.info(f'Got: {df_gps.shape[0]} different stations')

        fig, axes = plt.subplots(figsize=(10,6))
        axes.errorbar(df_gps['u_vel_orig'], df_gps['u_vel'], yerr=df_gps['u_sig'], color='k', fmt='o')
        axes.set_xlabel('Original Trend (mm/yr)', fontsize=12)
        axes.set_ylabel('Updated Trend (InSAR Time Period; mm/yr)', fontsize=12)

        tickn = 1 if ulim - llim < 11 else 2
        axes.set_xticks(np.arange(llim, ulim+tickn, tickn))
        axes.set_yticks(np.arange(llim, ulim+tickn, tickn))

        line = plt.Line2D([0, 1], [0, 1], color='red', linestyle='--', transform=axes.transAxes)
        axes.add_line(line)
        axes.grid(color='k', linestyle='--', alpha=0.1)
        return


    def plot_trends_diff(self):
        """ Plot the difference between the MIDAS full and truncated trend """
        df_gps = self.df_gps[self.df_gps.u_vel != self.df_gps.u_vel_orig]
        self.df_gps['resid'] = self.df_gps['u_vel'] - self.df_gps['u_vel_orig']
        bounds = [-3, -2, -1, 1, 2, 3]
        norm = mpl.colors.BoundaryNorm(bounds, 256)
        f, a = bbPlot.plot_dfll(self.df_gps, c=self.df_gps['resid'], use_dark=True,
                         adj_fact=1e-3, s=200, zoom=13, norm=norm)
        im  = a.get_children()[0]
        im.colorbar.ax.set_ylabel(f'$\\Delta$V$_{{Up}}$ (mm/yr)')
        return


### -------- misc functions for looking at time series evolution -----------

def make_all_st_en(lst_dates, nyrs=2, overlap=True):
    """ Get all possible nyr chunks between st and en """
    from dateutil.relativedelta import relativedelta
    delyrs = relativedelta(years=nyrs)

    lst_sten = []
    for i, dt in enumerate(lst_dates):
        st = dt
        ## for independent segments
        if i > 0 and not overlap and st < en:
            continue

        en = st+delyrs
        if en > lst_dates[-1]:
            break
        lst_sten.append((st, en))

    return lst_sten


def make_all_st_en_evolve(lst_dates, styrs=2, freq='m'):
    """ make a list of start end points for a timeseries that starts with 2 yrs
    and then adds a month over and over again
    """
    from dateutil.relativedelta import relativedelta
    delyrs = relativedelta(years=styrs)

    st = lst_dates[0]
    en = st + delyrs
    i  = 0

    lst_sten = []
    while en <= lst_dates[-1]:
        en  = en + relativedelta(months=i)
        # grab the actual date thats closest to (but before) the new end date
        nearest = lst_dates[lst_dates<=en][-1]
        lst_sten.append((st, nearest))
        i += 1

    return lst_sten


def calc_moving_trends_1d(ser, nyrs, overlap=True, plot=True):
    """ Make a list of nyr trends from a PANDAS SERIES with dates on index """
    lst_rates = []
    lst_sten  = make_all_st_en(ser.index, nyrs, overlap)

    for i, sten in enumerate(lst_sten):
        ser1 = ser[((ser.index>=sten[0]) & (ser.index<=sten[1]))]
        dt   = bbTS.date2dec(ser1.index)
        rate = np.polyfit(dt, ser1, 1)[0]
        lst_rates.append(rate)

    st = [lst_sten[i][0] for i in range(len(lst_sten))]
    en = [lst_sten[i][1] for i in range(len(lst_sten))]
    df = pd.DataFrame({'rate':lst_rates, 'st': st, 'en':en})

    df.rate *= 100 # convert to cm
    if plot:
        fig, axes = plt.subplots(figsize=(10,4))
        axes.plot(df.st, df.rate)
        axes.set_xlabel('Start Date')
        axes.grid(color='gray', linestyle = '--', linewidth=0.1)
        ax2 = axes.twiny()

        # Move twinned axis ticks and label from top to bottom
        ax2.xaxis.set_ticks_position("bottom")
        ax2.xaxis.set_label_position("bottom")

        # Offset the twin axis below the host
        ax2.spines["bottom"].set_position(("axes", -0.5))

        # Turn on the frame for the twin axis, but then hide all
        # but the bottom spine
        ax2.set_frame_on(True)
        ax2.patch.set_visible(False)

        for sp in ax2.spines.values():
            sp.set_visible(False)
        ax2.spines["bottom"].set_visible(True)

        ax2.set_xlabel('End Date')

        ax2.plot(df.en, df.rate)

        axes.set_ylabel(f'{nyrs} Year Trends (cm/yr)')

        fig.subplots_adjust(bottom=0.2)

    return


def calc_trend_evolve(ser, nyrs_st=2, plot=True):
    """ Calculate how the trend evolves in time starting with 2 yrs """
    lst_sten = make_all_st_en_evolve(ser.index, nyrs_st)

    lst_rates = []
    for i, sten in enumerate(lst_sten):
        ser1 = ser[((ser.index>=sten[0]) & (ser.index<=sten[1]))]
        dt   = bbTS.date2dec(ser1.index)
        rate = np.polyfit(dt, ser1, 1)[0]
        lst_rates.append(rate)

    st = [lst_sten[i][0] for i in range(len(lst_sten))]
    en = [lst_sten[i][1] for i in range(len(lst_sten))]
    df = pd.DataFrame({'rate':lst_rates, 'st': st, 'en':en})

    df.rate *= 1000 # convert to cm
    if plot:
        fig, axes = plt.subplots(figsize=(10,4))
        axes.plot(df.en, df.rate)
        axes.set_title(f'Trend Start: {df.st.iloc[0].date()}')
        axes.set_ylabel('Vertical Rate (mm/yr)')
        axes.set_xlabel('Trend End')
        axes.grid(color='gray', linestyle = '--', linewidth=0.1)


def analyze_one_station(sta='NRL1', st1=None, en1='20130101'):
    """ Carefully check out one station and update MIDAS """
    df0 = dl_ts(sta).sort_values('dates')
    df0['vel_u'] *= 1000
    df0['sig_u'] *= 1000
    sten0 = (df0.dates.iloc[0], df0.dates.iloc[-1])
    st1 = df0.dates.iloc[0] if st1 is None else pd.to_datetime(st1)
    en1 = df0.dates.iloc[-1] if en1 is None else pd.to_datetime(en1)

    fig, axes = plt.subplots(figsize=(10, 5), nrows=2, sharex=True, sharey=True)
    for i, (st, en) in enumerate([sten0, [st1, en1]]):
        dfi = df0[((df0.dates>=st) & (df0.dates<=en))].copy()
        dfi['vel_u'] = dfi['vel_u'] - dfi['vel_u'].iloc[0]
        ratei, unci = calc_MIDAS.midas_fortran_one(sta, st, en, n_sig=None, plot=False)
        ratei*=1000; unci*=1000
        modi = np.polyval([ratei, 0], dfi['decyr'])
        modi -= modi[0]

        axes[i].scatter(dfi.dates, dfi['vel_u'], c='k', alpha=0.5, s=2)
        l, = axes[i].plot(dfi.dates, modi, color='r', linestyle='--', label=f'{ratei:.2f}$\\pm${unci:.2f}')
        axes[i].fill_between(dfi.dates, modi-unci, modi+unci, color=l.get_color(), alpha=0.3)
        axes[i].legend()

        axes[i].set_ylabel('VLM (mm)', fontsize=13)
        axes[i].legend()
        axes[i].grid(color='k', alpha=0.1, linestyle='--')
        axes[i].set_ylim([-50, 50])

    axes[0].set_title(f'Full {sta} ({sten0[0].strftime("%Y/%m/%d")} to {sten0[1].strftime("%Y/%m/%d")})', fontsize=14)
    axes[1].set_title(f'Truncated {sta} ({st1.strftime("%Y/%m/%d")} to {en1.strftime("%Y/%m/%d")})', fontsize=14)
    return


def plot_all_ts(ExpI, stas=[]):
    """ Wrapped to plot all timeseries in a region

    Updates the actual ts on the fly, shows original MIDAS rate of downloaded
    """
    TSObj = GPS_TS_Analysis(ExpI)
    df_gps = TSObj.df_gps
    stas = TSObj.df_gps['sta'].tolist() if not any(stas) else stas
    df_ts = dl_tss(df_gps, stas).reset_index()
    TSObj.plot_ts(df_ts, stas)
    return


def possible_exps():
    """ Just to abstract a bit away """
    # ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20)
    # ExpBase(Miami_SRc, 'ERA5_PM_ex_Fast', 'ZMA1')
    # ExpBase(HR_SR, 'ERA5_SET_PM_ex_Fast', 'SPVA', neofs=25)
    # ExpBase(Kiribati_SR, 'ERA5_SET_ex_Fast', 'KRTM', neofs=5)
    # ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=10)
    # ExpBase(Kennedy_SR, 'Base_ex_Fast', 'CCV6', neofs=15)
    # return ExpBase(DC_SR, 'ERA5_SET_PM_ex_Fast', 'USN7', neofs=20)


if __name__ == '__main__':
    ExpObj = ExpBase(DC_SR, 'ERA5_SET_PM_ex_Fast_2017', 'USN8', neofs=15)
    plot_all_ts(ExpObj)
    plt.show()
    os.sys.exit()
    # analyze_one_station()

    ### For updating with new MIDAS rates -------------------------------------
    TSObj = GPS_TS_Analysis(ExpObj)
    # TSObj.plot_nearby_stations(100, 'USN7')

    ## it's okay to write this over but if you want to update need to redownload
    # TSObj.update_MIDAS(tthresh=0.0, show=True, write=True, lst_force=['KRTM', 0.0, 2e-4])
    # TSObj.update_MIDAS(tthresh=0.0, show=True, write=False)
    # TSObj.update_MIDAS(tthresh=0.01, start='20200107', end_insar=True, show=True, write=False)
    # TSObj.update_MIDAS(tthresh=0.01, end_insar=True, show=True, write=False)
    TSObj.update_MIDAS(tthresh=0.00, start='20161229', end_insar=True, show=True, write=False)

    # TSObj.update_MIDAS(tthresh=0.01, end_insar=True, show=True, write=False)

    ### For comparing MIDAS rates and InSAR rates------------------------------
    ## plot the comparisons with InSAR
    RateObj = GPS_Rate_Analysis(ExpObj)
    # RateObj.plot_gps_vs_insar(npix=2, unc_lsq=True)
    # RateObj.plot_gps_vs_insar_qq(npix=0)

    ### Compare Rates before/after new MIDAS-----------------------------------
    # RateObj.plot_trends_new() # compare before/after

    #### Save
    path_figs = op.join(op.expanduser('~'), 'Desktop', 'VLM', f'{ExpObj.reg}_2024')
    # bbPlot.savefigs(path_figs, True, True)
    plt.show()
