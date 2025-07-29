"""
Compare InSAR and Groundwater

Notes:
    There are no NGWMN in the NYC study area

"""
import xarray as xr
import h5py

from contrib import *
from contrib.experiments import *
from contrib.FRINGEBase import ExpBase
from BZ import bbTS
from mintpy.utils import readfile

from VLM.bzGW.GW_Analysis import *


def buff_wesn(wesn, lalo, buff=0.05):
    """ Expand wesn to include lalo """
    lat, lon   = lalo
    W, E, S, N = wesn
    if lat<S:
        S=lat-buff
    elif lat>N:
        N=lat+buff
    if lon<W:
        W=lon-buff
    elif lon>E:
        E=lon+buff

    return W,E,S,N


class LoadInSAR(ExpBase):
    def __init__(self, dct_exp, mp_exp='Base', ref_sta=None, neofs=''):
        super().__init__(dct_exp, mp_exp, ref_sta, neofs, v=True)
        self.ts, self.dates = self.load_ts_dt()
        self.decyr          = bbTS.date2dec(self.dates)
        self.ref_lalo       = self.set_ref_lalo()


    def load_ts_dt(self, mask=True):
        """ Load dates, timeseries and convert to mm """
        with h5py.File(self.path_ts_geo, 'r') as h5:
            arr_ts  = h5['timeseries'][:]
            # arr_tsa = h5['timeseries-ann'][:]
            arr_dt  = [dt.decode('utf-8') for dt in h5['date'][:]]

        dates = pd.to_datetime(arr_dt)
        arr0  = arr_ts*1000

        if mask:
            with h5py.File(self.path_mask_vup, 'r') as h5:
                mask = h5['mask'][:]
                arr0 *= np.where(np.isclose(mask, 0), np.nan, 1)
        return arr0, dates


    def load_iang(self):
        """ Load incidence angle as an array """
        with h5py.File(self.path_geom_mp_geo, 'r') as h5:
            arr_iang = h5['incidenceAngle'][:]
        return arr_iang


    def set_ref_lalo(self):
        ## average a bit around the reference point
        attrs    = readfile.read_attribute(self.path_ts_geo)
        return [float(attrs['REF_LAT']), float(attrs['REF_LON'])]



def InSAR_at_Wells(INObj, GWObj, rad=100, plot=True):
    ## convenience
    arr_vup     = INObj.ts / INObj.load_iang()
    df_wells    = GWObj.df_locs
    df_wells_ts = GWObj.df_ts

    ## change the column names
    if GWObj.name == 'NWIS':
        df_wells.rename(columns={'site_no':'SiteNo'}, inplace=True)

        df_wells_ts.rename(columns={'site_no':'SiteNo', 'lev_va': 'DTW',
                                            'dates': 'date'}, inplace=True)
        df_wells_ts['decyr'] = bbTS.date2dec(pd.to_datetime(df_wells_ts['date'].tolist()))

    df_wells.drop_duplicates('SiteNo', inplace=True)
    n_wells = len(df_wells_ts.SiteNo.unique())

    ## get the gps for converting to vertical
    df_gps  = prep_gps(INObj.path_gps, INObj.reg).set_index('sta')
    gps_vel = df_gps.loc[INObj.ref_sta, 'u_vel'] * 1000
    gps_sig = df_gps.loc[INObj.ref_sta, 'u_sig'] * 1000
    ts_gps  = np.polyval([gps_vel, 0], INObj.decyr)
    ts_gps -= ts_gps[0]


    ## for plotting
    da_insar  = xr.open_dataset(INObj.path_rate_msk_nc)['Band1'].assign_attrs(units='mm/yr')
    da_insar *=1000
    WESN      = [da_insar.lon.min(), da_insar.lon.max(), da_insar.lat.min(), da_insar.lat.max()]
    norm      = mpl.colors.TwoSlopeNorm(0, -3, 3)
    proj    = ccrs.PlateCarree()
    basemap = cimgt.GoogleTiles(url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
    # basemap   = cimgt.Stamen('terrain-background')

    n_skipped = 0
    dct_rates = {}
    for i, (site, df) in enumerate(df_wells_ts.groupby('SiteNo')):
        log.info('Site %s of %s', i, n_wells) if i % 5 == 0 else ''

        df_well  = df_wells[df_wells.SiteNo == site]
        if df_well.empty:
            continue

        lalo     = df_well.lat.item(), df_well.lon.item()
        ts_up    = avg_ts(arr_vup, lalo, rad, INObj.path_rate_nc)*1000
        if np.isnan(ts_up).all():
            n_skipped += 1
            continue

        ts_up   -= ts_up[0]
        ts_up    = ts_up + ts_gps # shift by GPS rate

        wesn     = buff_wesn(WESN, lalo)

        ## insar trend
        coeffs_up, cov_up = np.polyfit(INObj.decyr, ts_up, 1, cov=True)
        trend_up   = np.polyval(coeffs_up, INObj.decyr)
        sem        = np.sqrt(cov_up[0,0])

        ## gw trend (flipped from DTW to Head and mm)
        coeffs_gw  = np.polyfit(df.decyr, df.DTW, 1)
        rate_gw    = -1000 * coeffs_gw[0]

        ## rates
        dct_rates[site] = (coeffs_up[0], rate_gw, coeffs_up[0]-rate_gw)

        ## plot location and eventually timeseries below
        if plot:
            fig = plt.figure()
            gs = fig.add_gridspec(3, 3)
            axeMap    = fig.add_subplot(gs[0:2, :], projection=ccrs.PlateCarree())
            # fig, axes = plt.subplots(figsize=(12,9),
            #                         subplot_kw={'projection': basemap.crs})
            axeMap.set_extent(wesn)
            axeMap.add_image(basemap, 9)
            im   = axeMap.pcolormesh(da_insar.lon, da_insar.lat, da_insar, shading='auto',
                    norm=norm, transform=proj, cmap='cmc.roma_r')
            # bbPlot.cartopy_cbar(im, 'Vup (mm/yr)')

            axeMap.scatter(lalo[1], lalo[0], transform=proj, c='w', marker='v', s=100, zorder=20)
            gl  = axeMap.gridlines(draw_labels=True)
            gl.top_labels, gl.right_labels = False, False
            gl.xlines, gl.ylines = None, None

            axTS = fig.add_subplot(gs[2, :])
            axTS.scatter(INObj.dates, ts_up, color='k', marker='.', zorder=20)
            axTS.set_ylabel('Vup (mm)')

            axGW = axTS.twinx()
            axGW.scatter(df['date'], df['DTW'], color='darkgreen', marker='s', s=5)
            axGW.set_ylabel('GW Head (mm)', rotation=270, color='darkgreen', labelpad=10)
            plt.show()

    log.info('Skipped %s of %s wells outside of InSAR coverage', n_skipped, n_wells)
    if n_skipped == n_wells:
        log.critical('No %s wells within InSAR coverage!', GWObj.name)
        return

    ## move this elsewhere
    df_rates = pd.DataFrame.from_dict(dct_rates, orient='index', columns=['rate_Vup', 'rate_GW', 'resid'])
    fig, axes = plt.subplots(figsize=(10, 5))
    axes.scatter(df_rates['rate_Vup'], df_rates['rate_GW'], color='k')
    axes.plot([0,1],[0,1], 'r--', transform=axes.transAxes)
    axes.set_xlabel('InSAR Vup (mm/yr)')
    axes.set_ylabel('Groundwater Extraction (mm/yr)')
    fig.set_label(f'Vup_GW_{INObj.reg}')
    return


if __name__ == '__main__':
    # NWIS   = NWISBase(region, styr=2015, show=True, sub=True)

    region = 'NYC'
    # InExp  = LoadInSAR(Charleston_SR, 'Base', 'SCHA', neofs=15)
    InExp  = LoadInSAR(NYC_SR, 'Base', 'NYBP', neofs=15)
    NGWMN  = NGMWNBase(region, styr=2000, show=False)
    NWIS   = NWISBase(region, styr=2000, show=False, sub=True)

    # InSAR_at_Wells(InExp, NGWMN, plot=True)
    InSAR_at_Wells(InExp, NWIS, plot=True)
    bbPlot.savefigs(PATH_RES, True, True)
    plt.show()
