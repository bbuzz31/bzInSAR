from VLM.bzFRInGE import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from BZ import bbTS, bbPlot, bbGIS
from mintpy.utils import readfile, utils as ut
from mintpy.objects import timeseries
from mintpy.utils import time_func

import h5py, xarray as xr

DCT_STABLE = {
    ## coped from Stitch_ex_Fast
    'Base_ex_Fast' : {
        'LaGuardia_Stable': [(40.782755, -73.853617), 'NYC'],
        'Williams_Stable': [(40.712456, -73.919430), 'NYC'],
        'Ashe_Stable': [(40.75025, -73.825000), 'NYC'],
        'Woodside_Stable': [(40.742738, -73.901400), 'NYC'],
       'Jefferson_Memorial_Stable': [(38.8861, -77.0321), 'DC'], # 1/2 km
       'Anacostia_Stable': [(38.86777, -76.99155), 'DC'], # 1/2 km
        'Alexandria_Stable': [(38.81554, -77.05599), 'DC'],
        },

    'ERA5_SET_PM_ex_Fast' : {
        'Navy_Yard_Stable': [( 38.87693, -76.99905), 'DC'],
        # 'Anacostia_Stable': [(38.86749, -76.99155), 'DC'], # 675 m
        'Anacostia_Stable': [(38.84249, -76.99543), 'DC'], # 3km

        'NRL1_Stable': [(38.82138, -77.02460), 'DC'],
        'GSFC_Stable': [(38.99610, -76.85155), 'DC'], # 150 m
        'MD_Up_Stable': [(38.83999, -76.96127), 'DC'],
        'MD_Down_Stable': [(38.84082, -76.95377), 'DC'], # 600 m, 0.11
        'Lincoln_Memorial_Stable': [(38.89416, -77.04905), 'DC'], #500 m
        'Jefferson_Memorial_Stable': [(38.88610, -77.03016), 'DC'], #900 m
        'DC_TG_Stable': [(38.86971, -77.01821), 'DC'], # 550 m
        'Alexandria_Stable': [(38.77971, -77.06599), 'DC'], # 3 km!!
        'Auth_Village_Stable': [(38.81638, -76.90793), 'DC'],
        'Smithsonian_Stable': [(38.88832, -77.02071), 'DC'],
        },

    'ERA5_SET_PM_ex_Fast_2020' : {
        'Alexandria_Stable': [(38.80693, -77.06293), 'DC'], # 1.1 km
        'Auth_Village_Stable': [(38.81888, -76.90960), 'DC'], # 500 m
    },

    'ERA5_SET_PM_Stitch_ex_Fast' : {
           'LaGuardia_Stable': [(40.782755, -73.853617), 'NYC'],
           'Williams_Stable': [(40.712456, -73.919430), 'NYC'],
           'Ashe_Stable': [(40.75025, -73.825000), 'NYC'],
           'Woodside_Stable': [(40.742738, -73.901400), 'NYC'], # Paper, Coherence 0.88
        #    'Woodside_Stable': [(40.74303, -73.90166), 'NYC'] # ~500 m
        #    'Woodside_Stable': [(40.7199, -73.92305), 'NYC'] # ~ 3 km
           },

    'ERA5_SET_PM_Stitch_ex_Fast2017' : {
           'LaGuardia_Stable': [(40.773866, -73.843329), 'NYC'],
        #    'LaGuardia_Stable': [(40.792477, -73.849440), 'NYC'], #option2
           'Williams_Stable': [(40.714653, -73.915824), 'NYC'],
           'Ashe_Stable': [(40.773866, -73.843329), 'NYC'],
           'Woodside_Stable': [(40.742191, -73.898324), 'NYC']
           },

    'Base_ex_Fast_20150310_20180318' : {
        'Alexandria_Stable': [(38.80415,-77.05849), 'DC']
    },
    'Base_ex_Fast_20180318_20210326' : {
        'Alexandria_Stable': [(38.80193,-77.05266), 'DC']
    },
    'Base_ex_Fast_20210326_20240427' : {
        'Alexandria_Stable': [(38.80693,-77.06266), 'DC']
    },

    'ERA5_SET_PM_ex_Fast_20150310_20180318' : {
        'Alexandria_Stable': [(38.81388,-77.05488), 'DC']
    },

    'Base_ex_Fast_2017': {
       'Jefferson_Memorial_Stable': [(38.88554, -77.03516), 'DC'], # 1/2 km
       'Anacostia_Stable': [(38.86749, -76.99155), 'DC'], # 1/2 km
       'Alexandria_Stable': [(38.81971, -77.05516), 'DC'], # 2 km
    },
    'Base_ex_Fast_2017_20200822': {
       'Alexandria_Stable': [(38.80388, -77.05155), 'DC'], # 200 m
    },
    'Base_ex_Fast_20200822': {
       'Alexandria_Stable': [(38.81388, -77.05571), 'DC'], # 1.5 km
    },

    'ERA5_SET_PM_ex_Fast_2017': {
       'Jefferson_Memorial_Stable': [(38.88610, -77.03016), 'DC'], # 1 km
       'Anacostia_Stable': [(38.86721, -76.99127), 'DC'], # 1/2 km
       'Alexandria_Stable': [(38.82277, -77.02266), 'DC'], # 4 km
    },
    'ERA5_SET_PM_ex_Fast_2017_20200822': {
       'Alexandria_Stable': [(38.76721, -77.06016), 'DC'], # 4 km
    },
    'ERA5_SET_PM_ex_Fast_20200822': {
       'Alexandria_Stable': [(38.83721, -77.05016), 'DC'], # 4 km
    },

    'ERA5_SET_PM_ex_Fast_20180318_20210326' : {
        # 'Alexandria_Stable': [(38.80082,-77.05127), 'DC']
        # 'Alexandria_Stable': [(38.80082,-77.05183), 'DC']
        # 'Alexandria_Stable': [(38.80082,-77.04238), 'DC'] # 1km
        'Alexandria_Stable': [(38.77554,-77.05710), 'DC'] # 1km
    },
    'ERA5_SET_PM_ex_Fast_20210326_20240427' : {
            'Alexandria_Stable': [(38.82443, -77.04349), 'DC'] # 3km
        }
}


DCT_PTS = {
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
           'Navy_Yard': [( 38.874669, -76.994311), 'DC'],
           'Anacostia': [(38.867515, -76.985479), 'DC'],
           'Anacostia1': [(38.8675, -76.9853), 'DC'],
           'NRL1': [(38.82075, -77.0244), 'DC'],
           'GSFC': [(38.9949, -76.8523), 'DC'],
           'MD_Up': [(38.840696, -76.961077), 'DC'],
           'MD_Down': [(38.836976, -76.95), 'DC'],
           'GW_Cemetary': [(39.01317, -76.968211), 'DC'],
           'Lincoln_Memorial': [(38.88943, -77.04988), 'DC'],
        #    'Jefferson_Memorial': [(38.8814, -77.0365), 'DC'],
           'Jefferson_Memorial': [(38.8812, -77.0365), 'DC'],
           'DC_TG': [(38.873333,  -77.021667), 'DC'],
           'Alexandria': [(38.801951, -77.051810), 'DC'],
           'Auth_Village': [(38.821286, -76.906076), 'DC'],
           'Smithsonian': [(38.8860, -77.0214), 'DC'],
           'Ivy_City': [(38.9148, -76.9859), 'DC'],
           'Q': [(38.924389, -77.002139), 'DC'], #GW Well

           'Backliff': [(29.51284, -94.98520), 'Houston'],
           'BBayou': [(29.7372, -95.168658), 'Houston'],
           'La_Porte':[(29.652929, -95.011463), 'Houston'],
           'Galveston':[(29.294410, -94.788261), 'Houston'],
           'Galveston_Pier':[(29.31049, -94.79169), 'Houston'],
           'Friendship_Park': [(29.580402, -95.010165), 'Houston'],
           'Fresno': [(29.53886, -95.44744), 'Houston'],
           'Cusch_Terrace': [(29.737298, -95.019553), 'Houston'],
           'Clear_Lake': [(29.55756, -95.121729), 'Houston'],
           'Santa_Fe_Oil_Field': [(29.378951 -95.102045), 'Houston'],
           'Sugar_Land': [(29.598259, -95.623712), 'Houston'],
           'Teal_Run': [(29.523700, -95.467621), 'Houston'],
        #    'Texas_City': [(29.385771, -94.907259), 'Houston'],
           'Texas_City': [(29.370969, -94.903802), 'Houston'],
           'Channel_View': [(29.777341, -95.141080), 'Houston'],
           'Nonlinear': [(29.9225, -95.375), 'Houston'],
           }


def check_rate(rate_true, rate_ts, npix, error=False):
    """ Check rates calculate by hand (rate_ts) match those in velocity file """
    log.setLevel('DEBUG')
    if npix == 0:
        log.debug (f'"Real velocity": {rate_true:.3f}')
        log.debug (f'New velocity: {rate_ts:.3f}')
        log.info (f'Difference between new and "real" velocity: '\
            f'{np.abs(rate_ts - rate_true):.3f} mm/yr')

        if error:
            assert np.abs(rate_true - rate_ts) < 1e-3, \
            'Velocity estimated from time series does not match velocity from uvel'

    else:
        log.info (f'"Real velocity": {rate_true:.2f}')
        log.info (f'New velocity: {rate_ts:.2f}')
        log.info (f'Difference between new and "real" velocity: '\
               f'{np.abs(rate_ts - rate_true):.2f} mm/yr')
    return


def check_stable(exp):
    """ Wrapped """
    Obj = PlotTS1(exp.dct_exp, exp.mp_exp0, exp.ref_sta, neofs=exp.neofs)
    assert Obj.reg == 'NYC', print ('Only NYC supported')

    for pt in 'LaGuardia Williams Woodside Ashe'.split():
        Obj.plot_ts_vup(f'{pt}_Stable')
        # calculate the distance
        print (f'Distance to stable pt {pt}', np.round(
                        bbGIS.distance_pt2pt(DCT_PTS[pt][0],
                        DCT_STABLE[Obj.mp_exp0][f'{pt}_Stable'][0],
                        units='km'), 2), 'km')

    # plt.close('all')
    return


class PlotTS1(ExpBase):
    def __init__(self, dct_exp, mp_exp='Base_ex', ref_sta=None, neofs='', npix=0, sten=''):
        super(). __init__(dct_exp, mp_exp, ref_sta, neofs)
        self.npix = npix
        self.set_sten_paths(sten)
        self._set_stitched()
        self._set_timeseries()

        self.df_gps  = prep_gps(self.path_gps, self.reg).set_index('sta')
        self.gps_vel = self.df_gps.loc[self.ref_sta, 'u_vel']
        self.gps_sig = self.df_gps.loc[self.ref_sta, 'u_sig']

        self.mask = self._get_mask()

        arr_iang = readfile.read(self.path_geom_mp_geo, datasetName='incidenceAngle')[0]
        self.arr_iang = np.deg2rad(arr_iang).astype(np.float32)


    def calc_ts_vup_poly(self, gser, radius=500, corr_thresh=0.8, rm_atm=False):
        """ Calculate the average timeseries within a polygon """
        polygon = gser.geometry.to_numpy()
        log.debug('reading velocity and uncertainty from file: %s', self.path_vup_geo)
        meta = readfile.read(self.path_vup_geo)[1]
        mask = False if (self.mask == 1).all() else True
        da_vel0, da_unc0 = self.get_rate_unc_nc(mask)

        # spatial means of aoi
        vel_true =  da_vel0.rio.clip(polygon, all_touched=True).mean().item()
        unc = da_unc0.rio.clip(polygon, all_touched=True).mean().item()

        ## convert whole timeseries to Vup, subtract potential plate motion (mm/yr)
        arr_pm = self._rm_plate_motion() if 'PM' in self.mp_exp.upper() else 0
        arr  = 1000 * self.arr_ts / np.cos(self.arr_iang)
        arr  = arr - arr_pm

        ## get the GPS velocity everywhere
        gps_vel = meta.get('REF_RATE', self.gps_vel) * 1000
        gps_sig = meta.get('REF_UNC', self.gps_sig) * 1000

        # make a 'timeseries' from the GPS velocity and add everywhere
        ts_gps  = np.polyval([gps_vel, 0], self.decyr) #+ ref_ts
        ts_gps -= ts_gps[0]
        arr += ts_gps[:, np.newaxis, np.newaxis]

        ## convert to xarray for slicing
        arr_ts = np.fliplr(arr) if meta['ORBIT_DIRECTION'] == 'ASCENDING' else arr

        da_ts = xr.DataArray(arr_ts, name='timeseries', dims=['time', 'lat', 'lon'],
                            coords={'time': self.dt, 'lat': da_vel0.lat, 'lon': da_vel0.lon}
                            ).assign_attrs(units='m').rio.write_crs(4326)

        loc_ts = da_ts.rio.clip(polygon, all_touched=True).mean('lat lon'.split())
        ts_vup = (loc_ts - loc_ts[0])

        coeffs_vup, rss = np.polyfit(self.decyr-self.decyr[0], ts_vup, 1, full=True)[:2]
        rate_vup0 = coeffs_vup[0]
        rmse0 = np.sqrt(rss/(len(ts_vup)-2)).item()
        std0 = np.nanstd(ts_vup)
        check_rate(vel_true, rate_vup0, '', error=False)

        ## could use mintpy instead
        # preds, coeffs_vup, rss = time_func.estimate_time_func(
        #     {'polynomial':1}, list(self.dates), ts_vup)
        # coeffs_vup = coeffs_vup[::-1]

        ## not likely to work because of the averaging
        if rm_atm:
            ## get the pixel location for selecting atm
            lon_idx = da_ts.sel(lon=gser.centroid.x.item(), method='nearest').lon.values
            lat_idx = da_ts.sel(lat=gser.centroid.y.item(), method='nearest').lat.values
            x = da_ts.lon.values.tolist().index(lon_idx)
            y = da_ts.lat.values.tolist().index(lat_idx)
            npix = int((radius * 2 * np.sqrt(2))/30) # approximate number of 30 m pixels
            arr_crop = arr[:, y-npix:y+npix, x-npix:x+npix]
            ts_vup, coeffs_vup, unc_dd = self._rm_atm(arr_crop, ts_vup,
                                            (radius, corr_thresh), show_ts=False)

            rss = np.polyfit(self.decyr, ts_vup, 1, full=True)[1]
            # root mean difference around trend
            rmse = np.sqrt(rss/(len(ts_vup)-2)).item()
            unc = np.sqrt(unc**2 + unc_dd**2) # should figure out a better way

            log.info ('Original Polygon Rate, Timeseries Std, RMSE:, %.2f / %.1f / %.1f (mm)',
                        rate_vup0[0], std0, rmse0)

            log.info ('Double Diff Polygon Rate, Timeseries Std, RMSE: %.2f / %.1f / %.1f (mm)',
                    coeffs_vup[0], np.nanstd(ts_vup), rmse)

            ## the polyfit is too low
            # sem2 = np.sqrt(cov_vup[0,0]+gps_sig**2)

        del arr, arr_pm
        return ts_vup, coeffs_vup, unc


    def calc_ts_vup(self, loc, lalo=None):
        """ Calc the vertical displacement timeseries average (meters) around one loc

        Give a loc (from dict above) or specific lalo
        Give a ref_sta to use one other than the default
        """
        log.debug('reading velocity and uncertainty from file: %s', self.path_vup_geo)
        mask = False if (self.mask == 1).all() else True
        vel_true0, unc0, _, meta = self.get_rate_unc_h5(mask, self.path_vup_geo)
        vel_true0 *= 1000
        unc0 *= 1000

        # convert whole timeseries to Vup, subtract potential plate motion (mm/yr)
        arr = 1000 * self.arr_ts / np.cos(self.arr_iang)
        arr_pm = self._rm_plate_motion() if 'PM' in self.mp_exp.upper() else 0
        arr = arr - arr_pm

        ## get the GPS velocity everywhere
        gps_vel = meta.get('REF_RATE', self.gps_vel) * 1000
        gps_sig = meta.get('REF_UNC', self.gps_sig) * 1000

        # make a 'timeseries' from the GPS velocity and add everywhere
        ts_gps  = np.polyval([gps_vel, 0], self.decyr) #+ ref_ts
        ts_gps -= ts_gps[0]
        arr += ts_gps[:, np.newaxis, np.newaxis]

        if lalo is None:
            lalo = DCT_STABLE[self.mp_exp0][loc][0] if 'stable' in loc.lower() else DCT_PTS[loc][0]

        # average around the reference point (bbox or radius)
        # ref_lalo = [float(meta['REF_LAT']), float(meta['REF_LON'])]
        ref_y, ref_x  = int(meta['REF_Y']), int(meta['REF_X'])


        ## the pixel at the actual reference GPS is 0; can average around it a bit
        # ref_ts  = arr[:, ref_y, ref_x]
        # ref_ts  = np.nanmean(arr[:, ref_y-npix_ref:ref_y+npix_ref, ref_x-npix_ref:ref_x+npix_ref], axis=(1, 2)) # bbox
        # ref_ts1  = avg_ts(arr, ref_lalo, rad, self.path_rate_nc) # use a radius

        ## optionally average a bit around the target point
        coord    = ut.coordinate(meta, lookup_file=self.path_geom_mp_geo)
        y, x     = coord.geo2radar(lalo[0], lalo[1])[0:2]

        if self.npix == 0 or (loc is not None and 'stable' in loc.lower()):
            # always use 0 around stable
            loc_ts   = arr[:, y, x] # no real diff if surrounding points similar
            vel_true = vel_true0[y, x]
            unc      = unc0[y, x]
            assert not np.isnan(self.mask[y,x]), f'{loc} Pixel y={y} x={x} is masked!'
        else:
            loc_ts = arr[:, y-self.npix:y+self.npix, x-self.npix:x+self.npix]
            assert not np.isnan(loc_ts).all(), f'All {self.npix} pixels surrounding y={y} x={x} are masked!'
            loc_ts = np.nanmean(loc_ts, axis=(1,2))
            vel_true = np.nanmean(vel_true0[y-self.npix:y+self.npix, x-self.npix:x+self.npix])
            unc      = np.nanmean(unc0[y-self.npix:y+self.npix, x-self.npix:x+self.npix])

        if loc is not None and np.abs(vel_true) > 0.5 and 'stable' in loc.lower():
            raise Exception('Stable velocity is too high. Probably not correct stable point for')


        # DONT remove the timeseries at the GPS reference from the target point
            # remember timeseries is not shifted to GPS
            # but is referenced, so theoretically ref_ts should be 0
            # this mostly adds noise so I'M NOT DOING IT
        # ts_loc   = loc_ts + ref_ts

        try:
            npix_ref = int(meta['n_pix'])
            if npix_ref > 0:
                log.warning('New and calculated velocity will be off due to prior averaging during referencing')
                # gps_vel = np.nanmean(vel_true0[ref_y-npix_ref:ref_y+npix_ref, ref_x-npix_ref:ref_x+npix_ref]) # bbox
                # gps_sig = np.nanmean(unc0[ref_y-npix_ref:ref_y+npix_ref, ref_x-npix_ref:ref_x+npix_ref]) # bbox
        except:
            log.warning(f"Couldn't get number of pixels used for referencing from metada")


        ts_vup     = (loc_ts - loc_ts[0])

        ## use mintpy fitting for exact match (1e-4 mm/yr diff)
        # coeffs_vup, cov_vup = np.polyfit(decyr, ts_vup, 1, cov='unscaled')
        # rss = np.polyfit(decyr, ts_vup, 1, full=True)[1].item() # residual sum of squares
        # preds, coeffs_vup1 = bbTS.bzSeason.fit_lsq(decyr, ts_vup) # seasonal
        # trend_vup  = preds[:, 0] ## show seasonal

        preds, coeffs_vup, rss = time_func.estimate_time_func(
            {'polynomial':1}, list(self.dates), ts_vup)
        coeffs_vup[::-1] = coeffs_vup # reverse order for matching later
        rmse = np.sqrt(rss / (len(ts_vup)-2)) # mm

        log.info (f'Loc: {loc} y/x: %s/%s', y, x)
        npix = 0 if 'stable' in loc.lower() else self.npix
        check_rate(vel_true, coeffs_vup[0], npix, error=False)

        ## shown in double diff
        # log.info ('Rate, Timeseries Std, RMSE: %.2f / %.1f / %.1f (mm)',
        #             coeffs_vup[0], np.nanstd(ts_vup), rmse)

        ## the polyfit is too low
        # sem2 = np.sqrt(cov_vup[0,0]+gps_sig**2)
        sem2 = unc

        del arr, arr_pm
        return ts_vup, coeffs_vup, sem2


    def calc_ts_vup_dd_stable(self, loc_target, loc_stable=None, show_ts=False):
        """ Plot the vertical displacement timeseries average  around one loc

        Doubled differenced; plot 'loc_target' after removing 'loc_stable'
        Note that plate motion may give subtle differences in stable trend

        loc_target/stable can be a single name (gets lalo from dictionary) or a tuple with name / lalo
        """
        log.debug('reading velocity and uncertainty from file: %s', self.path_vup_geo)
        mask = False if (self.mask == 1).all() else True
        vel_true0, unc0, _, meta = self.get_rate_unc_h5(mask, self.path_vup_geo)
        vel_true0 *= 1000
        unc0 *= 1000

        arr_pm = self._rm_plate_motion() if 'PM' in self.mp_exp.upper() else 0

        ## convert whole timeseries to Vup, subtract potential plate motion (mm/yr)
        arr = 1000 * self.arr_ts / np.cos(self.arr_iang)
        arr = arr - arr_pm

        loc_stable = f'{loc_target}_Stable' if loc_stable is None else loc_stable
        ts_vup_target, coeffs_vup_target, unc_target = self.calc_ts_vup(loc_target)
        rate_target, bias_target = coeffs_vup_target
        rss0 = np.polyfit(self.decyr, ts_vup_target, 1, full=True)[1]
        rmse0 = np.sqrt(rss0/(len(ts_vup_target)-2)).item()

        ts_vup_stable, coeffs_vup_stable, unc_stable = self.calc_ts_vup(loc_stable)
        rate_stable, bias_stable = coeffs_vup_stable

        ## double difference
        ts_dd    = ts_vup_target - ts_vup_stable
        # ser = pd.Series(ts_dd, index=dates, color='')
        # import seaborn as sns
        # coeffs_dd = np.polyfit(decyr, ts_dd, 1)
        # trend_dd  = np.polyval(coeffs_dd, decyr)
        # sns.lmplot(data=df, x='decyr', y='ts', scatter=True,
                #    scatter_kws={'color': 'k'}, line_kws={'color': 'red'})

        ## show the two timeseries and diff
        if show_ts:
            fig, axes1 = plt.subplots(figsize=(10, 4))
            axes1.scatter(self.dt, ts_vup_target, label=f'{loc_target}', color='darkblue', s=15, alpha=0.25)
            axes1.scatter(self.dt, ts_vup_stable, label=f'{loc_stable}', color='darkgreen', s=15, alpha=0.25)
            axes1.scatter(self.dt, ts_dd, label=f'DDiff', color='black', s=15, alpha=0.75)
            axes1.set_ylabel('VLM (mm)', fontsize=13)
            # axes1.plot(dates, trend_dd, linestyle='--', label=f'DDiff: {coeffs_dd[0]:.2f} mm/yr', color='black')
            axes1.legend()
            axes1.set_title(f'Target/Stable Rate: {rate_target:.2f} / {rate_stable:.2f} mm/yr', fontsize=14)
            axes1.grid(color='k', linestyle='--', alpha=0.1)

        log.info (f'Vup Rate Target: {rate_target:.2f}, '\
                f'Vup Rate Stable: {rate_stable:.2f}')

        rate_dd, unc_dd, bias_dd = bbTS.bootstrap_1D(ts_dd, self.decyr, n_boots=1500)
        rss = np.polyfit(self.decyr, ts_dd, 1, full=True)[1]
        rmse = np.sqrt(rss/(len(ts_dd)-2)).item()

        # eventually need a better way of calculating the uncertinaty
        unc_gps = np.mean([unc_target**2, unc_stable**2])
        sem2       = np.sqrt(unc_dd**2+unc_gps)
        coeffs_dd, cov_vup = np.polyfit(self.decyr, ts_dd, 1, cov='unscaled')
        # check_ts_unc(dates, ts_dd) # other experiments

        # compute rho, remove trends and compute rho as additional check
        rho0 = np.corrcoef(ts_vup_target, ts_vup_stable)[0,1]
        ts_vup_target_detrend = ts_vup_target - np.polyval(coeffs_vup_target, self.decyr-self.decyr[0])
        ts_vup_stable_detrend = ts_vup_stable - np.polyval(coeffs_vup_stable, self.decyr-self.decyr[0])
        rho_detrend = np.corrcoef(ts_vup_target_detrend, ts_vup_stable_detrend)[0,1]
        show_corr = False
        if show_corr:
            fig, axes = plt.subplots(figsize=(10, 6), nrows=2, sharex=True)
            axes[0].scatter(self.dt, ts_vup_target, s=15, color='k')
            axes[0].scatter(self.dt, ts_vup_target_detrend, s=15, color='r')
            axes[0].set_title(loc_target)
            axes[1].scatter(self.dt, ts_vup_stable, s=15, color='k')
            axes[1].scatter(self.dt, ts_vup_stable_detrend, s=15, color='r')
            axes[1].set_title(loc_stable)
            for ax in axes:
                ax.set_ylabel('Vertical Displacement (mm)')
                ax.grid(color='k', alpha=0.1, linestyle='--')
        # using seasonal
        # preds, coeffs_vup = bbTS.bzSeason.fit_lsq(decyr, ts_dd)
        # coeffs_vup = coeffs_dd[-2:]
        # trend_vup = preds[:, 0]
        # rate_dd  = coeffs_dd[0]
        log.info ('Original Rate, Timeseries Std, RMSE:, %.2f / %.1f / %.1f (mm)',
                    coeffs_vup_target[0], np.nanstd(ts_vup_target), rmse0)

        log.info ('Double Diff Rate, Timeseries Std, RMSE: %.2f / %.1f / %.1f (mm)',
                    coeffs_dd[0], np.nanstd(ts_dd), rmse)

        log.info(f'Correlation between {loc_target} and {loc_stable}: {rho0:.2f}')
        log.info (f'Correlation between DETRENDED {loc_target} and {loc_stable}: {rho_detrend:.2f}')

        return ts_dd, coeffs_dd, sem2


    def calc_ts_vup_dd_nn(self, loc, show_ts=True):
        """ Calculate dd using neareast neighbor """
        mask = False if (self.mask == 1).all() else True
        vel_true0, unc0, _, meta = self.get_rate_unc_h5(mask, self.path_vup_geo)
        vel_true0 *= 1000
        unc0 *= 1000
        arr = 1000 * self.arr_ts / np.cos(self.arr_iang)
        arr_pm = self._rm_plate_motion() if 'PM' in self.mp_exp.upper() else 0
        arr = arr - arr_pm

        ## get the GPS velocity everywhere
        gps_vel = meta.get('REF_RATE', self.gps_vel) * 1000
        gps_sig = meta.get('REF_UNC', self.gps_sig) * 1000

        # make a 'timeseries' from the GPS velocity and add everywhere
        ts_gps  = np.polyval([gps_vel, 0], self.decyr) #+ ref_ts
        ts_gps -= ts_gps[0]
        arr += ts_gps[:, np.newaxis, np.newaxis]

        ## get the actual timeseries at the location of interest
        coord = ut.coordinate(meta, lookup_file=self.path_geom_mp_geo)
        lalo = DCT_PTS[loc][0]
        y, x = coord.geo2radar(lalo[0], lalo[1])[0:2]

        ## get the ts of interest
        if self.npix == 0 or (loc is not None and 'stable' in loc.lower()):
            # always use 0 around stable
            loc_ts   = arr[:, y, x] # no real diff if surrounding points similar
            vel_true = vel_true0[y, x]
            unc      = unc0[y, x]
            # loc_ts_pm = arr[:, y, x]
            # rate_pm, bias_pm = np.polyfit(self.decyr, loc_ts_pm, 1)
            assert not np.isnan(self.mask[y,x]), f'{loc} Pixel y={y} x={x} is masked!'

        else:
            loc_ts = arr[:, y-self.npix:y+self.npix, x-self.npix:x+self.npix]
            assert not np.isnan(loc_ts).all(), f'All {self.npix} pixels surrounding y={y} x={x} are masked!'
            loc_ts = np.nanmean(loc_ts, axis=(1,2))
            vel_true = np.nanmean(vel_true0[y-self.npix:y+self.npix, x-self.npix:x+self.npix])
            unc      = np.nanmean(unc0[y-self.npix:y+self.npix, x-self.npix:x+self.npix])

        loc_ts_stable = pd.Series()
        while loc_ts_stable.empty:
            for yy in range(-2, 2):
                for xx in range(-2,2):
                    if yy == 0 and xx == 0: # dont get same point
                        continue
                    if np.isnan(self.mask[y+yy, x+xx]):
                        continue
                    else:
                        loc_ts_stable = pd.Series(arr[:, y+yy, x+xx])
                        break
        assert not loc_ts_stable.empty, 'No nearby unmasked points!'

        loc_ts -= loc_ts[0]
        loc_ts_stable -= loc_ts_stable[0]
        ## remove the trends and then put back

        coeffs0, rss0 = np.polyfit(self.decyr, loc_ts, 1, full=True)[:2]
        ts_fit = np.polyval(coeffs0, self.decyr)
        rmse0 = np.sqrt(rss0/(len(loc_ts)-2)).item()
        check_rate(vel_true, coeffs0[0], self.npix, error=False)

        coeffs1, rss1 = np.polyfit(self.decyr, loc_ts_stable, 1, full=True)[:2]
        ts_fit_stable = np.polyval(coeffs1, self.decyr)
        rmse1 = np.sqrt(rss1/(len(loc_ts_stable)-2)).item()

        # first just try and remove the
        ts_dd = (loc_ts - (loc_ts_stable - ts_fit_stable)).to_numpy()
        ts_dd -= ts_dd[0]
        rate_dd, unc_dd, bias_dd= bbTS.bootstrap_1D(ts_dd, self.decyr, n_boots=2000)
        fit_dd = np.polyval([rate_dd, bias_dd], self.decyr)

        rss = np.polyfit(self.decyr, ts_dd, 1, full=True)[1]
        # root mean difference around trend
        rmse = np.sqrt(rss/(len(ts_dd)-2)).item()
        sem2 = np.sqrt(unc**2 + unc_dd**2) # should figure out a better way

        log.info ('Original Rate, Timeseries Std, RMSE:, %.2f / %.1f / %.1f (mm)',
                    coeffs0[0], np.nanstd(loc_ts), rmse0)

        log.info ('Double Diff Rate, Timeseries Std, RMSE: %.2f / %.1f / %.1f (mm)',
                    rate_dd, np.nanstd(ts_dd), rmse)

        if show_ts:
            fig, axes = plt.subplots(figsize=(10, 6), nrows=2, sharex=True)
            axes[0].scatter(self.dt, loc_ts, color='k', s=15)
            axes[0].plot(self.dt, ts_fit, color='k', linestyle='--',
                         label=f"'Target': {coeffs0[0]:.1f} (mm/yr)")
            axes[0].scatter(self.dt, loc_ts_stable, s=15, color='gray')
            axes[0].plot(self.dt, ts_fit_stable, color='gray', linestyle='--',
                         label=f"'Stable': {coeffs1[0]:.1f} (mm/yr)")

            axes[1].scatter(self.dt, ts_dd, s=15, color='darkred')
            axes[1].plot(self.dt, fit_dd, color='darkred', label=f"'DDiff': {rate_dd:.1f} (mm/yr)")
            for ax in axes:
                ax.set_ylabel('V$_{Up}$ (mm)', fontsize=12)
                ax.legend()
                ax.grid(color='k', linestyle='--', alpha=0.1)
            axes[0].set_title(loc, fontsize=14)

        # log.info(f'Correlation between {loc_target} and {loc_stable}: {rho0:.2f}')
        # log.info (f'Correlation between DETRENDED {loc_target} and {loc_stable}: {rho_detrend:.2f}')
        return ts_dd, [rate_dd, bias_dd], sem2


    def calc_ts_vup_dd_atm(self, loc, radius, corr_thresh=0.95, show_ts=True):
        """Calculate double-differenced vertical displacement time series by removing atmospheric effects

        Uses correlation-based filtering to identify and remove atmospheric signals from the time series
        at a given location by finding highly correlated pixels within a specified radius.

        Args:
            loc (str): Location name to analyze
            radius (float): Radius in meters to search for correlated pixels
            corr_thresh (float, optional): Correlation threshold for identifying atmospheric signals.
                                         Defaults to 0.95.
            show_ts (bool, optional): Whether to plot the time series. Defaults to True.

        Returns:
            tuple: (ts_dd, coeffs_dd, sem2) containing:
                - ts_dd: Double-differenced time series with atmospheric effects removed
                - coeffs_dd: Coefficients [rate, bias] of linear fit to ts_dd
                - sem2: Combined uncertainty from original and atmospheric correction
        """
        # convert whole timeseries to Vup, subtract potential plate motion (mm/yr)
        mask = False if (self.mask == 1).all() else True
        vel_true0, unc0, _, meta = self.get_rate_unc_h5(mask, self.path_vup_geo)
        vel_true0 *= 1000
        unc0 *= 1000
        arr = 1000 * self.arr_ts / np.cos(self.arr_iang)
        arr_pm = self._rm_plate_motion() if 'PM' in self.mp_exp.upper() else 0
        arr = arr - arr_pm
        # Multiply timeseries array by mask to handle invalid pixels
        arr = arr * self.mask[np.newaxis, :, :]

        ## get the GPS velocity everywhere
        gps_vel = meta.get('REF_RATE', self.gps_vel) * 1000
        gps_sig = meta.get('REF_UNC', self.gps_sig) * 1000

        # make a 'timeseries' from the GPS velocity and add everywhere
        ts_gps  = np.polyval([gps_vel, 0], self.decyr) #+ ref_ts
        ts_gps -= ts_gps[0]
        arr += ts_gps[:, np.newaxis, np.newaxis]
        ## write the adjusted timeseries to disk
        # import shutil
        # src = Path(self.path_ts_geo).parent / f'{Path(self.path_ts_geo).stem}_ITRF14-PM.h5'
        # shutil.copy(self.path_ts_geo, dst)
        # with h5py.File(dst, 'r+') as h5:
        #     del h5['timeseries']
        #     h5['timeseries'] = arr.astype(np.float32)/1000
        # print (f'Wrote: {dst}')

        ## get the actual timeseries at the location of interest
        coord = ut.coordinate(meta, lookup_file=self.path_geom_mp_geo)
        lalo = DCT_PTS[loc][0]
        y, x = coord.geo2radar(lalo[0], lalo[1])[0:2]

        if self.npix == 0 or (loc is not None and 'stable' in loc.lower()):
            # always use 0 around stable
            loc_ts   = arr[:, y, x] # no real diff if surrounding points similar
            vel_true = vel_true0[y, x]
            unc      = unc0[y, x]
            # loc_ts_pm = arr[:, y, x]
            # rate_pm, bias_pm = np.polyfit(self.decyr, loc_ts_pm, 1)
            assert not np.isnan(self.mask[y,x]), f'{loc} Pixel y={y} x={x} is masked!'
        else:
            loc_ts = arr[:, y-self.npix:y+self.npix, x-self.npix:x+self.npix]
            assert not np.isnan(loc_ts).all(), f'All {self.npix} pixels surrounding y={y} x={x} are masked!'
            loc_ts = np.nanmean(loc_ts, axis=(1,2))
            vel_true = np.nanmean(vel_true0[y-self.npix:y+self.npix, x-self.npix:x+self.npix])
            unc = np.nanmean(unc0[y-self.npix:y+self.npix, x-self.npix:x+self.npix])

        ts_vup_target = loc_ts - loc_ts[0]
        coeffs0, rss0 = np.polyfit(self.decyr, ts_vup_target, 1, full=True)[:2]
        rmse0 = np.sqrt(rss0/(len(ts_vup_target)-2)).item()
        check_rate(vel_true, coeffs0[0], self.npix, error=False)

        # approximate number of 30 m pixels on either side
        npix = int((radius * 2 * np.sqrt(2))/30/2)
        arr_crop = arr[:, y-npix:y+npix, x-npix:x+npix]
        ts_dd, coeffs_dd, unc_dd = self._rm_atm(
                    arr_crop, ts_vup_target, (radius, corr_thresh), show_ts)


        rss = np.polyfit(self.decyr, ts_dd, 1, full=True)[1]
        # root mean difference around trend
        rmse = np.sqrt(rss/(len(ts_dd)-2)).item()
        sem2 = np.sqrt(unc**2 + unc_dd**2) # should figure out a better way

        log.info ('Original Rate, Timeseries Std, RMSE:, %.2f / %.1f / %.1f (mm)',
                    coeffs0[0], np.nanstd(ts_vup_target), rmse0)

        log.info ('Double Diff Rate, Timeseries Std, RMSE: %.2f / %.1f / %.1f (mm)',
                    coeffs_dd[0], np.nanstd(ts_dd), rmse)

        r2 = bbTS.calc_r2(ts_dd, np.polyval(coeffs_dd, self.decyr))
        log.info(f'R^2={r2:.2f}')

        # log.info(f'Correlation between {loc_target} and {loc_stable}: {rho0:.2f}')
        # log.info (f'Correlation between DETRENDED {loc_target} and {loc_stable}: {rho_detrend:.2f}')
        return ts_dd, coeffs_dd, sem2


    def plot_ts_vup(self, loc, lalo=None, axes=None, show_fit=True):
        ts_vup, coeffs_vup, unc = self.calc_ts_vup(loc, lalo)
        trend_vup  = np.polyval(coeffs_vup, self.decyr-self.decyr[0])

        # show title if not in an inset
        show_title = True if axes is None else False
        if axes is None:
            fig, axes = plt.subplots(figsize=(10, 4))
        else:
            fig = axes.get_figure()

        col = 'r' if coeffs_vup[0] > 0 else 'b'
        axes.scatter(self.dt, ts_vup, c='k', s=15)#, label=f'$\\sigma$={np.std(ts_vup):.2f}')
        # axes.scatter(self.dt, ts_vup, alpha=0, label=f'$RMSE$={rmse:.2f}')

        lbl=f'{coeffs_vup[0]:.2f}$\\pm${unc:.2f} (mm/yr)'
        if show_fit:
            axes.plot(self.dt, trend_vup, color=col, linestyle='--', label=lbl)
            axes.fill_between(self.dt, trend_vup-unc, trend_vup+unc, color=col, alpha=0.3)
            axes.legend()
        axes.set_ylabel('Vertical Displacement (mm)', fontsize=13)
        axes.grid(color='gray', linestyle = '--', linewidth=0.1)

        if loc in DCT_PTS:
            axes.set_title(f'{loc}, {DCT_PTS[loc][1]}')
        elif loc in DCT_STABLE[self.mp_exp0]:
            axes.set_title(f'{loc}, {DCT_STABLE[self.mp_exp0][loc][1]}')
        else:
            axes.set_title(f'{loc}')

        if not show_title:
            axes.set_title('')

        fig.set_label(f'{loc}_ts')

        df_ts = pd.DataFrame({'ts': ts_vup, 'trend': trend_vup}, index=self.dt)
        return fig, axes, df_ts


    def plot_ts_vup_dd(self, loc_target, axes=None, radius=500, corr_thresh=0.90,
                       show_ts=False, show_fit=True):
        if isinstance(loc_target, (np.ndarray, pd.Series)):
            ts_dd = loc_target
            loc_target = 'Grad'

        log.info(f'---{loc_target} (npix={self.npix})---\n')
        loc_stable = f'{loc_target}_Stable'

        if self.reg == 'NYC':
            # For NYC, use stable reference method
            ts_dd, coeffs_dd, sem2 = self.calc_ts_vup_dd_stable(loc_target, loc_stable, show_ts)
        elif loc_target == 'Grad':
            # For gradient data, calculate rate using bootstrap
            rate1, unc1, bias1 = bbTS.bootstrap_1D(
               ts_dd, self.decyr, n_boots=2500)
            coeffs_dd = [rate1, bias1]
            sem2 = unc1
        else:
            # For other regions, use atmospheric correction method
            ts_dd, coeffs_dd, sem2 = self.calc_ts_vup_dd_atm(loc_target,
                                    radius, corr_thresh, show_ts=show_ts)

        rate_dd, bias_dd = coeffs_dd
        trend_vup = np.polyval(coeffs_dd, self.decyr)

        # center on mean instead of 0
        ts_dd -= np.nanmean(ts_dd)
        trend_vup -= np.nanmean(trend_vup)

        if axes is None:
            fig, axes = plt.subplots(figsize=(10, 3.5))
        else:
            fig = axes.get_figure()

        col = 'r' if rate_dd > 0 else 'b'
        col = 'k' if len(coeffs_dd) > 2 else col # nonlinear fits have >2 coeffs

        # Plot scatter points for dd time series
        axes.scatter(self.dt, ts_dd, c='k', s=15)

        # Optionally show trend line and uncertainty band
        if show_fit:
            axes.plot(self.dt, trend_vup, color=col, linestyle='--',
                    label=f'{rate_dd:.1f}$\\pm${sem2:.1f} (mm/yr)')
            axes.fill_between(self.dt, trend_vup-sem2, trend_vup+sem2, color=col, alpha=0.3)
            lloc = 'lower right' if loc_target == 'Williams' else 'best'
            axes.legend(prop ={"size": 13}, loc=lloc)

        axes.set_ylabel('Vertical Displacement (mm)', fontsize=13)
        axes.tick_params(axis='x', labelsize=13)

        ymin, ymax = int(np.floor(ts_dd.min()))-5, int(np.ceil(ts_dd.max())+5)
        ymin, ymax = (-10, 10) if ymin > -10 and ymax < 10 else (ymin, ymax)
        axes.set_ylim([ymin, ymax])

        axes.grid(color='k', linestyle = '--', alpha=0.1)
        fig.set_label(f'{loc_target}_dd_{self.npix}pix')
        df_ts = pd.DataFrame({'ts': ts_dd, 'trend': trend_vup}, index=self.dt)
        return fig, axes, df_ts


    def plot_ts_vup_piece2(self, en='20191231', loc='Anacostia',
                              double_diff=True, npix=1, axes=None):
        from dateutil import relativedelta
        if double_diff:
            if self.reg == 'NYC':
                ts_vup = self.calc_ts_vup_dd_stable(loc, show_ts=False)[0]
            if isinstance(loc, np.ndarray):
                ts_vup = loc
                loc = 'Grad'
            else:
                ts_vup = self.calc_ts_vup_dd_atm(loc, radius=500, corr_thresh=0.9, show_ts=False)[0]
        else:
            ts_vup = self.calc_ts_vup(loc)[0]

        df_ts = pd.DataFrame({'ts':ts_vup, 'decyr': self.decyr}, index=self.dt)

        if axes is None:
            fig, axes = plt.subplots(figsize=(10, 4))
            axes.set_title(f'{loc}-{loc}_Stable')
        else:
            fig = axes.get_figure()

        for i in range(2):
            if i == 0:
                df_ts1 = df_ts[df_ts.index<=en]
            elif i == 1:
                st = pd.to_datetime(en) + relativedelta.relativedelta(days=1)
                df_ts1 = df_ts[df_ts.index>=st]

            rate1, unc1, bias1 = bbTS.bootstrap_1D(
                df_ts1['ts'].to_numpy(), df_ts1['decyr'].to_numpy(), n_boots=2500)

            unc1 = np.sqrt(unc1**2 + (self.gps_sig*1000)**2)

            if np.abs(rate1) < unc1:
                col = 'dimgray'
            else:
                col = 'darkred' if rate1 > 0 else 'darkblue'

            trend_vup1 = np.polyval([rate1, bias1], df_ts1['decyr'])

            axes.scatter(df_ts1.index, df_ts1['ts'], c=col, s=15)
            axes.plot(df_ts1.index, trend_vup1, color=col, linestyle='--',
                    label=f'{rate1:.1f}$\\pm${unc1:.1f} (mm/yr)')
            axes.fill_between(df_ts1.index, trend_vup1-unc1, trend_vup1+unc1,
                              color=col, alpha=0.3)
        axes.legend(prop ={"size": 13}, loc='best')
        axes.set_ylabel('Vertical Displacement (mm)', fontsize=13)
        axes.tick_params(axis='x', labelsize=13)
        fig.set_label(f'{loc}_dd_{npix}pix')

        log.info(f'{loc} (npix={npix})')

        ymin, ymax = int(np.floor(ts_vup.min()))-5, int(np.ceil(ts_vup.max())+5)
        ymin, ymax = (-10, 10) if ymin > -10 and ymax < 10 else (ymin, ymax)
        axes.set_ylim([ymin, ymax])

        axes.grid(color='k', linestyle = '--', alpha=0.1)
        fig.set_label(f'tsdd_{loc}_{npix}')
        return fig, axes


    def plot_moving_rate(self, loc_target, axes=None, last_st='20220101'):
        """ Plot rates computed in sucessive months (starting at beginning) """
        ts_dd, coeffs_dd, sem2 = self.calc_ts_vup_dd_atm(loc_target, radius=500,
                                                         corr_thresh=0.9, show_ts=False)
        lst_rates, lst_uncs, lst_sts = [], [], []
        for i in range(len(self.decyr)):
            st = self.dt[i]
            if st >= pd.to_datetime(last_st):
                break
            ts = ts_dd[i:]
            decyr = self.decyr[i:]
            coeffs = np.polyfit(decyr, ts, 1)

            rate, unc, bias = bbTS.bootstrap_1D(ts, decyr, n_boots=1500)

            lst_rates.append(rate)
            lst_uncs.append(unc)
            lst_sts.append(self.dt[i])

        arr_rates, arr_uncs = np.array(lst_rates), np.array(lst_uncs)
        fig, axes = plt.subplots(figsize=(10, 4))
        # axes.scatter(lst_sts, lst_rates, color='k', s=15)
        l, = axes.plot(lst_sts, lst_rates, color='k', alpha=1.0)
        axes.fill_between(lst_sts, arr_rates-arr_uncs, arr_rates+arr_uncs,
                          color=l.get_color(), alpha=0.3)
        axes.set_ylabel('Vup (mm/yr)', fontsize=13)
        axes.set_xlabel('Start Date', fontsize=13)
        axes.grid(color='k', linestyle='--', alpha=0.1)
        axes.set_title(f'{loc_target} Moving Trends')
        fig.set_label(f'{loc_target}_moving')
        return fig, axes


    def plot_ts_models_1D(self, loc):
        """ Compute the fit of the raw data and check there is seasonality """
        assert not 'recon' in self.path_vup_geo, 'Dont use reconstructed timeseries'
        self._set_timeseries(rm_ann=False)
        mask = False if (self.mask == 1).all() else True

        vel_true0, unc0, _, meta = self.get_rate_unc_h5(mask, self.path_vup_geo)
        annAmp_true0 = readfile.read(self.path_vup_geo, datasetName='annualAmplitude')[0]*1000
        vel_true0 *= 1000
        unc0 *= 1000

        # convert whole timeseries to Vup, subtract potential plate motion (mm/yr)
        arr = 1000 * self.arr_ts# / np.cos(self.arr_iang)
        arr_pm = self._rm_plate_motion() if 'PM' in self.mp_exp.upper() else 0
        # arr = arr - arr_pm


        ## get the GPS velocity everywhere
        gps_vel = meta.get('REF_RATE', self.gps_vel) * 1000
        gps_sig = meta.get('REF_UNC', self.gps_sig) * 1000

        # make a 'timeseries' from the GPS velocity and add everywhere
        # ts_gps  = np.polyval([gps_vel, 0], self.decyr) #+ ref_ts
        # ts_gps -= ts_gps[0]
        # arr += ts_gps[:, np.newaxis, np.newaxis]

        lalo = DCT_PTS[loc][0]

        ## get target point
        coord    = ut.coordinate(meta, lookup_file=self.path_geom_mp_geo)
        y, x     = coord.geo2radar(lalo[0], lalo[1])[0:2]

        # always use 0 pixels around target
        loc_ts   = arr[:, y, x] # no real diff if surrounding points similar
        vel_true = vel_true0[y, x]
        unc      = unc0[y, x]
        annAmp   = annAmp_true0[y, x]
        assert not np.isnan(self.mask[y,x]), f'{loc} Pixel y={y} x={x} is masked!'

        ts_vup  = (loc_ts - loc_ts[0])

        ## use mintpy fitting for exact match (1e-4 mm/yr diff)
        tis = 'Full Model', 'Linear Fit', 'Linear + Annual', 'Linear + semiAnnual'
        fig, axes = plt.subplots(figsize=(10, 6), nrows=4, sharex=True, sharey=True)
        for i, mod in enumerate([[1.0, 0.5], [], [1.0], [0.5]]):
            G, coeffs_vup, rss = time_func.estimate_time_func(
                {'polynomial':1, 'periodic': mod}, list(self.dates), ts_vup)
            rmsei = np.sqrt(rss / (len(ts_vup)-(2+len(mod)))) # mm

            log.info (f'Loc: {loc} y/x: %s/%s', y, x)
            npix = 0 if 'stable' in loc.lower() else self.npix
            check_rate(vel_true, coeffs_vup[1], npix, error=False)
            predsi = np.dot(G, coeffs_vup)

            ax = axes[i]
            ax.scatter(self.dt, ts_vup, color='k', s=15)
            ax.plot(self.dt, predsi, color='k', label=f'RMSE={rmsei:.1f} (mm)')
            ax.grid(color='k', linestyle='--', alpha=0.1)
            ax.set_ylabel('V$_{{Up}}$', fontsize=12)
            ax.set_title(tis[i], fontsize=14)
            ax.legend()

            if mod:
                G1, coeffs_vup1, rss1 = time_func.estimate_time_func(
                    {'polynomial':1, 'periodic': mod}, list(self.dates), ts_vup)
                ## annual amplitude; wont match mintpy
                annAmp1 = np.sqrt(coeffs_vup1[2]**2 + coeffs_vup1[3]**2)
                if i == 0:
                    breakpoint()
                    log.info(f'{tis[i]} Seasonal Amplitude={annAmp1:.2f} (mm)')

        fig.subplots_adjust(hspace=0.6)
        return


    def plot_ts_models_2D(self):
        """ Compare RMSE of linear fit vs period fits everywhere """
        from VLM.bzFRInGE.plotting.plotExp import PlotExp
        assert not 'recon' in self.path_vup_geo, 'Dont use reconstructed timeseries'
        self._set_timeseries(rm_ann=False)
        mask = False if (self.mask == 1).all() else True

        vel_true0, unc0, _, meta = self.get_rate_unc_h5(mask, self.path_vup_geo)
        vel_true0 *= 1000
        unc0 *= 1000

        # convert whole timeseries to Vup, subtract potential plate motion (mm/yr)
        arr = 1000 * self.arr_ts / np.cos(self.arr_iang)
        arr_pm = self._rm_plate_motion() if 'PM' in self.mp_exp.upper() else 0
        arr = arr - arr_pm

        ## get the GPS velocity everywhere
        gps_vel = meta.get('REF_RATE', self.gps_vel) * 1000
        gps_sig = meta.get('REF_UNC', self.gps_sig) * 1000

        # make a 'timeseries' from the GPS velocity and add everywhere
        ts_gps  = np.polyval([gps_vel, 0], self.decyr) #+ ref_ts
        ts_gps -= ts_gps[0]
        arr += ts_gps[:, np.newaxis, np.newaxis]
        arr1 = np.where(np.isnan(arr), 0, arr).reshape(len(self.dates), -1)

        PObj = PlotExp(self)
        PObj.pparms['norm'] = mpl.colors.Normalize(0, 15)
        PObj.pparms['cmap'] = 'cmo.amp'

        ## use mintpy fitting for exact match (1e-4 mm/yr diff)
        tis = 'Full Model', 'Linear Fit', 'Linear + Annual', 'Linear + semiAnnual'

        fig, axes = plt.subplots(figsize=(10, 10), ncols=2, nrows=2,
                                 sharex=True, sharey=True,
                                 subplot_kw={'projection': PObj.basemap.crs})

        for i, mod in enumerate([[1.0, 0.5], [], [1.0], [0.5]]):
            G, coeffs_vup, rss = time_func.estimate_time_func(
                {'polynomial':1, 'periodic': mod}, list(self.dates), arr1)
            rmse = np.sqrt(rss / (len(self.dates)-(2+len(mod)))).reshape(vel_true0.shape)
            rmse = np.where(np.isnan(vel_true0), np.nan, rmse)
            da_rmse = PObj.da.copy()
            da_rmse.data = np.flipud(rmse) if meta['ORBIT_DIRECTION'] == 'ASCENDING' else rmse
            ax = axes.ravel()[i]
            ax.add_image(PObj.basemap, PObj.zoom)
            im   = ax.pcolormesh(da_rmse.lon, da_rmse.lat, da_rmse, shading='nearest',
                            alpha=PObj.alpha, transform=PObj.proj, **PObj.pparms)
            bbPlot.cartopy_cbar(im, ylabel=f'{tis[i]} RMSE (mm)', fontsize=PObj.fs1, pad=0.15)#, extend='min')
            log.info (f'{tis[i]} RMSE: {da_rmse.mean():.2f} +/- {da_rmse.std():.2f}')

        fig.subplots_adjust(hspace=0.25, wspace=0.25)
        fig.set_label('Tmod_compare')
        return


    def _rm_atm(self, arr_crop, ts_vup_target, alg_parms, show_ts=True):
        """Remove atmospheric signals from a time series using correlation-based filtering.

        This function identifies and removes atmospheric signals by:
        1. Detrending the target time series
        2. Finding highly correlated pixels within a specified radius
        3. Computing the mean atmospheric signal from correlated pixels
        4. Subtracting the atmospheric signal from the original time series

        Args:
            arr_crop (ndarray): 3D array of time series data around target point (referenced to GPS)
            ts_vup_target (ndarray): Target time series to remove atmosphere from
            alg_parms (tuple): (radius, corr_thresh) parameters for correlation filtering
                radius (float): Search radius in meters around target point
                corr_thresh (float): Correlation threshold for identifying atmospheric signals
            show_ts (bool, optional): Whether to plot diagnostic figures. Defaults to True.

        Returns:
            tuple: (ts_dd, [rate_dd, bias_dd], unc_dd)
                ts_dd (ndarray): Atmosphere-corrected time series
                rate_dd (float): Linear rate of corrected time series
                bias_dd (float): Y-intercept of corrected time series
                unc_dd (float): Uncertainty of corrected time series rate

        Notes:
            - Detrending before correlation has minimal impact compared to not detrending
            - Only keeps pixels with correlation >= corr_thresh and < 1.0 (excludes self)
            - Uses bootstrap to estimate uncertainty of corrected time series
        """
        radius, corr_thresh = alg_parms

        ## detrend the actual timeseries
        coeffs0 = np.polyfit(self.decyr, ts_vup_target, 1)
        trend0 = np.polyval(coeffs0, self.decyr)
        ts_vup_target_detrend = ts_vup_target - trend0

        ## subset around the target location
        arr_crop1 = arr_crop.reshape(len(self.decyr), -1)
        coeffs = np.polyfit(self.decyr, arr_crop1, 1)
        trend = np.polyval(coeffs, self.decyr[:, None])
        ## detrend everything or dont
        arr_crop1_detrend = (arr_crop1 - trend).T

        # eventually could do this and fit seasonal after keeping in EOF
        # for ts in arr_crop1.T:
        #     pred, coeffs1 = bbTS.bzSeason.fit_lsq(self.decyr, ts)
        lst_corrs = []
        for ts in arr_crop1_detrend:
            lst_corrs.append(np.corrcoef(ts_vup_target_detrend, ts)[0, 1])

        ser_corrs = pd.Series(lst_corrs).dropna()
        ser_corrs = ser_corrs[ser_corrs>=corr_thresh]
        ser_corrs = ser_corrs[ser_corrs<1.0] # drop itself
        assert not ser_corrs.empty, 'No points found with high enough correlation ' \
                        f'{corr_thresh} in {radius} m radius! Try lowering '\
                        'thresh or increasing radius.'
        log.info (f'Computing average of {len(ser_corrs)} timeseries with rho>={corr_thresh} in {radius} m radius')

        # compute the mean over the local well correlated pixels
        arr_atm_ts = np.nanmean(arr_crop1_detrend[ser_corrs.index], 0)
        # subtract the atmosphere from the original (with trend) timeseries
        ts_dd = ts_vup_target - arr_atm_ts
        ts_dd -= ts_dd[0]
        ## compute a new trend and uncertainty
        rate_dd, unc_dd, bias_dd = bbTS.bootstrap_1D(ts_dd, self.decyr, n_boots=1500)

        # show the original, detrended, and atmosphere
        if show_ts:
            trend_dd = np.polyval([rate_dd, bias_dd], self.decyr)
            fig, axes = plt.subplots(figsize=(10, 7), nrows=4, sharex=True)
            axes[0].scatter(self.dt, ts_vup_target, color='k', s=15)
            axes[0].plot(self.dt, trend0, color='k', linestyle='--',
                    label=f'Orig: {coeffs0[0]:.2f}$\\pm${unc_dd:.2f} (mm/yr)')
            axes[1].scatter(self.dt, ts_vup_target_detrend, color='black', s=15, label='Detrended', marker='s')
            axes[1].scatter(self.dt, arr_atm_ts, color='darkgreen', s=15, label='Atmosphere')
            ax1a = axes[1].twinx()
            ax1a.set_ylabel(f'$\\rho$={np.corrcoef(ts_vup_target_detrend, arr_atm_ts)[0,1]:.2f}',
                          rotation=270, labelpad=10)

            ## show the residual between the double diff and the linear model
            axes[2].set_title('Residuals (TS-Trends)', fontsize=15)
            axes[2].scatter(self.dt, ts_vup_target - trend0, color='darkgreen', s=15)
            axes[2].set_ylabel('Orig')
            ax2a = axes[2].twinx()
            ax2a.set_ylabel(f'$\\mu$={(ts_vup_target-trend0).mean():.2f}', rotation=270, labelpad=10)

            axes[3].scatter(self.dt, ts_dd-trend_dd, color='darkgreen', s=15)
            axes[3].set_ylabel('DDiff')
            # axes[3].set_ylim(axes[2].get_ylim())
            ax3a = axes[3].twinx()
            ax3a.set_ylabel(f'$\\mu$={(ts_dd-trend_dd).mean():.2f}', rotation=270, labelpad=10)

            for ax in [ax1a, ax2a, ax3a]:
                ax.set_yticks([])
                ax.set_yticklabels([])

            for i, ax in enumerate(axes):
                if i < 2:
                    ax.set_ylabel('V$_{Up}$ (mm)', fontsize=12)
                    ax.legend()
                if i >= 2:
                    ax.axhline(0, color='k', linewidth=1)
                ax.grid(color='k', linestyle='--', alpha=0.1)
            fig.subplots_adjust(wspace=0.5)

        return ts_dd, [rate_dd, bias_dd], unc_dd


    def _rm_plate_motion(self):
        if hasattr(self, 'arr_pm') and self.arr_pm is not None:
            return self.arr_pm
        log.debug('removing plate motion from timeseries')
        arr_pm, meta1 = readfile.read(op.join(self.path_mp_exp_geo, 'ITRF14.h5'))
        arr_pm = arr_pm.astype(np.float32)
        assert 'REF_X' in meta1.keys(), 'Need to reference the plate motion'

        arr_pm0 = (1000 * arr_pm / np.cos(self.arr_iang)).reshape(-1).astype(np.float32)
        inters = np.zeros_like(arr_pm0, dtype=np.float32)

        arr_pm1 = (arr_pm0[:, np.newaxis] * self.decyr.astype(np.float32)) + inters[:, np.newaxis]
        arr_pm = arr_pm1.transpose(1, 0).reshape(self.arr_ts.shape)
        del arr_pm1, inters
        self.arr_pm = arr_pm.astype(np.float32)

        return self.arr_pm


    def _get_mask(self):
        with h5py.File(self.path_mask_vup, 'r') as h5:
            key = list(h5.keys())[0]
            mask = np.where(np.isclose(h5[key][:], 0), np.nan, 1)
        return mask


    def _set_timeseries(self, path_ts_geo=None, rm_ann=False):
        if path_ts_geo is not None:
            # ignore the rest of the if/else
            pass
        elif 'stitch' in self.mp_exp.lower():
            path_ts_geo  = self.path_ts_geo_stitch
        elif 'recon' in str(self.path_ts_geo).lower():
            log.info(f'Using reconstructed time series at {self.path_ts_geo}...')
            path_ts_geo  = self.path_ts_geo
        elif 'base' in str(self.path_ts_geo).lower():
            log.info('NOT using reconstructed time series or removing annual...')
            path_ts_geo  = self.path_ts_geo
        else:
            if rm_ann:
                log.info('Using time series with annual cycle removed...')
                path_ts_geo  = self.path_ts_geo_ann
            else:
                log.info(f'Using time series WITH annual cycle at {self.path_ts_geo}...')
                path_ts_geo  = self.path_ts_geo

        assert Path(path_ts_geo).exists(), f'Timeseries {path_ts_geo} does not exist!'

        # read the timeseries file with mintpy
        tsObj  = timeseries(path_ts_geo)
        arr_ts, dates = self._adjust_ts_sten(tsObj)

        path_excl = self.path_mp_exp / 'exclude_date.txt'
        if path_excl.exists() and not 'recon' in str(path_ts_geo).lower():
            lst_excl = pd.read_csv(path_excl, header=None).squeeze().astype(str).tolist()
            log.info(f'Excluding {len(lst_excl)} dates.')
            keep_ix = [i for i, dt in enumerate(dates) if not dt in lst_excl]
            arr_ts, dates = arr_ts[keep_ix], list(np.array(dates)[keep_ix])

        self.arr_ts = arr_ts.astype(np.float32)
        self.dates = dates
        self.dt = pd.to_datetime(self.dates)
        self.decyr  = bbTS.date2dec(self.dt)
        return


    def _adjust_ts_sten(self, tsobj):
        """ Adjust the full timeseries to the st/end of the experiment """
        assert self.path_vup_geo, f'Cannot find: {self.path_vup_geo}'
        meta = readfile.read_attribute(self.path_vup_geo)
        st0, en0 = meta['START_DATE'], meta['END_DATE']
        st0, en0 = pd.to_datetime([st0, en0])

        arr_ts = tsobj.read(print_msg=True)
        dates  = np.array(tsobj.get_date_list())
        dt = pd.to_datetime(dates)
        n0 = len(dates)
        ix = (dt>=st0) & (dt<=en0)
        n1 = len(dates[ix])
        if n0 != n1:
            log.info(f'Dropped {n0-n1} dates outside of st/en range.')

        return arr_ts[ix, :, :], dates[ix].tolist()


    def find_stable(self, lalo, mask=True, alpha=0.8, loc=None, plot=True):
        """ Find a stable point around another one. Loc just for naming

        Radius in meters around target point
        """
        import geopandas as gpd
        from shapely.geometry import Point, Polygon
        import xarray as xr
        import warnings
        warnings.filterwarnings("ignore", category=UserWarning)
        # get the insar rates
        epsg    = 'epsg:4326'
        path_nc = self.path_rate_msk_nc_stitch if 'stitch' in self.mp_exp.lower() \
            else self.path_rate_msk_nc
        da_rate = xr.open_dataset(path_nc)['Band1'].rio.write_crs(4326).rename('velocity')

        # if mask:
        #     mask = np.flipud(readfile.read(Exp.path_mask_vup)[0])

        #     da_rate *= mask
        #     da_rate  = da_rate.where(~np.isclose(da_rate, 0), np.nan)

        lat, lon = lalo
        # look around the target point
        gser = gpd.GeoSeries(gpd.points_from_xy([lon], [lat]), crs=epsg)
        rad  = 4e3 # 3 km search radius
        poly = gser.to_crs(4087).buffer(rad, cap_style=1).to_crs(epsg).geometry
        da_pt = da_rate.rio.clip(poly, epsg, drop=True)

        # first cut of mostly high values; can cut again later
        thresh = 0.25 / 1000 # keep only values +/- thresh mm
        da_pt  = da_pt.where(np.abs(da_pt) < thresh, np.nan)
        da_pt *= 1000 # convert to mm
        # da_pt = da_pt.where(np.abs(da_pt) < 0.001, np.nan)
        if da_pt.isnull().all():
            log.error ('No valid points! Try increasing search radius or values')
            return

        df = da_pt.to_dataframe().reset_index().abs().dropna()
        df.lon *= -1
        gdf = bbGIS.df2gdf(df, 'all')

        # get distance from target point to potential
        gdf_pt0 = gpd.GeoDataFrame(geometry=gpd.points_from_xy([lon], [lat], crs=4326))
        gdf_m  = gpd.sjoin_nearest(gdf.to_crs(4087), gdf_pt0.to_crs(4087),
                distance_col='distance').drop(columns='spatial_ref index_right'.split()).to_crs(4326)
        center = gdf_m.dissolve().centroid

        ## get time series noise ---
        if 'stitch' in self.mp_exp.lower():
            if self.reg != 'NYC':
                raise Exception ('Stitching not working for multiple stations (ok for NYC)')
            log.info('Using stitched time series and velocity...')
            path_ts_geo = self.path_ts_geo_stitch
        else:
            path_ts_geo = self.path_ts_geo

        meta = readfile.read_attribute(path_ts_geo)
        coord = ut.coordinate(meta, lookup_file=self.path_geom_mp_geo)

        ts_std = []
        for ix, row in gdf.iterrows():
            y, x = coord.geo2radar(row.lalo.y, row.lalo.x)[0:2]
            ts_std.append(np.nanstd(self.arr_ts[:, y, x]*1000))

        gdf_m['ts_std'] = ts_std
        gdf_m = gdf_m.sort_values('distance velocity'.split())
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

        print (gdf_m.head(10))
        breakpoint()

        return gdf_m


    def _set_stitched(self):
        if 'stitch' in self.mp_exp.lower():
            if self.reg != 'NYC':
                raise Exception ('Stitching not working for multiple stations (ok for NYC)')
            log.info('Using stitched time series and velocity...')
            self.path_vup_geo = self.path_vup_geo_stitch
            self.path_ts_Geo = self.path_ts_geo_stitch
        return


if __name__=='__main__':
    # Obj = PlotTS1(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20, npix=0)
    # Obj = PlotTS1(DC_SR, 'ERA5_SET_PM_ex_Fast_2017', 'USN8', neofs=20, npix=1)#, sten='_20161229-20200822')
    # Obj = PlotTS1(DC_SR, 'ERA5_SET_PM_ex_Fast_2017', 'USN8', neofs=0, npix=0)
    # Obj = PlotTS1(DC_SR, 'Base_ex_Fast_2017', 'USN8', neofs=0, npix=0)
    # Obj.calc_ts_vup_dd_nn('Anacostia', show_ts=True)


    # Obj.calc_ts_vup()
    # Obj.plot_ts_models_1D('Q')
    # Obj.plot_ts_models_1D('Anacostia')
    # Obj.plot_ts_models_1D('Anacostia')

    # Obj.plot_ts_models_2D()
    # Obj.plot_ts_vup('Anacostia')

    # Obj.plot_ts_vup_dd('Jefferson_Memorial', show_ts=True)
    # Obj.plot_ts_vup_dd('Anacostia', show_ts=True)
    # Obj.plot_ts_vup_dd('Alexandria', show_ts=True)

    # Obj.calc_ts_vup_dd_nn('Anacostia', show_ts=True)
    # Obj.plot_ts_vup_dd('Anacostia', show_ts=True)

    # Obj.plot_ts_vup_dd('MD_Up')
    # Obj.plot_moving_rate('MD_Up')

    # Obj.plot_ts_vup_dd('MD_Down')
    # Obj.plot_moving_rate('MD_Down')

    # for loc in 'Alexandria GSFC Smithsonian MD_Up MD_Down'.split():
        # Obj.plot_ts_vup_dd(loc, show_ts=False)

    Obj = PlotTS1(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=5, npix=2)
    df_ts = Obj.calc_ts_vup_dd_atm('Texas_City', radius=2.5e3, corr_thresh=0.95)

    # Obj.plot_ts_vup_dd('BBayou')
    Obj.plot_ts_vup_dd('BBayou')
    # Obj.plot_ts_vup_dd('Galveston_Pier')
    # Obj.plot_ts_vup_dd('Backliff')
    # Obj.plot_ts_vup_dd('BBayou')
    # bbPlot.savefigs(Path.home() / 'VLM', True, True)

    plt.show()
