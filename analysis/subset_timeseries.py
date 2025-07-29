"""
Cut a full timeseries into parts and compute the ratemaps and uncertainties
"""
from VLM.bzFRInGE import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from BZ import bbTS, bbPlot, bbGIS
from mintpy.utils import readfile, utils as ut
from mintpy.objects import timeseries
from mintpy.utils import time_func

import h5py, xarray as xr




def SliceTS(ExpBase):
    def __init__(self, dct_exp, mp_exp='Base_ex', ref_sta=None, neofs=''):
        super(). __init__(dct_exp, mp_exp, ref_sta, neofs)
        self.npix = npix
        self._set_stitched()
        self._set_timeseries()

        self.df_gps  = prep_gps(self.path_gps, self.reg).set_index('sta')

        self.mask = self._get_mask()

        arr_iang = readfile.read(self.path_geom_mp_geo, datasetName='incidenceAngle')[0]
        self.arr_iang = np.deg2rad(arr_iang)


    def __call__(self):
        pass


    def _rm_plate_motion(self):
        log.debug ('removing plate motion from timeseries')
        arr_pm, meta1 = readfile.read(op.join(self.path_mp_exp_geo, 'ITRF14.h5'))
        arr_pm = arr_pm.astype(np.float32)
        assert 'REF_X' in meta1.keys(), 'Need to reference the plate motion'
        arr_pm0= (1000 * arr_pm / np.cos(self.arr_iang).astype(np.float32)).reshape(-1)
        inters = np.zeros_like(arr_pm0, dtype=np.float32)
        # arr_pm = np.polyval([arr_pm.reshape(-1), np.zeros_like(arr_pm.reshape(-1))], decyr)
        arr_pm1 = (arr_pm0[:, np.newaxis] * self.decyr) + inters[:, np.newaxis]
        arr_pm  = arr_pm1.transpose(1, 0).reshape(self.arr_ts.shape)
        del arr_pm1, inters
        return arr_pm


    def _get_mask(self):
        with h5py.File(self.path_mask_vup, 'r') as h5:
            key = list(h5.keys())[0]
            mask = np.where(np.isclose(h5[key][:], 0), np.nan, 1)
        return mask


    def _set_timeseries(self, path_ts_geo=None):
        if path_ts_geo is not None:
            # ignore the rest of the if/else
            pass
        elif 'stitch' in self.mp_exp.lower():
            path_ts_geo  = self.path_ts_geo_stitch
        elif 'recon' in self.path_ts_geo.lower():
            log.info('Using reconstructed time series...')
            path_ts_geo  = self.path_ts_geo
        elif 'base' in self.path_ts_geo.lower():
            log.info('NOT using reconstructed time series or removing annual...')
            path_ts_geo  = self.path_ts_geo
        else:
            log.info('Using time series with annual cycle removed...')
            path_ts_geo  = self.path_ts_geo_ann
            # path_ts_geo  = self.path_ts_geo

        assert Path(path_ts_geo).exists(), f'Timeseries {path_ts_geo} does not exist!'

        # read the timeseries file with mintpy
        tsObj  = timeseries(path_ts_geo)
        arr_ts, dates = self._adjust_ts_sten(tsObj)

        path_excl = self.path_mp_exp / 'exclude_date.txt'
        if path_excl.exists() and not 'recon' in path_ts_geo:
            lst_excl = pd.read_csv(path_excl, header=None).squeeze().astype(str).tolist()
            log.info(f'Excluding {len(lst_excl)} dates.')
            keep_ix = [i for i, dt in enumerate(dates) if not dt in lst_excl]
            arr_ts, dates = arr_ts[keep_ix], list(np.array(dates)[keep_ix])

        self.arr_ts = arr_ts
        self.dates = dates
        self.dt = pd.to_datetime(self.dates)
        self.decyr  = bbTS.date2dec(self.dt)
        return


    def _adjust_ts_sten(self, tsobj):
        """ Adjust the full timeseries to the st/end of the experiment """
        meta = readfile.read_attribute(self.path_vlos_geo)
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