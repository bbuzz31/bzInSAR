"""
Plot NYC Maps in Python for Paper
"""
from matplotlib.colors import TwoSlopeNorm, BoundaryNorm, Normalize
import matplotlib.patheffects as pe
from mintpy.utils import readfile

from BZ import bbLogger, bbGIS, bbPlot
from BZ.bbGIS import city_shape

from VLM import *
from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from VLM.bzFRInGE.analysis import analyzeGPS

DCT_CITIES    = {'NYC': ['Manhattan', 'Brooklyn', 'Queens', 'Bronx',
                         'Hudson', 'Essex', 'Staten Island', 'Bergen County'],
                'HR'  :['Norfolk', 'Portsmouth', 'Suffolk', 'Hampton',
                        'Virginia Beach', 'Newport News'],
                'NNJ' : [],
                'Charleston' : [],
                'Savannah' : [],
                'Miami' : [],
                }

# rate continunous, rate discrete, unc cont, unc discrete
DCT_NORMS = {'NYC': [TwoSlopeNorm(0, -5, 5), BoundaryNorm(np.arange(-5, 6, 1), 256),
                    Normalize(0.7, 1.1), BoundaryNorm(np.arange(0, 2.5, 0.5), 256)],
            'DC': [TwoSlopeNorm(0, -4, 4), BoundaryNorm(np.arange(-4, 5, 1), 256),
                    Normalize(0.65, 0.85), BoundaryNorm(np.arange(0, 1, 0.1), 256)],
                    # Normalize(0, 3), BoundaryNorm(np.arange(0, 1, 0.1), 256)],
            'NJ': [TwoSlopeNorm(0, -3, 3), BoundaryNorm(np.arange(-3, 3.25, 0.25), 256),
                    Normalize(0.65, 0.85), BoundaryNorm(np.arange(0, 1, 0.1), 256)],
                    # Normalize(0, 3), BoundaryNorm(np.arange(0, 1, 0.1), 256)],
            'HR': [TwoSlopeNorm(0, -4, 4), BoundaryNorm(np.arange(-4, 5, 1), 256),

                    Normalize(0, 2.0), BoundaryNorm(np.arange(0, 2.5, 0.5), 256)],
            'Philly': [TwoSlopeNorm(0, -3, 3), BoundaryNorm(np.arange(-3, 3.5, 0.5), 256),
                    Normalize(0, 2.0), BoundaryNorm(np.arange(0, 2.5, 0.5), 256)],
            'Kiritimati': [TwoSlopeNorm(0, -5, 5), BoundaryNorm(np.arange(-5, 6, 1), 256),
                    Normalize(0.3, 3.0), BoundaryNorm(np.arange(0, 3.5, 0.5), 256)],
            'Houston': [Normalize(-16, 0), BoundaryNorm(np.arange(-16, 1, 1), 256),
                        Normalize(0.5, 1.8), BoundaryNorm(np.arange(0.5, 2.0, 0.5), 256)],
                        # Normalize(0, 5), BoundaryNorm(np.arange(0, 6, 1), 256)], # tvar
            'Kennedy': [TwoSlopeNorm(0, -4, 4), BoundaryNorm(np.arange(-4, 5, 1), 256),
                    Normalize(0.7, 1.4), BoundaryNorm(np.arange(0, 1.5, 0.5), 256)],
            }

## placeholders
DCT_NORMS['Kiribati'] = DCT_NORMS['Kiritimati']

GS  = 13 # map tick size
CFS = 13 # colorbar fontlabel size
FS1 = 9 # GNSS station name size

class PlotExp(ExpBase):
    def __init__(self, exp, kind='rate', show_gps=True, continuous=True, sten=''):
        super().__init__(exp.dct_exp, exp.mp_exp0, exp.ref_sta, exp.neofs)
        self.show_gps     = show_gps
        self.continuous   = continuous
        self.sten = sten # not sure if this is used
        self.use_stitched = True if 'stitch' in exp.mp_exp0.lower() else False
        self.set_sten_paths(sten)

        # self.use_stitched = False
        self.WESN   = [*self.SNWE[2:], *self.SNWE[0:2]]
        self.df_gps = prep_gps(self.path_gps, self.reg, units='mm')
        self.kind   = self._set_kind(kind)
        self.da     = self.get_da()

        self.path_figs = op.join(PATH_RES, f'{exp.reg}_{datetime.today().year}')
        self._set_plot_parms()


    def _set_kind(self, kind='rate'):
        return kind.title()


    def _set_plot_parms(self):
        ## Plot Parms
        self.gs = GS
        self.fs1 = FS1
        self.cfs = CFS
        self.basemap = cimgt.GoogleTiles(
            url='https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}')
            # url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
            # url='https://server.arcgisonline.com/ArcGIS/rest/services/' \
            # 'Elevation/World_Hillshade_Dark/MapServer/tile/{z}/{y}/{x}.jpg')

        # self.basemap = cimgt.Stamen('toner-background')
        self.proj    = ccrs.PlateCarree()
        self.alpha   = 1.0
        self.zoom    = 13

        self.gs = GS
        self.fs1 = FS1
        self.gps_s = 100 # gps marker size
        self.gps_fs = 12 #FS1 # gps text size

        i = 0 if self.continuous else 1
        if self.kind == 'Rate':
            norm = DCT_NORMS.get(self.reg, DCT_NORMS['NYC'])[i]
            cmap = 'cmc.vik' if self.reg != 'Kiribati' else 'cmc.lajolla'
            self.pparms = {'cmap': cmap, 'norm': norm}
            self.cbar_lbl = 'Vertical Rate (mm/yr)'
            self.show_gps_names = True

        elif self.kind == 'Unc':
            norm = DCT_NORMS.get(self.reg, DCT_NORMS['NYC'])[i+2]
            # self.pparms = {'cmap': 'cmc.bilbao', 'norm': norm}
            # self.pparms = {'cmap': 'cmc.lajolla_r', 'norm': norm}
            self.pparms = {'cmap': 'afmhot_r', 'norm': norm}
            self.cbar_lbl = 'Vertical Rate Uncertainty (mm/yr)'

        elif self.kind == 'Tvar':
            norm = mpl.colors.Normalize(0, 6)
            self.pparms = {'cmap': 'cmo.amp', 'norm': norm}
            self.cbar_lbl = 'Temporal Variability (mm/yr)'

        else:
            raise Exception (f'Incorrect kind: {self.kind}. Use rate or unc')

        return


    def get_da(self, kind=None):
        kind = self.kind if kind == None else kind.title()
        """ Get the masked netcdf of the data to plot """
        path = self.path_rate_msk_nc if kind == 'Rate' else self.path_std_msk_nc
        units = '?'
        if kind == 'Tvar': # Houston, from Marin
            log.info ('Got temporal variability')
            path = Path(path).parent / 'geo_tstd_SR_recon10_masked.nc'
            self.cbar_lbl = 'Temporal Variability (mm/yr)'

        if self.use_stitched:
            path1 = self.path_rate_msk_nc_stitch if kind == 'Rate' else self.path_std_msk_nc_stitch
            if op.exists(path1):
                path = path1
            else:
                log.warning('Stitched experiment, but stitched da doesnt exist.')
        if 'ARIA' in self.dct_exp['root']:
            log.warning('Hacking aria to use the rate nc')
            path = op.join(op.dirname(op.dirname(path)), 'geo_rate_ARIA.nc')

        # get the netcdf
        ds = xr.open_dataset(path)
        da = ds['Band1'].rename(kind)
        da = da if 'mm' in da.units else da*1000
        # da *= 1000
        # log.critical('Forcing to mm')
        log.info('Data source: %s', path)
        return da


    def plot_basic(self, axes=None, da=None, grid_left=True, grid_bott=True):
        """ Plot rate or uncertainty in a single panel w/wout GPS """
        da = self.da.copy() if da is None else da
        if axes is None:
            fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
            axes.set_extent(self.WESN, crs=self.proj)
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, left=grid_left, bottom=grid_bott, size=self.gs)

        else:
            fig = axes.get_figure()

        axes.add_image(self.basemap, self.zoom)

        im   = axes.pcolormesh(da.lon, da.lat, da, shading='nearest', alpha=self.alpha,
                        transform=self.proj, **self.pparms)
        bbPlot.cartopy_cbar(im, ylabel=self.cbar_lbl, fontsize=self.cfs, pad=0.15)
        fig.set_label(f'Vup_{self.reg}_{self.kind}')

        ## add the GPS
        self.plot_gps(axes) if self.show_gps else None

        qH = 0.975 # 0.975 0.997
        qL = 1 - qH
        l, h = np.nanquantile(da, qL), np.nanquantile(da, qH)
        m    = np.nanmedian(da)
        # l, m, h = da.min(), da.mean(), da.max()
        l, m, h = da.mean()-2*da.std(), da.mean(), da.mean()+2*da.std()
        log.critical (f'InSAR Vup ({qH*100}%) {self.kind}: {m:.2f} [{l:.2f}, {h:.2f}]')
        # log.critical(f'InSAR Vup {self.kind}: {da.mean().item():.2f} +/- {da.std().item():.2f} mm/yr')

        interactive = False
        if interactive:
            import hvplot.xarray
            ## otherwise super slow (LaGuardia)
            # da1 = da.sel(lon=slice(-73.88460830568238, -73.8486756943176),
            #                     lat=slice(40.75743769431761, 40.79337030568239))
            del kind
            da1 = da.sel(lon=slice(-74,-73.96), lat=slice(40.6930, 41.7168))
            plot = da1.hvplot(cmap='coolwarm', clim=(0, 1.2))

                    #  coastline='10m',
            hvplot.show(plot)

        # axes.set_frame_on(False)
        # gl = axes.gridlines(draw_labels=False, linewidth=0)
        # bbPlot.fmt_gridlines(gl, left=False, bottom=False, size=0)

        return fig, axes


    def plot_gps(self, axes=None):
        """ Plot the MIDAS GPS locations (+ rates and names if axes=None) """
        col  = 'u_vel' if self.kind == 'Rate' else 'u_sig'
        msty = dict(marker='o', s=self.gps_s, linewidth=1.0, transform=self.proj,
                    alpha=1, zorder=30, edgecolors='w', c=self.df_gps[col],
                    **self.pparms)

        if axes is None:
            fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
            axes.add_image(self.basemap, self.zoom)
            s = axes.scatter(self.df_gps['lon'], self.df_gps['lat'], **msty)

            bbPlot.cartopy_cbar(s, 'V$_{up}$ (mm/yr)')
            axes.set_title(f'MIDAS GPS {self.reg}')
            fig.set_label(f'Midas_{self.reg}')
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, bottom=True, size=self.gs)

        else:
            fig = axes.get_figure()
            s   = axes.scatter(self.df_gps['lon'], self.df_gps['lat'], **msty)

        # Add a square reference pixel(s)
        ref_stas = self.get_ref_stas()
        ref_sta  = self.df_gps[self.df_gps.sta.isin(ref_stas)]

        msty    = {**msty, 'marker':'s', 'c':ref_sta[col], 's': self.gps_s+25, 'zorder':100}
        s = axes.scatter(ref_sta['lon'], ref_sta['lat'],  **msty)

        # dont show station name for overlapping stations
        dct_exl_plot = {'DC': 'USN7 USN9 WDC5 WDC6'.split()}
        # Add text x pixels of the station for rates only
        if self.kind == 'Rate' and self.show_gps_names:
            # sty_name = dict(facecolor='sandybrown', alpha=0.5, boxstyle='round')
            sty_name = dict(facecolor='lightgray', alpha=0.5, boxstyle='round')
            names = self.df_gps['sta'].astype('str').values
            names = [name for name in names if not name in dct_exl_plot.get(self.reg, [])]
            for i, name in enumerate(names):
                # y = 30 if name == 'NYBK' else -20
                y = 0#-10
                x = 10
                # halign = 'right' if name in 'NJI2 NYOB'.split() else 'left'
                halign = 'left'
                txt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')])

                geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=y, x=x)
                axes.text(self.df_gps['lon'].iloc[i], self.df_gps['lat'].iloc[i],
                        name.strip(), transform=text_transform, zorder=80,
                        verticalalignment='center', horizontalalignment=halign,
                        # size=self.gps_fs, color='lightgray', **txt_parms) # bbox=sty_name
                        size=self.gps_fs, color='white', **txt_parms) # bbox=sty_name

        print('Stations for MIDAS mean:', ', '.join(self.df_gps.sta.tolist()))
        qH = 0.997
        qL = 1 - qH
        m2 = self.df_gps[col].median()
        s2 = self.df_gps[col].std()
        l2, h2 = self.df_gps.u_vel.quantile(qL), self.df_gps.u_vel.quantile(qH)
        nn = f'{m2:.2f} +/- {s2:.2f}, [{l2:.2f}, {h2:.2f}] mm/yr'
        log.info(f'MIDAS ({qH*100}%) {self.kind.lower()} and spread: {nn}')

        axes.set_extent(self.WESN, crs=self.proj)
        return fig, axes


    def add_pts(self, axes, df, names=None, **plot_kws):
        """ Add scatter points (e.g., oil well point locations) """
        sty = dict(marker='v', s=100, facecolor='none', edgecolor='red', linewidth=1.25)
        sty.update(**plot_kws)
        try:
            lons, lats = df['lon'], df['lat']
        except:
            lons, lats = df.geometry.x, df.geometry.y

        s   = axes.scatter(lons, lats, transform=ccrs.PlateCarree(), **sty)
        if names is not None:
            sty_name = dict(facecolor='lightgray', alpha=0.5, boxstyle='round')
            for i, name in enumerate(names):
                # y = 30 if name == 'NYBK' else -20
                y = 0#-10
                x = 10
                # halign = 'right' if name in 'NJI2 NYOB'.split() else 'left'
                halign = 'left'
                txt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')])

                geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=y, x=x)
                axes.text(lons[i], lats[i],
                        name.strip(), transform=text_transform, zorder=80,
                        verticalalignment='center', horizontalalignment=halign,
                        size=self.fs1, color='lightgray', **txt_parms) # bbox=sty_name


    def plot_transect_map(self, transects, color='w', fontsize=16):
        """ Plot the Transects on the rate map with the GPS """
        import geopandas as gpd
        ## for setting grid labels manually
        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent(self.WESN, crs=self.proj)
        self.plot_basic(axes=axes, show_pts=False)

        for i, transect in enumerate(transects):
            path = self.path_wd / 'Transects' / f'{transect}.GeoJSON'
            assert path.exists(), f'Cannot find transect: {path}'
            gdf = gpd.read_file(path)
            gdf.plot(ax=axes, transform=self.proj, color=color, linewidth=3)

            gdf_pts = gdf.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
            df_pts = pd.DataFrame(gdf_pts.to_numpy()[0]).iloc[:, :2]
            # df_pts.columns = 'lon lat'.split()
            st, en = df_pts.iloc[0].to_numpy(), df_pts.iloc[-1].to_numpy()
            text_parms = dict(path_effects=[pe.withStroke(linewidth=3, foreground='w')],
                            color='k', transform=ccrs.PlateCarree(), fontsize=fontsize, zorder=500)

            letter     = string.ascii_uppercase[i]
            # if letter == 'A':
                # en[0] -= 0.005
                # en[1] -= 0.005

            if letter == 'B' and self.reg == 'NYC':
                text_parms['va'] = 'bottom'
                en[1] -= 0.005
                # st[1] -= 0.005

            axes.text(st[0], st[1], letter, **text_parms)
            axes.text(en[0], en[1], f"{letter}'", **text_parms)

        fig.set_label(f'{fig.get_label()}_transects')
        return fig, axes


    def plot_hist(self, ax=None):
        """ Histogram of rate/unc"""
        fig, axes  = plt.subplots(figsize=(8,8)) if ax is None else ax.get_figure(), ax
        m, l, h     = np.nanmedian(self.da), self.da.quantile(0.003), \
                            self.da.quantile(0.997)

        self.da.plot.hist(ax=axes, color='dimgray')
        axes.axvline(m, color='darkred', label=f'{m:.2f} +/- {(l+h)/2:.2f}')
        axes.axvline(l, linestyle='--', color='darkred', alpha=0.75)
        axes.axvline(h, linestyle='--',color='darkred', alpha=0.75)
        axes.legend()
        axes.set_title(f'{self.kind} (mm/yr)')

        # fig.suptitle(f'Correction Residuals (Raw - Corrected)', y=0.925)

        # print (f'2x Std: {resid.std():.2f}')
        # print (f'[{resid.quantile(0.025):.2f}, {resid.quantile(0.975):.2f}]')
        # print (f'{(resid.quantile(0.025) + resid.quantile(0.975))/2:.2f}')

        return fig, axes


    def plot_polys(self):
        """ Calculate average rate for hotspots """
        import geopandas as gpd
        epsg = 4326
        self.da.rio.write_crs(epsg, inplace=True)
        for f in os.listdir(self.path_crops):
            if f.endswith('GeoJSON') and not 'NYC' in f:
                path_crop = op.join(self.path_crops, f)
                gdf_crop     = gpd.read_file(path_crop)
                da_crop = self.da.rio.clip(gdf_crop.geometry, epsg, drop=True)
                # da_crop.plot()
                mu  = da_crop.mean()
                med = da_crop.median()
                log.info(f'{op.splitext(f)[0]} median {self.kind}: {med:.1f}')
        return


    def plot_city(self, city='Manhattan', loc=0, qH=0.997, da=None, crop=False,
                  show=True, alpha=0.0):
        """ Plot the rates with city boundaries and calculate the mean

        Manhattan, Brooklyn, Queens, Staten Island, NJ
        """
        import rioxarray as xrr
        import geopandas as gpd
        if isinstance(city, str):
            city_wgs, city_utm = city_shape(city, loc)
        else:
            assert isinstance(city, gpd.GeoDataFrame), 'city must be a GeoDataFrame'
            city_wgs = city
            city = city_wgs['CNTY_NM'].item() # hack for Texas

        da = self.da.copy() if da is None else da
        da.rio.set_spatial_dims(x_dim='lon', y_dim='lat', inplace=True)
        da.rio.write_crs(4326, inplace=True)

        da_reg = da.rio.clip(city_wgs.geometry, city_wgs.crs,  all_touched=False)
        W, S, E, N = da_reg.rio.bounds()


        if show:
            # plt.figure()
            sty  = dict(zorder=50, alpha=0.9, facecolor='none', edgecolor='white',
            linestyle=':', linewidth=1.0)
            pparms = {**self.pparms, 'shading':'nearest',
                      'transform':self.proj, 'zorder':11}
            fig, axes = plt.subplots(figsize=(10,10),
                                    subplot_kw={'projection': self.basemap.crs})
            WESN = [W-1e-3, E+1e-3, S-1e-3, N+1e-3] if crop else self.WESN
            axes.set_extent(WESN, crs=self.proj)
            axes.add_image(self.basemap, self.zoom)

            ## add the interpolated insar rates
            im   = axes.pcolormesh(da_reg.lon, da_reg.lat, da_reg, **pparms)
            bbPlot.cartopy_cbar(im, ylabel=self.cbar_lbl)

            # add the outline of the city
            axes.add_geometries(city_wgs.geometry, crs=self.proj, **sty)
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, bottom=True, size=self.gs)
            fig.set_label(f'Vup_{self.reg}_{self.kind}')
        else:
            fig, axes = None, None

        qL = 1-qH
        lbls = ['', 'Subsidence Only', 'Uplift Only']
        da_unc = self.get_rate_unc_nc()[1]
        for i, lbl in enumerate(lbls):
            if i == 0:
                dai = da_reg
                da_unci = da_unc
            elif i == 1:
                dai = da_reg.where(da_reg<0)
                da_unci = da_unc.where(da_reg<0)
            elif i == 2:
                dai = da_reg.where(da_reg>0)
                da_unci = da_unc.where(da_reg>0)

            l, h = np.nanquantile(dai, qL), np.nanquantile(dai, qH)
            m    = np.nanmedian(dai)

            l_unc, h_unc = np.nanquantile(da_unci, qL), np.nanquantile(da_unci, qH)
            m_unc = np.nanmedian(da_unci)

            # l, m, h = dai.min(), dai.mean(), dai.max()
            # l_unc, m_unc, h_unc = da_unci.min(), da_unci.mean(), da_unci.max()

            log.critical (f'{city} {lbl} ({100*qH}%) {self.kind}: {m:.2f} [{l:.2f}, {h:.2f}]')
            log.critical (f'{city} {lbl} ({100*qH}%) Uncertainty: {m_unc:.2f} [{l_unc:.2f}, {h_unc:.2f}]')

            # log.critical(f'InSAR Vup {self.kind}: {dai.mean().item():.2f} +/- {dai.std().item():.2f} mm/yr')
            # log.critical(f'InSAR Vup Unc: {da_unci.mean().item():.2f} +/- {da_unci.std().item():.2f} mm/yr')
            # print (f'{city} range: {h-l:.2f} mm/yr')

        return fig, axes


    def plot_rm_background(self, kind='median'):
        """ Remove the background trend and re-plot """
        kind = kind.lower()
        assert kind in 'mean median ice6g caron garner'.split(), 'Incorrect kind'
        if kind == 'mean':
            self.da -= self.da.mean()

        elif kind == 'median':
            self.da -= self.da.median()

        elif kind in 'ice6g caron'.split():
            from VLM.bzFRInGE.plotting.plotVLM import load_GIA
            ds_re = load_GIA(kind).load().interp_like(self.da, 'cubic')
            da_re = ds_re['VLM'] # VLM_unc
            self.da -= da_re
            kind = kind.upper() if kind == 'ice6g' else kind.title()

        elif kind == 'garner':
            self.da -= -1.5
            kind     = 'GIA'

        # self.da = self.da.where(np.abs(self.da)>0.3, np.nan)
        # self.pparms['norm'] = BoundaryNorm(np.arange(-5, 6, 2), ncolors=256)

        self.cbar_lbl = f'Vertical Rate-{kind} (mm/yr)'
        self.show_gps = False

        f, a = self.plot_basic()
        # a.set_title(kind)
        f.set_label(f'Vup_{self.reg}_{self.kind}-{kind}')


def get_tay22(reg, unc=False):
    """ Read the dataset from tay22 and format for use in plot_basic e.g."""
    import rioxarray as xrr
    if reg == 'NYC':
        src = op.join(PATH_RES, 'NYC_2023', 'Tay22.tif')
    elif reg == 'DC':
        f = 'Tay_Washington,D.C.tif' if not unc else 'Tay_Washington,D.C_std.tif'
        src = PATH_RES / 'DC_2024' / f
    else:
        raise Exception(f'{reg} not supported.')

    assert op.exists(src), 'Download from https://www.nature.com/articles/s41893-022-00947-z#Sec15'

    da_tay = 1000*xrr.open_rasterio(src, band_as_variable=True)['band_1'
            ].assign_attrs(units='mm/y', region=reg).rename(x='lon', y='lat')
    return da_tay


def plot_together(Exp1, show_gps=True, continuous=True):
    PObjR = PlotExp(Exp1, 'rate', show_gps, continuous)
    PObjU = PlotExp(Exp1, 'unc', show_gps, continuous)

    fig, axes = plt.subplots(figsize=(18,18), ncols=2,
                             subplot_kw={'projection': PObjR.basemap.crs})

    PObjR.plot_basic(axes[0])
    PObjU.plot_basic(axes[1])
    return


def make_numbers(Exp1, PlotClass, qH=0.997):
    PObjR = PlotClass(Exp1, kind='rate')
    PObjU = PlotClass(Exp1, kind='unc')
    qL = 1-qH
    attrs = readfile.read_attribute(PObjR.path_vup_geo)
    log.critical(f'Experiment: {PObjR.mp_exp} (neofs={PObjR.neofs})')
    log.critical('%s to %s', attrs['START_DATE'], attrs['END_DATE'])

    l, h = np.nanquantile(PObjR.da, qL), np.nanquantile(PObjR.da, qH)
    m    = np.nanmedian(PObjR.da)

    l1, h1 = np.nanquantile(PObjU.da, qL), np.nanquantile(PObjU.da, qH)
    m1   = np.nanmedian(PObjU.da)

    ## overall InSAR mean
    log.info (f'InSAR Vup {100*qH}%: {m:.2f} [{l:.2f}, {h:.2f}]')
    log.info (f'InSAR Unc {100*qH}%: {m1:.2f} [{l1:.2f}, {h1:.2f}]')

    l, m, h = PObjR.da.min(), PObjR.da.mean(), PObjR.da.max()
    l1, m1, h1 = PObjU.da.min(), PObjU.da.mean(), PObjU.da.max()

    log.info (f'InSAR Vup 100%: {m:.2f} [{l:.2f}, {h:.2f}]')
    log.info (f'InSAR Unc 100%: {m1:.2f} [{l1:.2f}, {h1:.2f}]')


    ## overall GPS mean
    print ('Stations for mean:', ', '.join(PObjR.df_gps.sta.tolist()))
    m2     = PObjR.df_gps.u_vel.median()
    l2, h2 = PObjR.df_gps.u_vel.quantile(qL), PObjR.df_gps.u_vel.quantile(qH)

    m3     = PObjR.df_gps.u_sig.median()
    l3, h3 = PObjR.df_gps.u_sig.quantile(qL), PObjR.df_gps.u_sig.quantile(qH)

    log.info(f'MIDAS rate {100*qH}%: {m2:.2f}, [{l2:.2f}, {h2:.2f}] mm/yr')
    log.info(f'MIDAS unc {100*qH}%: {m3:.2f} +/- [{l3:.2f}, {h3:.2f}] mm/yr')

    l2, m2, h2 = PObjR.df_gps.u_vel.min(), PObjR.df_gps.u_vel.mean(), PObjR.df_gps.u_vel.max()
    l3, m3, h3 = PObjR.df_gps.u_sig.min(), PObjR.df_gps.u_sig.mean(), PObjR.df_gps.u_sig.max()
    log.info(f'MIDAS rate 100%: {m2:.2f}, [{l2:.2f}, {h2:.2f}] mm/yr')
    log.info(f'MIDAS unc 100%: {m3:.2f} +/- [{l3:.2f}, {h3:.2f}] mm/yr')

    ## percent significant
    # get number of nonnan pix, aka the actual rates we are considering
    n_pix  = PObjR.da.where(np.isnan(PObjR.da), 1).sum()

    # get a count of where pixels are >= sigma
    da_sig  = xr.where(np.abs(PObjR.da)>=PObjU.da, 1, 0)
    per_sig = 100*(da_sig.sum() / n_pix)
    log.info (f'Percent significant at 1sig: {per_sig:.1f}%')


    ## RMSE
    # PObjR.plot_rmse(show=False)

    ## cities
    for city in DCT_CITIES[Exp1.reg]:
        loc = 4 if city == 'Essex' else 0
        PObjR.plot_city(city, loc, show=False)
        # PObjU.plot_city(city, loc, show=False)
    return



if __name__ == '__main__':
    # Exp0      = ExpBase(Kiribati_SR, 'ERA5_SET_ex_Fast', neofs=5)
    # Exp0      = ExpBase(Kennedy_SR, 'Base_ex_Fast', neofs=15)
    Exp0      = ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', neofs=10)
    PlotExp(Exp0, 'rate', show_gps=True).plot_basic()
    # PlotExp(Exp0, 'unc', show_gps=False).plot_basic()

    # plot_together(Exp0)

    # GPSObj = analyzeGPS.GPS_Rate_Analysis(Exp0)
    # GPSObj.plot_gps_vs_insar()

    path_figs = op.join(PATH_RES, f'{Exp0.reg}_2024')
    # bbPlot.savefigs(path_figs, False, True, dpi=300)

    plt.show()