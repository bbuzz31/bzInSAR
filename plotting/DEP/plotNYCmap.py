"""
Plot NYC Maps in Python for Paper
"""
import string
import xarray as xr
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
from cartopy.io import img_tiles as cimgt

from mintpy.utils import readfile
from BZ import bbLogger, bbGIS, bbPlot
from BZ.bbGIS import city_shape

from VLM.bzFRInGE.contrib import *
from VLM.bzFRInGE.contrib.experiments import *
from VLM.bzFRInGE.contrib.FRINGEBase import ExpBase

## for inset plots
from VLM.bzFRInGE.contrib.plotting.plotTS import *

DCT_TRANSECTS = {'NYC': 'NJ1 Manhattan1 Brooklyn1'.split()}

GS = 13 # map tick size
FS1 = 9 # GNSS station name size

# plt.switch_backend('Qt5Agg')

## plot parms
class PlotNYC(ExpBase):
    def __init__(self, exp, kind='rate', show_gps=True, continuous=True):
        super().__init__(exp.dct_exp, exp.mp_exp0, exp.ref_sta, exp.neofs)
        self.show_gps     = show_gps
        self.continuous   = continuous
        self.use_stitched = True if 'stitch' in exp.mp_exp0.lower() else False
        self.WESN   = [*self.SNWE[2:], *self.SNWE[0:2]]
        self.df_gps = prep_gps(self.path_gps, self.reg, units='mm')
        self.kind   = self._set_kind(kind)
        self.da     = self.get_da()
        self._set_plot_parms()


    def _set_kind(self, kind='rate'):
        return kind.title()


    def _set_plot_parms(self):
        ## Plot Parms
        self.basemap = cimgt.GoogleTiles(
            url='https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}')
        self.proj    = ccrs.PlateCarree()
        self.alpha   = 1.0
        self.zoom    = 13

        if self.kind == 'Rate':
            norm = mpl.colors.TwoSlopeNorm(0, -5, 5) if self.continuous else \
                mpl.colors.BoundaryNorm(np.arange(-5, 6, 1), 256)
            self.pparms = {'cmap': 'cmc.vik', 'norm': norm}
            # self.cbar_lbl = 'V$_{Up}$ (mm/yr)'
            self.cbar_lbl = 'Vertical Rate (mm/yr)'
        elif self.kind == 'Unc':
            norm = mpl.colors.Normalize(0, 2.0) if self.continuous else \
                mpl.colors.BoundaryNorm(np.arange(0, 2.5, 0.5), 256)
            self.pparms = {'cmap': 'cmc.lajolla', 'norm': norm}
            self.cbar_lbl = 'Vertical Rate Uncertainty (mm/yr)'

        else:
            raise Exception ('Incorrect kind: {self.kind}. Use rate or unc')

        return


    def get_da(self, kind=None):
        kind = self.kind if kind == None else kind
        """ Get the masked netcdf of the data to plot """
        path = self.path_rate_msk_nc if kind == 'Rate' else self.path_std_msk_nc
        if self.use_stitched:
            path1 = self.path_rate_msk_nc_stitch if kind == 'Rate' else self.path_std_msk_nc_stitch
            if op.exists(path1):
                path = path1
            else:
                log.warning('Stitched experiment, but stitched da doesnt exist.')

        # get the netcdf
        ds = xr.open_dataset(path)
        da = ds['Band1'].rename(kind)
        da = da*1000 if 'm' in da.units else da
        return da


    def plot_basic(self, axes=None):
        """ Plot rate or uncertainty in a single panel w/wout GPS """
            # path = self.path_rate_nc if self.kind == 'Rate' else self.path_std_nc
        da = self.da.copy()
        # ds1 = xr.open_dataset('/Users/buzzanga/Downloads/dataverse_files/NewYork_TA033/NewYork_TA033_velocity_InSAR.nc')
        # da = ds1['Band1']*1000
        # breakpoint()

        # da = da.where(da<=-5, np.nan)
        # print (da.size - np.isnan(da).sum())

        if axes is None:
            fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
            axes.set_extent(self.WESN, crs=self.proj)
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, bottom=True, size=GS)
        else:
            fig = axes.get_figure()

        axes.add_image(self.basemap, self.zoom)

        im   = axes.pcolormesh(da.lon, da.lat, da, shading='nearest', alpha=self.alpha,
                        transform=self.proj, **self.pparms)

        bbPlot.cartopy_cbar(im, ylabel=self.cbar_lbl)
        fig.set_label(f'Vup_{self.reg}_{self.kind}')

        ## add the GPS
        if self.show_gps:
            self.plot_gps(axes)


        # print(f'Vup {self.kind}: {da.mean().item():.2f} +/- {da.std().item():.2f} mm/yr')
        l, h = np.nanquantile(da, 0.003), np.nanquantile(da, 0.997)
        m    = np.nanmedian(da)
        log.info (f'InSAR Vup {self.kind}: {m:.2f} [{l:.2f}, {h:.2f}]')

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


        return fig, axes


    def plot_basic_box(self, pts='LaGuardia Williams'.split()):
        """ Plot basic plus a 2000 m bbox around point to show inset locations """
        fig, axes  = self.plot_basic()

        for pt in pts:
            lalo = dct_pts[pt][0] #  use plot_ts to make a box around the point
            SNWE = buffer_point(*lalo, 2000)

            sty  = dict(zorder=50, alpha=0.9, facecolor='none', edgecolor='white',
                        linestyle='--', linewidth=1.5)
            bbox = bbGIS.bbox2poly(SNWE=SNWE)
            axes.add_geometries([bbox], crs=self.proj, **sty)

            ## add the point in the middle
            pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                            color='w', transform=ccrs.PlateCarree(), s=30, marker='^')
            axes.scatter(lalo[1], lalo[0], **pt_parms)

        ## add city borders and names
        # 'Bronx': (-73.8145, 40.8055), 'Bergen County': (-74.13, 40.8),  'Staten Island': (-74.111, 40.527)
        dct_names = {'Bronx': (-73.8145, 40.7950), 'Brooklyn': (-73.9762, 40.625),
                     'Bergen County': (-74.13, 40.79), 'Hudson': (-74.069, 40.650),
                    'Manhattan': (-73.990, 40.770), 'Queens': (-73.8634, 40.695),
                    'Staten Island': (-74.111, 40.5155), 'Essex': (-74.19036, 40.780)}

        cities = 'Brooklyn Bronx Queens Manhattan Hudson Essex'.split()
        for city in cities+['Staten Island', 'Bergen County']:
            loc = 4 if city == 'Essex' else 0
            city_wgs, city_utm = city_shape(city, loc)
            sty  = dict(zorder=50, alpha=0.9, facecolor='none', edgecolor='white',
                        linestyle=':', linewidth=1.0)
            axes.add_geometries(city_wgs.geometry, crs=self.proj, **sty)

            if self.kind == 'Rate':
                # y = 30 if name == 'NYBK' else -20
                sty_name = dict(facecolor='lightgray', alpha=0.5, boxstyle='round')
                y, halign = 0, 'left'
                rot = 61 if city == 'Manhattan' else 0
                rot = 65 if city == 'Hudson' else rot
                x, y, txt = *dct_names[city], city.replace('County', '')
                geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=y)
                axes.text(x, y, txt, transform=text_transform, zorder=20, rotation=rot,
                        verticalalignment='center', horizontalalignment=halign,
                        bbox=sty_name, size=11.5, color='black')

        ## add the Battery TG
        lalo = (40.7010692,  -74.0143231)
        # pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='w')],
        pt_parms = dict(edgecolors='white', facecolor='none',
                        transform=ccrs.PlateCarree(), s=70, marker='v')
        axes.scatter(lalo[1], lalo[0], **pt_parms)


        fig.set_label(f'{fig.get_label()}_{"_".join(pts)}')
        return


    def plot_gps(self, axes=None):
        """ Plot the MIDAS GPS locations (+ rates and names if axes=None) """
        col  = 'u_vel' if self.kind == 'Rate' else 'u_sig'
        msty = dict(marker='o', s=50, linewidth=1.5, transform=self.proj,
                    alpha=1, zorder=30, edgecolors='black', c=self.df_gps[col],
                    **self.pparms)

        if axes is None:
            fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
            axes.add_image(self.basemap, self.zoom)
            s = axes.scatter(self.df_gps['lon'], self.df_gps['lat'], **msty)

            bbPlot.cartopy_cbar(s, 'V$_{up}$ (mm/yr)')
            axes.set_title(f'MIDAS GPS {self.reg}')
            fig.set_label(f'Midas_{self.reg}')
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, bottom=True, size=GS)

        else:
            fig = axes.get_figure()
            s   = axes.scatter(self.df_gps['lon'], self.df_gps['lat'], **msty)

        # Add a square reference pixel(s)
        ref_stas = self.get_ref_stas()
        ref_sta  = self.df_gps[self.df_gps.sta.isin(ref_stas)]
        msty    = {**msty, 'marker':'s', 'c':ref_sta[col], 's': 75}
        s = axes.scatter(ref_sta['lon'], ref_sta['lat'],  **msty)

        # Add text x pixels of the station for rates only
        if self.kind == 'Rate':
            # sty_name = dict(facecolor='sandybrown', alpha=0.5, boxstyle='round')
            sty_name = dict(facecolor='lightgray', alpha=0.5, boxstyle='round')
            names = self.df_gps['sta'].astype('str').values
            for i, name in enumerate(names):
                # y = 30 if name == 'NYBK' else -20
                y = -10
                # halign = 'right' if name in 'NJI2 NYOB'.split() else 'left'
                halign = 'left'
                txt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')])

                geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=y)
                axes.text(self.df_gps['lon'].iloc[i], self.df_gps['lat'].iloc[i],
                        name.strip(), transform=text_transform, zorder=80,
                        verticalalignment='center', horizontalalignment=halign,
                        size=FS1, color='lightgray', **txt_parms) # bbox=sty_name

        log.info(
            f'MIDAS mu: {self.df_gps.u_vel.mean():.2f} +/- {self.df_gps.u_sig.mean():.2f} mm/yr')

        axes.set_extent(self.WESN, crs=self.proj)
        return fig, axes


    def plot_transect(self):
        """ Plot the Transects on the rate map with the GPS """
        import geopandas as gpd
        import string
        fig, axes  = self.plot_basic()
        for i, transect in enumerate(DCT_TRANSECTS[self.reg]):
            path = op.join(self.path_wd, 'Transects', f'Transect_{transect}.GeoJSON')
            gdf = gpd.read_file(path)
            gdf.plot(ax=axes, transform=self.proj, color='w')
            gdf_pts = gdf.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
            df_pts = pd.DataFrame(gdf_pts.to_numpy()[0]).iloc[:, :2]
            # df_pts.columns = 'lon lat'.split()
            st, en = df_pts.iloc[0].to_numpy(), df_pts.iloc[-1].to_numpy()
            text_parms = dict(path_effects=[pe.withStroke(linewidth=3, foreground='w')],
                            color='k', transform=ccrs.PlateCarree(), fontsize=16)

            letter     = string.ascii_uppercase[i]
            # if letter == 'A':
                # en[0] -= 0.005
                # en[1] -= 0.005

            if letter == 'B':
                text_parms['va'] = 'bottom'
                en[1] -= 0.005
                # st[1] -= 0.005

            axes.text(st[0], st[1], letter, **text_parms)
            axes.text(en[0], en[1], f"{letter}'", **text_parms)

        ## add the Battery TG
        lalo = (40.7010692,  -74.0143231)
        pt_parms = dict(edgecolors='white', facecolor='none',
                        transform=ccrs.PlateCarree(), s=70, marker='v')
        axes.scatter(lalo[1], lalo[0], **pt_parms)

        # axes = bbPlot._scale_bar(axes, 10, (.680, 0.9625))
        axes = bbPlot._scale_bar(axes, 10, (0.85, 0.015))
        # axes.grid(False)
        fig.set_label(f'{fig.get_label()}_transects')
        return


    def plot_hist(self):
        """ Histogram of rate/unc"""
        fig, axes   = plt.subplots(figsize=(8,8))
        mu, sig = self.da.mean(), self.da.std()

        self.da.plot.hist(ax=axes, color='dimgray')
        axes.axvline(mu, color='darkred', label=f'{mu:.2f} +/- {sig:.2f} (mm/yr)')
        ax.axvline(mu-sig, linestyle='--', color='darkred', alpha=0.75)
        ax.axvline(mu+sig, linestyle='--',color='darkred', alpha=0.75)
        ax.legend()

        fig.suptitle(f'Correction Residuals (Raw - Corrected)', y=0.925)

        # print (f'2x Std: {resid.std():.2f}')
        # print (f'[{resid.quantile(0.025):.2f}, {resid.quantile(0.975):.2f}]')
        # print (f'{(resid.quantile(0.025) + resid.quantile(0.975))/2:.2f}')

        return


    def plot_insets_col(self, pts='LaGuardia Williams'.split(), npix=3):
        """ Plot a close up of two locations, column wise """
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig  = plt.figure(figsize=(13, 13))
        gs   = GridSpec(nrows=2, ncols=2, figure=fig, wspace=0.25, hspace=0.01, height_ratios=[2.5,1])

        ## plot the map with s6, j3, tg
        axeMap0 = plt.subplot(gs[0, 0], projection=self.basemap.crs)
        axeMap1 = plt.subplot(gs[0, 1], projection=self.basemap.crs)
        axeMaps = [axeMap0, axeMap1]

        axeTS0    = fig.add_subplot(gs[1, 0])
        axeTS1    = fig.add_subplot(gs[1, 1])
        axeTSs    = [axeTS0, axeTS1]

        for i, pt in enumerate(pts):
            lalo = dct_pts[pt][0] #  use plot_ts to make a box around the point
            S, N, W, E = buffer_point(*lalo, 2000)

            self.zoom  = 15
            axes  = self.plot_basic(axes=axeMaps[i])[1]
            im    = axes.get_children()[0]
            im.set_alpha(0.75)
            axes.set_extent([W, E, S, N])

            gl = axes.gridlines(draw_labels=True)
            gparms = dict(left=False, right=True, size=GS-1) if i == 1 else {}
            bbPlot.fmt_gridlines(gl, bottom=True, **gparms)

            axeTS = plot_ts_vup1_dd(self, pts[i], f'{pts[i]}_Stable', npix=npix, axes=axeTSs[i])[1]
            axeTS.set_title('')
            if i == 1:
                axeTS.set_ylabel('')

            # axeTSs[i].plot(np.random.rand(180))
        # axeTS0.plot(np.random.rand(180))
        # axeTS1    = fig.add_subplot(gs[1, 1])#, sharey=axeTS0)
        # axeTS1.plot(np.random.rand(180))

        cbar_ax = fig.add_axes([0.5, 0.377, 0.02, 0.457]) # left, bot, wid, hgt
        cbar_ax.set_title('V$_{Up}$')
        cbar_ax.set_xlabel('(mm yr$^{-1}$)')
        cbar    = plt.colorbar(im, cax=cbar_ax,  spacing='uniform')


    def plot_inset(self, pt='LaGuardia', npix=3):
        """ Plot one point, row wise """
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig = plt.figure(figsize=(18, 10), constrained_layout=True)
        gs  = GridSpec(2, 3, width_ratios=[1.3, 1, 1], figure=fig)
        axes = []
        axes.append(plt.subplot(gs[0, 0], projection=self.basemap.crs))
        axes.append(plt.subplot(gs[0, 1:]))

        lalo = dct_pts[pt][0] #  use plot_ts to make a box around the point
        S, N, W, E = buffer_point(*lalo, 2000)

        self.zoom  = 15

        ax  = self.plot_basic(axes=axes[0])[1]
        im  = ax.get_children()[0]
        im.set_alpha(0.75)
        ax.set_extent([W, E, S, N])
        gl = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 190, 0.01),
                    ylocs=np.arange(-90, 100, 0.01))

        gparms = dict(left=False, right=True)# if i == 1 else {}
        bbPlot.fmt_gridlines(gl, bottom=True, size=GS-1)#, **gparms)

        axeTS = plot_ts_vup1_dd(self, pt, f'{pt}_Stable', npix=npix, axes=axes[1])[1]
        # axes[1].plot(np.random.rand(180))
        axeTS.set_title('')
        axeTS.set_ylim([-40, 20]) if pt == 'LaGuardia' else ''

        ## make colorbar
        divider = make_axes_locatable(axes[0])
        cbar_ax = divider.new_horizontal(size="5%", pad=0.85, axes_class=plt.Axes, pack_start=True)
        fig.add_axes(cbar_ax)

        cbar    = fig.colorbar(im, cbar_ax, label='Vertical Rate (mm/yr)')
        # cbar_ax.set_xlabel('mm/yr', labelpad=10)
        cbar.ax.yaxis.set_label_position('left')
        cbar_ax.yaxis.set_ticks_position('left')
        # leg = axes[1].legend(loc='lower left', frameon=True, shadow=True, ncol=1, prop={'size':'8'})


    def plot_insets(self, pts='LaGuardia Williams'.split()):
        """ Plot two points, row wise """
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        self.show_gps = False
        fig = plt.figure(figsize=(14, 9), constrained_layout=False)
        gs  = GridSpec(2, 3, width_ratios=[1, 1, 1], hspace=0.33, wspace=0.0, figure=fig)
        axeMaps = []
        axeTSs  = []
        axeMaps.append(plt.subplot(gs[0, 0], projection=self.basemap.crs))
        axeMaps.append(plt.subplot(gs[1, 0], projection=self.basemap.crs))
        axeTSs.append(plt.subplot(gs[0, 1:]))
        axeTSs.append(plt.subplot(gs[1, 1:]))
        dct_npix = {'LaGuardia': 3, 'Williams': 3, 'Ashe': 3, 'Woodside': 3}
        for i, pt in enumerate(pts):
            lalo = dct_pts[pt][0] #  use plot_ts to make a box around the point
            npix = dct_npix[pt]
            S, N, W, E = buffer_point(*lalo, 2000)

            self.zoom  = 15

            ax  = self.plot_basic(axes=axeMaps[i])[1]
            # ax.text(-0.01, 1.01, string.ascii_uppercase[i], transform=ax.transAxes,
            ax.text(1.0, 0.999, string.ascii_uppercase[i], transform=ax.transAxes,
                    verticalalignment='top', fontsize=20, color='w',
                    bbox=dict(boxstyle='square,pad=0.05', facecolor='black'))
            im  = ax.get_children()[0]
            im.set_alpha(0.75)
            im.colorbar.remove()
            ax.set_extent([W, E, S, N])
            # bbPlot._scale_bar(ax, 1, (0.85, 0.015))
            bbPlot._scale_bar(ax, 1, (0.20, 0.015))

            pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                            color='w', transform=ccrs.PlateCarree(), s=30)
            ax.scatter(lalo[1], lalo[0], **pt_parms)

            gl = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 190, 0.01),
                        ylocs=np.arange(-90, 100, 0.01))

            bbPlot.fmt_gridlines(gl, left=True, right=False, bottom=True, size=GS-1)

            if pt == 'Woodside':
                axeTS = plot_woodside_dd(self, axes=axeTSs[i])[1]

            else:
                axeTS = plot_ts_vup1_dd(self, pt, f'{pt}_Stable', npix=npix, axes=axeTSs[i])[1]
                # axeTSs[i].plot(np.random.rand(180))
                # axeTS = axeTSs[i]

            # pos = axeTS.get_position().bounds
            # pos[0] = pos[0]-0.1
            # axeTS.set_position(mpl.transforms.Bbox.from_bounds(*pos))

            axeTS.set_title('')
            # axeTS.set_ylim([-40, 10]) if pt == 'LaGuardia' else ''

            axeTS.tick_params(labelbottom=False) if i == 0 else ''
            axeTS.tick_params(axis='both', which='major', labelsize=13)

            if i == 1:
                axeTS.set_ylabel(axeTS.get_ylabel(), labelpad=12)

            # ## make colorbars
            # divider = make_axes_locatable(axeMaps[i])
            # cbar_ax = divider.new_horizontal(size="5%", pad=0.55, axes_class=plt.Axes, pack_start=True)
            # fig.add_axes(cbar_ax)


            # cbar    = fig.colorbar(im, cbar_ax, label='Vertical Rate (mm/yr)')
            # cbar.ax.yaxis.set_label_position('left')
            # cbar_ax.yaxis.set_ticks_position('left')

        # make a single colorbar across
        cbar_ax = fig.add_axes([0.152, 0.475, 0.160, 0.0160]) # left, bot, wid, hgt
        cbar_ax.set_title('V$_{Up}$  (mm yr$^{-1}$)', fontsize=11)#, labelpad=1v#, rotation=270, )
        # cbar_ax.yaxis.set_label_position('right')
        cbar    = plt.colorbar(im, cax=cbar_ax, orientation='horizontal', spacing='uniform')
        # gs.update(left=0.07, right=0.23)#, bottom=0.1, top=0.9)
        # gs.update(hspace=-0.5)#, wspace=10.00)
        fig.set_label(f'Insets_{"_".join(pts)}')
        fig.constrained_layout=True
        return


    def plot_cities(self, city='NJ', loc=0, show=True):
        """ Plot the rates with city boundaries and calculate the mean

        Manhattan, Brooklyn, Queens, Staten Island, NJ
        """
        import rioxarray as xrr
        from shapely.geometry import mapping
        city_wgs, city_utm = city_shape(city, loc)

        da = self.da.copy()
        da.rio.set_spatial_dims(x_dim='lon', y_dim='lat', inplace=True)
        da.rio.write_crs(4326, inplace=True)

        da_reg = da.rio.clip(city_wgs.geometry.apply(mapping))

        if show:
            plt.figure()
            fig, axes = plt.subplots(figsize=(10,10),
                                    subplot_kw={'projection': self.basemap.crs})
            axes.set_extent(self.WESN, crs=self.proj)
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, bottom=True, size=GS)

            axes.add_image(self.basemap, self.zoom)

            ## add the interpolated insar rates
            im   = axes.pcolormesh(da_reg.lon, da_reg.lat, da_reg, zorder=20 **self.pparms)

            bbPlot.cartopy_cbar(im, ylabel=self.cbar_lbl)
            fig.set_label(f'Vup_{self.reg}_{self.kind}')

        l, h = np.nanquantile(da_reg, 0.003), np.nanquantile(da_reg, 0.997)
        m    = np.nanmedian(da_reg)

        log.info (f'{city} Vup {self.kind}: {m:.1f} [{l:.1f}, {h:.1f}]')
        print (f'{city} range: {h-l:.2f} mm/yr')
        return


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

    ## inherit this from analyzeGPS
    def plot_rmse(self, npix=0, show=True):
        print ('NO')
        return
        excl    = DCT_EXCL[self.reg]
        if npix == 0:
            excl += [self.ref_sta]

        df_rmse = gps_rmse_npix(self, npix=npix, exclude=excl, verbose=False)[0]


        print ('Stations for RMSE:', ', '.join(df_rmse.index.tolist()))
        log.info(f'RMSE: {1000*df_rmse.abs().mean().item():.1f} mm/yr')

        if show:
            ## plot the RMSE and actual resid for the two comparisons
            fig, axes = plt.subplots(figsize=(10,4))
            resid = df_rmse.resid*1000
            mu = resid.abs().mean()
            axes.scatter(df_rmse.index, resid, color='k')
            axes.plot(df_rmse.index, np.tile(mu, len(df_rmse.index)),
                    color='k', label=f'RMSE: {mu:.1f} (mm/yr)', linestyle='--')

            axes.grid(color='gray', linestyle = '--', linewidth=0.1)
            axes.legend(ncol=1)
            axes.set_ylabel('Residual (mm/yr)')
            axes.set_xlabel('GNSS Station')

            if (resid < 1).all() and (resid>0).all():
                axes.set_ylim([0,1])

            fig.set_label(f'{self.mp_exp0}_RMSE')
        return


def plot_box_ref_sta():
    import seaborn as sns
    """ Plot smoothed histograms of different ref_sta experiments """
    Exp0      = ExpBase(NYC_SRc, 'ERA5_SET_PM_Fast', 'NJHT', neofs=20)
    Exp1      = ExpBase(NYC_SRc, 'ERA5_SET_PM_Fast', 'NYBK', neofs=20)
    Exp2      = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_Fast', 'NYBK', neofs=20)

    lst_sers = []
    lbls     = 'NJHT NYBK NJHT+NYBK'.split()
    for i, Expi in enumerate([Exp0, Exp1, Exp2]):
        path = Expi.path_vup_geo_stitch if 'stitch' in \
            Expi.mp_exp0.lower() else Expi.path_vup_geo
        arr0 = readfile.read(path)[0] * 1000
        ser  = pd.Series(arr0.reshape(-1)).dropna().rename(lbls[i])
        lst_sers.append(ser)

    df = pd.concat(lst_sers, axis=1)
    ## drop any point with a high value
    df = df[~(df.abs()>=5).any(axis=1)]

    ax = sns.boxplot(data=df, color='dimgray', whis=[1.0, 99.0], fliersize=0)
    ax.set_ylim([-4, 2])
    ax.grid(color='gray', linestyle = '--', linewidth=0.1)
    plt.show()
    # breakpoint()


def plot_hist_ref_sta():
    import seaborn as sns
    """ Plot smoothed histograms of different ref_sta experiments """
    Exp0      = ExpBase(NYC_SRc, 'ERA5_SET_PM_ex_Fast', 'NJHT', neofs=20)
    Exp1      = ExpBase(NYC_SRc, 'ERA5_SET_PM_ex_Fast', 'NYBK', neofs=20)
    Exp2      = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20)

    lbls     = 'NJHT NYBK NJHT+NYBK'.split()
    fig, axes = plt.subplots(figsize=(13, 5), ncols=3, sharey=True)
    for i, Expi in enumerate([Exp0, Exp1, Exp2]):
        path = Expi.path_rate_msk_nc_stitch if 'stitch' in \
            Expi.mp_exp0.lower() else Expi.path_rate_msk_nc
        arr0 = xr.open_dataarray(path) * 1000
        ser  = arr0.to_dataframe().reset_index()['Band1'].dropna().rename(lbls[i])
        ser  = ser[ser.abs()<=5]
        ax = sns.kdeplot(data=ser, fill=True, color='k', ax=axes[i], legend=False)
        # ax.set_ylim([-4, 2])
        ## add the mean and standard deviation
        mu, l, h = ser.median(), ser.quantile(0.003), ser.quantile(0.997)
        lh = np.mean([l,h])
        ax.axvline(mu, ymax=0.65, color='r',
                   linestyle='--', label=f'{mu:.1f} $\pm$ {lh:.1f}')

        ax.axvspan(mu-lh, mu+lh, ymax=0.65, color='darkred', alpha=0.2)

        # ax.axvline(mu-sig, ymax=0.65, color='darkred', linestyle=':')
        # ax.axvline(mu+sig, ymax=0.65, color='darkred', linestyle=':')
        ax.grid(color='gray', linestyle = '--', linewidth=0.1)
        ax.set_xlabel(lbls[i])

        ax.legend()

    fig.set_label('Exp_Hists')
    # breakpoint()


def plot_hist_corr():
    import seaborn as sns
    """ Plot smoothed histograms of correction vs none experiments """
    Exp1      = ExpBase(NYC_SRc, 'Base_Stitch_ex_Fast', 'NYBK', neofs=0)
    Exp2      = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20)

    rate_diff = xr.open_dataarray(Exp1.path_rate_msk_nc_stitch) - \
                xr.open_dataarray(Exp2.path_rate_msk_nc_stitch)
    unc_diff = xr.open_dataarray(Exp1.path_std_msk_nc_stitch) - \
                xr.open_dataarray(Exp2.path_std_msk_nc_stitch)

    lbls     = 'Rate Uncertainty'.split()
    fig, axes = plt.subplots(figsize=(13, 5), ncols=2)
    for ax, arr, lbl in zip(axes, [rate_diff, unc_diff], lbls):
        ser  = 1000*arr.to_dataframe().reset_index()['Band1'].dropna().rename(lbl)
        # ser  = ser[ser.abs()<=5]
        cut   = 0 if lbl == 'Uncertainy' else 3
        ax = sns.kdeplot(data=ser, fill=True, color='k', ax=ax, legend=False, cut=cut)
        # ax.set_ylim([-4, 2])
        ## add the mean and standard deviation
        mu, l, h = ser.median(), ser.quantile(0.003), ser.quantile(0.997)
        print (f'{mu:.2f}, [{l:.2f}, {h:.2f}]')
        print (ser.max())

        lh = np.mean([l,h])
        ax.axvline(mu, ymax=0.65, color='r',
                   linestyle='--', label=f'{mu:.1f} $\pm$ {lh:.1f}')

        ax.axvspan(0, mu+lh, ymax=0.65, color='darkred', alpha=0.2)

        # ax.axvline(mu-sig, ymax=0.65, color='darkred', linestyle=':')
        # ax.axvline(mu+sig, ymax=0.65, color='darkred', linestyle=':')
        ax.grid(color='gray', linestyle = '--', linewidth=0.1)
        ax.set_xlabel(lbl)

        ax.legend()

    fig.set_label('Corr_Hists')


def plot_together(Exp1, show_gps=True, continuous=True):
    PObjR = PlotNYC(Exp1, 'rate', show_gps, continuous)
    PObjU = PlotNYC(Exp1, 'unc', show_gps, continuous)

    fig, axes = plt.subplots(figsize=(18,18), ncols=2,
                             subplot_kw={'projection': PObjR.basemap.crs})

    PObjR.plot_basic(axes[0])
    PObjU.plot_basic(axes[1])
    return


## ----------- Wrappers
def make_main(Exp1):
    PObj = PlotNYC(Exp1, kind='rate')
    PObj.plot_basic_box(pts='LaGuardia Williams'.split())
    # PObj.plot_insets(pts='LaGuardia Williams'.split())
    # PObj.plot_cities('Essex', 4)


def make_supp(Exp1):
    from VLM.bzFRInGE.contrib.analysis import analyzeGPS
    # PObjU= PlotNYC(Exp1, kind='unc')
    # PObjU.plot_gps()
    # PObjU.plot_basic_box(pts='Ashe Woodside'.split())

    PObjR = PlotNYC(Exp1, kind='rate')
    # PObjR.plot_insets(pts='Ashe Woodside'.split())
    PObjR.plot_transect()

    # GPSObj = analyzeGPS.GPS_Rate_Analysis(Exp1)
    # GPSObj.plot_gps_vs_insar()
    # GPSObj.plot_trends_new()
    # analyzeGPS.plot_trends_new(Exp1)

    # plot_hist_ref_sta()
    # plot_hist_corr()

    # PObjU.plot_polys()


def make_numbers(Exp1):
    PObj  = PlotNYC(Exp1, kind='rate')
    PObj1  = PlotNYC(Exp1, kind='unc')

    attrs = readfile.read_attribute(PObj.path_vup_geo)
    log.critical(f'Experiment: {PObj.mp_exp} (neofs={PObj.neofs})')
    log.critical('%s to %s', attrs['START_DATE'], attrs['END_DATE'])

    l, h = np.nanquantile(PObj.da, 0.003), np.nanquantile(PObj.da, 0.997)
    m    = np.nanmedian(PObj.da)

    l1, h1 = np.nanquantile(PObj1.da, 0.003), np.nanquantile(PObj1.da, 0.997)
    m1   = np.nanmedian(PObj1.da)

    ## overall InSAR mean
    log.info (f'InSAR Vup: {m:.2f} [{l:.2f}, {h:.2f}]')
    log.info (f'InSAR Unc: {m1:.2f} [{l1:.2f}, {h1:.2f}]')

    ## overall GPS mean
    print ('Stations for mean:', ', '.join(PObj.df_gps.sta.tolist()))
    m2     = PObj.df_gps.u_vel.median()
    l2, h2 = PObj.df_gps.u_vel.quantile(0.003), PObj.df_gps.u_vel.quantile(0.997)

    m3     = PObj.df_gps.u_sig.median()
    l3, h3 = PObj.df_gps.u_sig.quantile(0.003), PObj.df_gps.u_sig.quantile(0.997)

    log.info(f'MIDAS rate: {m2:.2f}, [{l2:.2f}, {h2:2f}] mm/yr')
    log.info(f'MIDAS unc: {m3:.2f} +/- [{l3:.2f}, {h3:.2f}] mm/yr')

    ## percent significant
    PObj1  = PlotNYC(Exp1, kind='unc')

    # get number of nonnan pix, aka the actual rates we are considering
    n_pix  = PObj.da.where(np.isnan(PObj.da), 1).sum()

    # get a count of where pixels are >= sigma
    da_sig  = xr.where(np.abs(PObj.da)>=PObj1.da, 1, 0)
    per_sig = 100*(da_sig.sum() / n_pix)
    log.info (f'Percent significant at 1sig: {per_sig:.1f}%')


    ## RMSE
    PObj.plot_rmse(show=False)

    ## cities
    for city in 'Manhattan Brooklyn Queens Hudson Essex'.split() + \
        ['Staten Island', 'Bergen County']:
        loc = 4 if city == 'Essex' else 0

        PObj.plot_cities(city, loc, show=False)
        # PObj1.plot_cities(city, loc, show=False)
    return


def diagnostics(Exp1):
    from VLM.bzFRInGE.contrib.analysis import analyzeGPS
    PObjR = PlotNYC(Exp1, kind='rate')
    PObjU = PlotNYC(Exp1, kind='unc')
    PObjR.plot_basic()

    n_pix  = PObjR.da.where(np.isnan(PObjR.da), 1).sum()

    # get a count of where pixels are >= sigma
    da_sig  = xr.where(np.abs(PObjR.da)>=PObjU.da, 1, 0)
    per_sig = 100*(da_sig.sum() / n_pix)
    log.critical (f'Percent significant at 1sig: {per_sig:.1f}%')

    GPSObj = analyzeGPS.GPS_Rate_Analysis(Exp1)
    GPSObj.plot_gps_vs_insar()
    # GPSObj.plot_trends_new()
    return

if __name__ == '__main__':
    # Exp       = ExpBase(NYC_SRc, 'ERA5_SET_PM_Fast', 'NYBK', neofs=20)
    # Exp       = ExpBase(NYC_SRc, 'ERA5_SET_PM_NJHT_Fast', 'NYBK', neofs=20)
    Exp0      = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20)
    # make_main(Exp0)
    # make_supp(Exp0)
    # make_numbers(Exp0)



    path_figs = op.join(op.expanduser('~'), 'Desktop', 'VLM', 'NYC_2023')
    # bbPlot.savefigs(path_figs, False, True)

    plt.show()