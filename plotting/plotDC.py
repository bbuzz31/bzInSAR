from VLM import *
from VLM.bzFRInGE.plotting.plotExp import *
from VLM.bzFRInGE.plotting.plot_TS1 import PlotTS1, DCT_PTS
from mintpy.utils import utils as ut
import geopandas as gpd
from shapely.geometry import Point
from matplotlib.gridspec import GridSpec#, GridSpecFromSubplotSpec

# DCT_CITIES['DC'] = ["Washington DC", "Prince George's", 'Montgomery, MD',
#             'Alexandria, VA', 'Arlington, VA', 'Fairfax County', 'Falls Church']
DCT_CITIES['DC'] = ["Washington DC", 'Alexandria, VA']#, 'Anacostia', 'Fairlawn, Washington D.C.']

class PlotDC(PlotExp):
    def __init__(self, exp, kind='rate', show_gps=True, continuous=True, sten=''):
        super().__init__(exp, kind, show_gps, continuous, sten)
        self._set_plot_parms()
        self.get_rate_unc_h5()
        # self.update_unc()


    def _set_plot_parms(self):
        super()._set_plot_parms()
        self.gs = 15
        if self.kind == 'Rate':
            self.show_gps_names=True

        elif self.kind == 'Unc' and self.sten:
            self.pparms['norm'] = Normalize(1.0, 2.0)
        return


    def plot_basic_DC(self, show_cities=True, show_names=False, show_tg=True, show_fl=True, axes=None):
        self.gps_s = 100 # increase GPS sizes
        self.gps_fs = 14
        self.fs1 = 16
        # temp to crop
        # gdf = gpd.read_file(Path.home() / 'Desktop' / 'Test_DC_Crop.GeoJSOn')
        gdf = make_crop()
        self.da.rio.write_crs(4326, inplace=True)
        self.da = self.da.rio.clip(gdf.geometry, all_touched=True, drop=True)

        fig, axes  = self.plot_basic(axes)
        if show_tg:
            lalo = (38.873333,  -77.021667)
            txt_parms = dict(path_effects=[pe.withStroke(linewidth=3.0, foreground='k')])
            pt_parms = dict(edgecolors='yellow', facecolor='none',
                            transform=ccrs.PlateCarree(), s=200, marker='v')
            axes.scatter(lalo[1], lalo[0], **pt_parms, **txt_parms)

        if show_cities:
            lst_polys = [city_shape(city, 0)[0] for city in DCT_CITIES[self.reg]]
            gdf_cities = pd.concat(lst_polys).set_index('name')
            # plot dc in a different color
            # gser_dc = gdf_cities.loc[['Washington']]
            # gdf_cities.drop('Washington', inplace=True)
            sty  = dict(zorder=50, alpha=1.0, facecolor='none', transform=self.proj,
                        edgecolor='w', linestyle='--', linewidth=2.0)
            gdf_cities.plot(ax=axes, **sty)
            # gser_dc.plot(ax=axes, **{**sty, 'edgecolor':'deeppink'})

            # not implemented
            if show_names:
                for city in DCT_CITIES:
                    # y = 30 if name == 'NYBK' else -20
                    sty_name = dict(facecolor='lightgray', alpha=0.5, boxstyle='round')
                    y, halign = 0, 'left'
                    rot = 61 if city == 'Manhattan' else 0
                    rot = 65 if city == 'Hudson' else rot
                    x, y, txt = *DCT_NAMES[city], city.replace('County', '')
                    geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                    text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=y)
                    axes.text(x, y, txt, transform=text_transform, zorder=20, rotation=rot,
                            verticalalignment='center', horizontalalignment=halign,
                            bbox=sty_name, size=11.5, color='black')
        if show_fl:
            txt_parms = dict(path_effects=[pe.withStroke(linewidth=1.5, foreground='k')])
            gdf_fl = gpd.read_file(self.path_shapes / 'Fall_Line.GeoJSON')
            gdf_fl.plot(ax=axes, transform=self.proj, color='sandybrown',
                        zorder=30, linewidth=2, alpha=1)
            axes.text(-77.19, 38.79, 'Fall Line', **txt_parms,
                    fontsize=18, color='sandybrown', zorder=30,
                    transform=self.proj, rotation=55)

        return fig, axes


    def plot_transect_map(self):
        self.gps_s = 150 # increase GPS sizes
        self.gps_fs = 18
        self.show_gps_names = False

        transects = 'West_Subsidence East_Uplift'.split()
        fig, axes = super().plot_transect_map(transects, fontsize=24)

        ## move the colorbar to the right side
        # cbar = axes.im.get_cbar()
        im  = axes.get_children()[0]
        im.colorbar.remove()
        dct_ps = dict(ylabel='', size='5%', pad=0.15, pack_start=True)
        divider = make_axes_locatable(im.axes)
        cbar_ax = divider.new_horizontal(**dct_ps, axes_class=plt.Axes)
        cbar_ax = fig.add_axes(cbar_ax)

        cbar  = plt.colorbar(im, cax=cbar_ax)

        cbar.set_label('Vertical Rate (mm/yr)', rotation=90,
                       labelpad=0, fontsize=18)
        cbar.ax.yaxis.set_ticks_position('left')
        cbar.ax.yaxis.set_label_position('left')
        cbar.ax.tick_params(labelsize=20)

        ## close up
        axes.set_extent([-77.115, -76.965, 38.75, 38.96])

        self._add_ART(axes, fontsize=24, linewidth=5)

        ## add the tide gauge
        pt_parms = dict(edgecolors='white', facecolor='none',
                        transform=ccrs.PlateCarree(), s=70, marker='v')
        axes.scatter(-77.021667, 38.873333, **pt_parms)

        ## scale bar
        axes = bbPlot._scale_bar(axes, 5, (0.77, 0.95), fs=16, outline=2.0)

        ## gridlines
        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, left=False, bottom=True, size=17)

        return


    def get_well_locs(self):
        gdf_wells = gpd.read_file(self.path_shapes / 'DC_GW_rates.GeoJSON')

        #get the closest reference GPS station; for stitching not really needed
        # ref_stas = 'NASA UHDT CFJV SESG TXTG NASA CSTE DMFB TXB6 TXLQ'.split()
        # gdf_refs  = bbGIS.df2gdf(self.df_gps[self.df_gps.sta.isin(ref_stas)].set_index('sta'), 'all')
        # gdf_wells = gdf_wells.to_crs(4087).sjoin_nearest(gdf_refs.to_crs(4087), distance_col='dist').to_crs(4326)
        # gdf_wells['dist'] /= 1000
        return gdf_wells


    def add_bbox(self, axes, SNWE):
        """ Add a bounding box on top of the main map and a scale bar"""
        sty  = dict(zorder=50, alpha=0.9, facecolor='none', edgecolor='white',
                    linestyle='--', linewidth=1.5)
        bbox = bbGIS.bbox2poly(SNWE=SNWE)
        axes.add_geometries([bbox], crs=self.proj, **sty)

        ## save the bbox for QGIS
        # gdf = gpd.GeoDataFrame(geometry=[bbox], crs='epsg:4326')
        # gdf.to_file(f'{self.path_wd}/Houston_box2.GeoJSON')

        return axes


    def plot_inset(self, loc='Navy_Yard', npix=0, double_diff=True):
        """ Plot one point, row wise """
        log.warning('Not double differencing...') if not double_diff else ''
        self.show_gps_names = False
        fig = plt.figure(figsize=(18, 10), constrained_layout=True)
        gs  = GridSpec(2, 3, width_ratios=[1.3, 1, 1], figure=fig)
        axes = []
        axes.append(plt.subplot(gs[0, 0], projection=self.basemap.crs))
        axes.append(plt.subplot(gs[0, 1:]))

        ## Plot the Map
        lalo = DCT_PTS[loc][0] #  use plot_ts to make a box around the point
        S, N, W, E = buffer_point(*lalo, 2000)

        self.zoom  = 15

        ax  = self.plot_basic(axes=axes[0])[1]
        im  = ax.get_children()[0]
        im.colorbar.remove()
        im.set_alpha(0.75)

        ## add actual point
        pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                        color='w', transform=ccrs.PlateCarree(), s=50)
        ax.scatter(lalo[1], lalo[0], **pt_parms)
        ## not using this kind of stable anymore
        # if show_stable:
        #     lalo1 = DCT_STABLE[self.mp_exp0][f'{loc}_Stable'][0]
        #     ax.scatter(lalo1[1], lalo1[0], marker='s', **pt_parms)

        # else:
        ax.set_extent([W, E, S, N], crs=self.proj)

        if loc == 'Navy_Yard':
            bbPlot._scale_bar(ax, 1, (0.80, 0.015))
        else:
            bbPlot._scale_bar(ax, 1, (0.20, 0.015))

        ## make colorbar
        divider = make_axes_locatable(axes[0])
        cbar_ax = divider.new_horizontal(size="5%", pad=0.60, axes_class=plt.Axes, pack_start=True)
        fig.add_axes(cbar_ax)

        gl = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 190, 0.01),
                    ylocs=np.arange(-90, 100, 0.01))

        # gparms = dict(left=False, right=True)# if i == 1 else {}
        bbPlot.fmt_gridlines(gl, left=False, right=True, bottom=True, size=self.gs-1)#, **gparms)

        cbar = fig.colorbar(im, cbar_ax)
        cbar.set_label('Vertical Rate (mm/yr)', fontsize=13)
        # cbar_ax.set_xlabel('mm/yr', labelpad=10)
        cbar.ax.yaxis.set_label_position('left')
        cbar_ax.yaxis.set_ticks_position('left')
        # leg = axes[1].legend(loc='lower left', frameon=True, shadow=True, ncol=1, prop={'size':'8'})

        ## Plot the Timeseries
        TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs, npix)
        if double_diff:
            axeTS = TSObj.plot_ts_vup_dd(loc, axes=axes[1])[1]
        else:
            axeTS = TSObj.plot_ts_vup(loc, axes=axes[1])[1]
        axeTS.set_title('')
        # axeTS.set_ylim([-40, 20]) if loc == 'LaGuardia' else ''
        # axeTS.set_ylim([-20, 20])# if loc == 'LaGuardia' else ''
        ## this is done already
        # ymin, ymax = int(np.floor(ts_dd.min()))-5, int(np.ceil(ts_dd.max())+5)
        # ymin, ymax = (-10, 10) if ymin > -10 and ymax < 10 else (ymin, ymax)
        # axes.set_ylim([ymin, ymax])

        axeTS.set_ylabel('Vertical Displacment (mm)', fontsize=13, rotation=270, labelpad=20)
        axeTS.yaxis.set_label_position('right')
        axeTS.yaxis.set_ticks_position('right')
        fig.set_label(f'Inset_{loc}')
        return


    def plot_insets(self, pts='Jefferson_Memorial Anacostia'.split(), double_diff=True, npix=1):
        """ Plot two points, row wise """
        log.warning('Not double differencing...') if not double_diff else ''
        self.show_gps = False
        fig = plt.figure(figsize=(14, 9), constrained_layout=False)
        gs  = GridSpec(2, 3, width_ratios=[1, 1, 1], hspace=0.33, wspace=0.0, figure=fig)
        axeMaps = []
        axeTSs  = []
        axeMaps.append(plt.subplot(gs[0, 0], projection=self.basemap.crs))
        axeMaps.append(plt.subplot(gs[1, 0], projection=self.basemap.crs))
        axeTSs.append(plt.subplot(gs[0, 1:]))
        axeTSs.append(plt.subplot(gs[1, 1:]))
        for i, pt in enumerate(pts):
            lalo = DCT_PTS[pt][0] #  use plot_ts to make a box around the point
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

            # gl = ax.gridlines(draw_labels=True)#, xlocs=np.arange(-180, 190, 0.01),
                        # ylocs=np.arange(-90, 100, 0.01))

            # bbPlot.fmt_gridlines(gl, left=True, right=False, bottom=True, size=GS-1)

            TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs, npix)
            show_fit = True
            # if pt == 'Anacostia':
            #     axeTS = TSObj.plot_ts_vup_piece2('20191231', pt, double_diff, npix, axes=axeTSs[i])[1]
            # else:
            #     show_fit = False if i == 1 else True
            if double_diff:
                axeTS = TSObj.plot_ts_vup_dd(pts[i], axes=axeTSs[i], show_fit=show_fit)[1]
            else:
                axeTS = TSObj.plot_ts_vup(pts[i], axes=axeTSs[1])[1]

            # pos = axeTS.get_position().bounds
            # pos[0] = pos[0]-0.1
            # axeTS.set_position(mpl.transforms.Bbox.from_bounds(*pos))

            axeTS.set_title('')
            ylims = [-11, 6] if i == 0 else [-1, 21]
            axeTS.set_ylim(ylims)
            axeTS.tick_params(labelbottom=False) if i == 0 else ''
            axeTS.tick_params(axis='both', which='major', labelsize=13)

            if i == 1:
                axeTS.set_ylabel(axeTS.get_ylabel(), labelpad=12)

            axeTS.set_ylabel('Vertical Displacment (mm)', fontsize=13,
                             rotation=270, labelpad=20)
            axeTS.yaxis.set_label_position('right')
            axeTS.yaxis.set_ticks_position('right')

        # make a single colorbar across
        cbar_ax = fig.add_axes([0.165, 0.475, 0.160, 0.0160]) # left, bot, wid, hgt
        cbar_ax.set_title('V$_{Up}$  (mm yr$^{-1}$)', fontsize=11)#, labelpad=1v#, rotation=270, )
        # cbar_ax.yaxis.set_label_position('right')
        cbar    = plt.colorbar(im, cax=cbar_ax, orientation='horizontal', spacing='uniform')
        # gs.update(left=0.07, right=0.23)#, bottom=0.1, top=0.9)
        # gs.update(hspace=-0.5)#, wspace=10.00)
        fig.set_label(f'Insets_{"_".join(pts)}')
        fig.constrained_layout=True
        return


    def plot_inset_grad(self, md_locs=['MD_Up', 'MD_Down'], npix=0, double_diff=True):
        """ Plot one point, row wise """
        log.warning('Not double differencing...') if not double_diff else ''
        self.show_gps = False
        self.show_gps_names = False

        fig = plt.figure(figsize=(18, 10), constrained_layout=True)
        gs  = GridSpec(2, 3, width_ratios=[1.3, 1, 1], figure=fig)
        axes = []
        axes.append(plt.subplot(gs[0, 0], projection=self.basemap.crs))
        axes.append(plt.subplot(gs[0, 1:]))

        ## Plot the Map
        # lalo = DCT_PTS[loc][0] #  use plot_ts to make a box around the point
        # S, N, W, E = buffer_point(*lalo, 2000)

        self.zoom  = 15

        ax  = self.plot_basic(axes=axes[0])[1]
        im  = ax.get_children()[0]
        im.colorbar.remove()
        # im.set_alpha(0.75)

        ## add actual point
        WESN_crop = [-77.075, -76.85, 38.75, 38.920] # [-77.075, -76.9, 38.745, 38.865]
        ax.set_extent(WESN_crop, crs=self.proj)

        ## make colorbar
        divider = make_axes_locatable(axes[0])
        cbar_ax = divider.new_horizontal(size="5%", pad=0.60, axes_class=plt.Axes, pack_start=True)
        fig.add_axes(cbar_ax)

        gl = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 190, 0.01),
                    ylocs=np.arange(-90, 100, 0.01))

        # gparms = dict(left=False, right=True)# if i == 1 else {}
        bbPlot.fmt_gridlines(gl, left=False, right=False, bottom=False, size=self.gs-2)#, **gparms)

        cbar = fig.colorbar(im, cbar_ax)
        cbar.set_label('Vertical Rate (mm/yr)', fontsize=13)
        cbar.ax.yaxis.set_label_position('left')
        cbar_ax.yaxis.set_ticks_position('left')

        ## add the wells to the map
        msty = dict(marker='o', s=200, linewidth=1.5, transform=self.proj,
                    alpha=1, zorder=30, edgecolors='black', c='white')
        gdf_wells = gpd.read_file(self.path_vlm / 'GW' / 'DC_GW_rates.GeoJSON')
        gdf_wells = bbGIS.in_dfll(gdf_wells, WESN=WESN_crop)
        s = axes[0].scatter(gdf_wells.geometry.x, gdf_wells.geometry.y, **msty)

        txt_parms = dict(path_effects=[pe.withStroke(linewidth=1.75, foreground='k')],
                            fontsize=18, color='white')
        names = gdf_wells['name'].astype('str').values
        for i, name in enumerate(names):
            axes[0].text(gdf_wells.geometry.x.iloc[i], gdf_wells.geometry.y.iloc[i],
                    name.strip(), transform=self.proj, zorder=80,
                    verticalalignment='center', **txt_parms)

        self._add_hydrogeo(axes[0])

        pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                        color='yellow', transform=self.proj, s=175)

        ## lbl and scalebar
        axes[0].text(0, 1.01, 'A', transform=axes[0].transAxes,
                    verticalalignment='top', fontsize=30, color='w',
                    bbox=dict(boxstyle='square,pad=0.05', facecolor='black'))
        bbPlot._scale_bar(axes[0], 5, (0.85, 0.015), fs=26, outline=2.0)

        TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs, npix)
        # TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, neofs=0, npix=npix)
        lst_ts, lst_coeffs, lst_uncs = [], [], []

        for loc in md_locs:
            lalo = DCT_PTS[loc][0]
            marker = '^' if 'Up' in loc else 'v'
            ax.scatter(lalo[1], lalo[0], marker=marker, **pt_parms)

            if double_diff:
                ts, coeffs, unc = TSObj.calc_ts_vup_dd_atm(loc, 500, 0.95, show_ts=False)
            else:
                ts, coeffs, unc = TSObj.calc_ts_vup(loc)
            lst_ts.append(ts)
            lst_coeffs.append(coeffs)
            lst_uncs.append(unc)

        ## compute the difference of the two and a new trend / unc
        ts_diff = lst_ts[0] - lst_ts[1]
        axeTS = axes[-1]
        # TSObj.plot_ts_vup_piece2('20201231', ts_diff, double_diff, npix=npix, axes=axeTS)
        TSObj.plot_ts_vup_dd(ts_diff, axes=axeTS, show_fit=False)

        ## i removed all this full linear trend and plotting in favor of pieces
        # coeffs_diff, rss_diff = np.polyfit(TSObj.decyr, ts_diff, 1, full=True)[:2]
        # rmse_diff = np.sqrt(rss_diff/(len(ts_diff)-2)).item()

        # rate_diff = coeffs_diff[0]
        # fit_diff = np.polyval(coeffs_diff, TSObj.decyr)
        # unc_diff = np.sqrt(lst_uncs[0]**2 + lst_uncs[1]**2)

        ## plot the timeseries
        ## use plot_gradient_ts_atm to see the individual ts
        # col = 'r' if rate_diff > 0 else 'b'
        # col = 'k' if len(coeffs_diff) > 2 else col # nonlinear fits
        # axeTS.scatter(TSObj.dt, ts_diff, c='k', s=15)
        # if show_fit: # show a line; remove to see nonlinearity better
        #     axeTS.plot(TSObj.dt, fit_diff, color=col, linestyle='--',
        #             label=f'{rate_diff:.1f}$\\pm${unc_diff:.1f} (mm/yr)')
        #     axeTS.fill_between(TSObj.dt, fit_diff-unc_diff, fit_diff+unc_diff, color=col, alpha=0.3)
        #     axeTS.legend(prop ={"size": 13}, loc='best')
        axeTS.tick_params(axis='both', labelsize=13)
        # axeTS.set_ylim([-1, 15])
        axeTS.set_ylim([-10, 10])
        # axeTS.set_ylim([-20, 20])
        axeTS.grid(color='k', linestyle = '--', alpha=0.1)
        axeTS.set_ylabel('$\\Delta$ Vertical Displacment (mm)', fontsize=13, rotation=270, labelpad=20)
        axeTS.yaxis.set_label_position('right')
        axeTS.yaxis.set_ticks_position('right')
        axeTS.text(-0.01, 1.01, 'B', transform=axeTS.transAxes,
                    verticalalignment='top', fontsize=30, color='w',
                    bbox=dict(boxstyle='square,pad=0.05', facecolor='black'))
        ext = loc[-1] if loc[-1].isdigit() else ''
        fig.set_label(f'Inset_grad{ext}')
        return


    def plot_gradient_map(self, buff=500, show_box=True, lbl='A', axes=None):
        """ Plot maps close up along the gradient with the GW contour

        Optionally show the box where you'll average and compute ts_diff
        Optionally show wells
        """
        self.show_gps = False
        # WESN_crop = [-77.075, -76.85, 38.745, 38.880] # [-77.075, -76.9, 38.745, 38.865]
        WESN_crop = [-77.075, -76.85, 38.75, 38.920] # [-77.075, -76.9, 38.745, 38.865]

        fig, axes = self.plot_basic_DC(False, False, False, False, axes=axes)
        axes.set_extent(WESN_crop)

        im  = axes.get_children()[1]
        self._add_hydrogeo(axes)

        if show_box:
            md_locs = ['MD_Up', 'MD_Down']
            pts = [Point(DCT_PTS[loc][0][::-1]) for loc in md_locs]
            gser = gpd.GeoSeries(pts, crs=4326)
            gser_buff = gser.to_crs(4087).buffer(buff, cap_style='round').to_crs(4326)
            sty  = dict(zorder=50, alpha=1.0, facecolor='none', edgecolor='k',
                        linestyle=':', linewidth=3.5)
            sty2 = {**sty, 'edgecolor':'yellow', 'linestyle':'-', 'linewidth': 4.0}

            ##  the point in the middle
            pt_parms1 = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                            color='yellow', transform=self.proj, s=150, marker='X')
            # pt_parms2 = {**pt_parms1, 'marker': 'X'}
            pt_parms2 = pt_parms1

            if buff:
                for poly in gser_buff:
                    w, s, e, n = poly.bounds
                    # white solid
                    gser_buff.plot(ax=axes, transform=self.proj, **sty2)
                    # black dashed
                    gser_buff.plot(ax=axes, transform=self.proj, **sty)

            else:
                # add middle point
                axes.scatter(gser.geometry.x[0], gser.geometry.y[0],  **pt_parms1)
                axes.scatter(gser.geometry.x[1], gser.geometry.y[1],  **pt_parms2)

        if lbl == 'A':
            msty = dict(marker='o', s=200, linewidth=1.5, transform=self.proj,
                        alpha=1, zorder=30, edgecolors='black', c='white')
            gdf_wells = gpd.read_file(self.path_vlm / 'GW' / 'DC_GW_rates.GeoJSON')
            gdf_wells = bbGIS.in_dfll(gdf_wells, WESN=WESN_crop)
            s = axes.scatter(gdf_wells.geometry.x, gdf_wells.geometry.y, **msty)

            sty_name = dict(facecolor='lightgray', alpha=0.5, boxstyle='round')
            txt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                                fontsize=30, color='white')

            names = gdf_wells['name'].astype('str').values
            for i, name in enumerate(names):
                axes.text(gdf_wells.geometry.x.iloc[i], gdf_wells.geometry.y.iloc[i],
                        name.strip(), transform=self.proj, zorder=80,
                        verticalalignment='center', **txt_parms)

        # if lbl == 'B':
            # self._add_ART(axes, fontsize=24, linewidth=10)


        ## add the label
        pos = 0 #if lbl == 'A' else 1
        axes.text(pos, 1.01, lbl, transform=axes.transAxes,
                    verticalalignment='top', fontsize=30, color='w',
                    bbox=dict(boxstyle='square,pad=0.05', facecolor='black'))

        bbPlot._scale_bar(axes, 5, (0.85, 0.015), fs=26, outline=2.0)
        ext = md_locs[0][-1]
        fig.set_label(f'{self.mp_exp}_gradmap_{lbl}_unc_MD{ext}')
        return fig, axes


    def plot_gradient_ts(self, buff=500, axes=None, show_each=True):
        from BZ.bbTS import bootstrap_1D

        locs = 'MD_Up', 'MD_Down'
        pts = [Point(DCT_PTS[loc][0][::-1]) for loc in locs]
        gser = gpd.GeoSeries(pts, crs=4326)
        gser_buff = gser.to_crs(4087).buffer(buff, cap_style='round').to_crs(4326)
        sty = dict(zorder=50, alpha=0.9, facecolor='none', edgecolor='white',
                    linestyle='--', linewidth=1.5)

        TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs)

        lst_ts, lst_coeffs, lst_uncs = [], [], []
        for i, _ in enumerate(gser_buff):
            gser_poly = gser_buff.iloc[[i]]
            ## add triangle
            ts, coeffs, unc = TSObj.calc_ts_vup_poly(gser_poly)
            lst_ts.append(ts)
            lst_coeffs.append(coeffs)
            lst_uncs.append(unc)

        ## compute the difference of the two and a new trend / unc
        ts_diff = lst_ts[0] - lst_ts[1]
        coeffs_diff = np.polyfit(TSObj.decyr, ts_diff, 1)
        fit_diff = np.polyval(coeffs_diff, TSObj.decyr)
        unc_diff = np.sqrt(lst_uncs[0]**2 + lst_uncs[1]**2)

        ## plot the individual timeseries
        if show_each:
            f, a = plt.subplots(figsize=(10, 6), sharex=True, nrows=3)
            for i, loc in enumerate(locs):
                fit = np.polyval(lst_coeffs[i], TSObj.decyr-TSObj.decyr[0])
                rate, unc = lst_coeffs[i][0], lst_uncs[i]
                a[i].scatter(TSObj.dt, lst_ts[i], color='k', s=15)
                l, = a[i].plot(TSObj.dt, fit, color='darkred', linestyle='--',
                                label=f'{rate:.1f}$\\pm${unc:.1f}')
                a[i].fill_between(TSObj.dt, fit-unc, fit+unc, color=l.get_color(), alpha=0.3)
                a[i].set_title(f'{locs[i]} (Buff={buff} m)', fontsize=14)

            a[2].scatter(TSObj.dt, ts_diff, color='k', s=15)
            l, = a[2].plot(TSObj.dt, fit_diff, color='darkred', linestyle='--',
                            label=f'{coeffs_diff[0]:.1f}$\\pm${unc_diff:.1f}')
            a[2].fill_between(TSObj.dt, fit_diff-unc_diff, fit_diff+unc_diff, color=l.get_color(), alpha=0.3)
            a[2].set_title(f'{locs[0]} + {locs[1]}', fontsize=14)
            for ax in a:
                ax.set_ylabel('V$_{Up}$ (mm)', fontsize=13)
                ax.grid(color='k', alpha=0.1, linestyle='--')
                # ax.set_ylim([-10, 10])
                ax.legend()
            f.subplots_adjust(hspace=0.5)
            f.set_label(f'{locs[0]}-{locs[1]}_ts')

        ## now the moving
        lst_rates, lst_uncs, lst_sts = [], [], []
        en = TSObj.dt[-1]
        for i in range(len(TSObj.decyr)):
            st = TSObj.dt[i]
            # require a 3 year timeseries
            if (en-st).days < (365.25 * 2):
                break
            ts = ts_diff[i:].data
            decyr = TSObj.decyr[i:]
            coeffs = np.polyfit(decyr, ts, 1)

            rate, unc, bias = bootstrap_1D(ts, decyr, n_boots=2500)

            lst_rates.append(rate)
            lst_uncs.append(unc)
            lst_sts.append(TSObj.dt[i])

        # filter out uncertainties
        df_mt = pd.DataFrame({'rate': lst_rates, 'unc': lst_uncs}, index=lst_sts)
        df_mt['1sig'] = df_mt['rate'].abs() >= df_mt['unc']
        df_mt['2sig'] = df_mt['rate'].abs() >= 2*df_mt['unc']

        print (f"delRate falls below 0.8 at {df_mt[(df_mt['rate']+df_mt['unc'])<0.8].index[0]}")
        # this justifies next statement
        # df_mt[df_mt.index>pd.to_datetime('20191201')].head(10)
        print (f"delRate mostly insignificantly different than 0 after 202002")

        df_plot = df_mt # df_mt[df_mt['1sig']] # no point unless showing scatters

        log.info (f'Last start date for trend computation: {st}')
        log.info (f'Dropped {df_mt.shape[0] - df_plot.shape[0]} insignificant dates')

        if axes is None:
            fig, axes = plt.subplots(figsize=(10, 3))
        else:
            fig = axes.get_figure()

        axes.scatter(lst_sts, lst_rates, color='k', s=15)
        l, = axes.plot(df_plot.index, df_plot['rate'], color='k', alpha=1.0)
        axes.fill_between(df_plot.index, df_plot['rate']-df_plot['unc'], df_plot['rate']+df_plot['unc'],
                          color=l.get_color(), alpha=0.3)
        axes.set_ylabel('$\\Delta$ Vup (mm/yr)', fontsize=14)
        axes.set_xlabel('Start Date', fontsize=14)
        # axes.set_title(f'{locs[0]}-{locs[1]} Moving Trends', fontsize=16)
        axes.grid(color='k', linestyle='--', alpha=0.1)

        lbl = 'C'
        pos = 0
        axes.text(pos, 1.01, lbl, transform=axes.transAxes,
                    verticalalignment='top', fontsize=20, color='w',
                    bbox=dict(boxstyle='square,pad=0.05', facecolor='black'))
        fig.set_label(f'{locs[0]}_{locs[1]}_moving')
        return fig, axes


    def plot_gradient_ts_atm(self, axes=None, show_each=True):
        from BZ.bbTS import bootstrap_1D
        locs = 'MD_Up', 'MD_Down'
        # pts = [Point(DCT_PTS[loc][0][::-1]) for loc in locs]
        TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs)
        lst_ts, lst_coeffs, lst_uncs = [], [], []
        for loc in locs:
            ts, coeffs, unc = TSObj.calc_ts_vup_dd_atm(loc, 500, 0.95, show_ts=False)
            lst_ts.append(ts)
            lst_coeffs.append(coeffs)
            lst_uncs.append(unc)

        ## compute the difference of the two and a new trend / unc
        ts_diff = lst_ts[0] - lst_ts[1]
        coeffs_diff = np.polyfit(TSObj.decyr, ts_diff, 1)
        fit_diff = np.polyval(coeffs_diff, TSObj.decyr)
        unc_diff = np.sqrt(lst_uncs[0]**2 + lst_uncs[1]**2)

        ## plot the individual timeseries
        if show_each:
            f, a = plt.subplots(figsize=(10, 6), sharex=True, nrows=3)
            for i, loc in enumerate(locs):
                fit = np.polyval(lst_coeffs[i], TSObj.decyr)
                rate, unc = lst_coeffs[i][0], lst_uncs[i]
                a[i].scatter(TSObj.dt, lst_ts[i], color='k', s=15)
                l, = a[i].plot(TSObj.dt, fit, color='darkred', linestyle='--',
                                label=f'{rate:.1f}$\\pm${unc:.1f}')
                a[i].fill_between(TSObj.dt, fit-unc, fit+unc, color=l.get_color(), alpha=0.3)
                a[i].set_title(f'{locs[i]}$_{{-Atm}}$', fontsize=14)

            a[2].scatter(TSObj.dt, ts_diff, color='k', s=15)
            l, = a[2].plot(TSObj.dt, fit_diff, color='darkred', linestyle='--',
                            label=f'{coeffs_diff[0]:.1f}$\\pm${unc_diff:.1f}')
            a[2].fill_between(TSObj.dt, fit_diff-unc_diff, fit_diff+unc_diff, color=l.get_color(), alpha=0.3)
            a[2].set_title(f'{locs[0]} - {locs[1]}', fontsize=14)
            for ax in a:
                ax.set_ylabel('V$_{Up}$ (mm)', fontsize=13)
                ax.grid(color='k', alpha=0.1, linestyle='--')
                # ax.set_ylim([-10, 10])
                ax.legend()
            f.subplots_adjust(hspace=0.5)
            f.set_label(f'{locs[0]}-{locs[1]}_ts')

        ## now the moving
        lst_rates, lst_uncs, lst_sts = [], [], []
        en = TSObj.dt[-1]
        for i in range(len(TSObj.decyr)):
            st = TSObj.dt[i]
            # require a 3 year timeseries
            if (en-st).days < (365.25 * 2):
                break
            ts = ts_diff[i:]
            decyr = TSObj.decyr[i:]
            coeffs = np.polyfit(decyr, ts, 1)

            rate, unc, bias = bootstrap_1D(ts, decyr, n_boots=2500)

            lst_rates.append(rate)
            lst_uncs.append(unc)
            lst_sts.append(TSObj.dt[i])

        # filter out uncertainties
        df_mt = pd.DataFrame({'rate': lst_rates, 'unc': lst_uncs}, index=lst_sts)
        df_mt['1sig'] = df_mt['rate'].abs() >= df_mt['unc']
        df_mt['2sig'] = df_mt['rate'].abs() >= 2*df_mt['unc']

        print (f"delRate falls below 0.8 in {df_mt[(df_mt['rate']+df_mt['unc'])<0.8].index[0]}")
        # this justifies next statement
        # df_mt[df_mt.index>pd.to_datetime('20191201')].head(10)
        print (f"delRate mostly insignificantly differnt than 0 after 202002")

        df_plot = df_mt # df_mt[df_mt['1sig']] # no point unless showing scatters

        log.info (f'Last start date for trend computation: {st}')
        log.info (f'Dropped {df_mt.shape[0] - df_plot.shape[0]} insignificant dates')

        if axes is None:
            fig, axes = plt.subplots(figsize=(10, 3))
        else:
            fig = axes.get_figure()

        axes.scatter(lst_sts, lst_rates, color='k', s=15)
        l, = axes.plot(df_plot.index, df_plot['rate'], color='k', alpha=1.0)
        axes.fill_between(df_plot.index, df_plot['rate']-df_plot['unc'], df_plot['rate']+df_plot['unc'],
                          color=l.get_color(), alpha=0.3)
        axes.set_ylabel('$\\Delta$ Vup (mm/yr)', fontsize=14)
        axes.set_xlabel('Start Date', fontsize=14)
        # axes.set_title(f'{locs[0]}-{locs[1]} Moving Trends', fontsize=16)
        axes.grid(color='k', linestyle='--', alpha=0.1)

        lbl = 'C'
        pos = 0
        axes.text(pos, 1.01, lbl, transform=axes.transAxes,
                    verticalalignment='top', fontsize=20, color='w',
                    bbox=dict(boxstyle='square,pad=0.05', facecolor='black'))
        fig.set_label(f'{locs[0]}_{locs[1]}_moving')
        return fig, axes


    def plot_tsgrad(self, npix=0):
        """" Plot the DIFFERENCE between the two double differenced timeseries across grad

        Difference so that it goes towards 0 as they become the same
        """
        locs = 'MD_Up', 'MD_Down'
        TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs, npix)
        ts_dd1, coeffs_dd1, sem21 = TSObj.calc_ts_vup_dd(locs[0])
        rate1, bias1 = coeffs_dd1
        fit1 = np.polyval(coeffs_dd1, TSObj.decyr)

        ts_dd2, coeffs_dd2, sem22 = TSObj.calc_ts_vup_dd(locs[1])
        rate2, bias2 = coeffs_dd2
        fit2 = np.polyval(coeffs_dd2, TSObj.decyr)

        ts_diff = ts_dd1 - ts_dd2
        coeffs_diff = np.polyfit(TSObj.decyr, ts_diff, 1)
        fit_diff = np.polyval(coeffs_diff, TSObj.decyr)
        sem2_diff = np.sqrt(sem21**2 + sem22**2)

        fig, axes = plt.subplots(figsize=(10, 6), sharex=True, nrows=3)
        axes[0].scatter(TSObj.dt, ts_dd1, color='k', s=15)
        l, = axes[0].plot(TSObj.dt, fit1, color='darkred', linestyle='--',
                          label=f'{rate1:.1f}$\\pm${sem21:.1f}')
        axes[0].fill_between(TSObj.dt, fit1-sem21, fit1+sem21, color=l.get_color(), alpha=0.3)
        axes[0].set_title(f'{locs[0]}' , fontsize=14)

        axes[1].scatter(TSObj.dt, ts_dd2, color='k', s=15)
        l, = axes[1].plot(TSObj.dt, fit2, color='darkred', linestyle='--',
                          label=f'{rate2:.1f}$\\pm${sem22:.1f}')
        axes[1].fill_between(TSObj.dt, fit2-sem22, fit2+sem22, color=l.get_color(), alpha=0.3)
        axes[1].set_title(f'{locs[1]}' , fontsize=14)
        axes[2].scatter(TSObj.dt, ts_diff, color='k', s=15)
        l, = axes[2].plot(TSObj.dt, fit_diff, color='darkred', linestyle='--',
                          label=f'{coeffs_diff[0]:.1f}$\\pm${sem2_diff:.1f}')
        axes[2].fill_between(TSObj.dt, fit_diff-sem2_diff, fit_diff+sem2_diff, color=l.get_color(), alpha=0.3)
        axes[2].set_title(f'{locs[0]} + {locs[1]}', fontsize=14)
        for ax in axes:
            ax.set_ylabel('V$_{Up}$ (mm)', fontsize=13)
            ax.grid(color='k', alpha=0.1, linestyle='--')
            ax.set_ylim([-10, 10])
            ax.legend()
        fig.subplots_adjust(hspace=0.5)
        fig.set_label(f'{locs[0]}+{locs[1]}')
        return


    def plot_tsgrad_moving(self, locs='MD_Up MD_Down'.split(), npix=0):
        from BZ.bbTS import bootstrap_1D
        TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs, npix)
        ts_dd1, coeffs_dd1, sem21 = TSObj.calc_ts_vup_dd(locs[0])
        rate1, bias1 = coeffs_dd1
        fit1 = np.polyval(coeffs_dd1, TSObj.decyr)

        ts_dd2, coeffs_dd2, sem22 = TSObj.calc_ts_vup_dd(locs[1])
        rate2, bias2 = coeffs_dd2
        fit2 = np.polyval(coeffs_dd2, TSObj.decyr)

        ts_diff = ts_dd1 - ts_dd2
        # calculate trends over shorter and shorter record (moving start date back)
        lst_rates, lst_uncs, lst_sts = [], [], []
        en = TSObj.dt[-1]
        for i in range(len(TSObj.decyr)):
            st = TSObj.dt[i]
            if (en-st).days < (365.25 * 3):
                break
            ts = ts_diff[i:]
            decyr = TSObj.decyr[i:]
            coeffs = np.polyfit(decyr, ts, 1)

            rate, unc, bias = bootstrap_1D(ts, decyr, n_boots=1500)

            lst_rates.append(rate)
            lst_uncs.append(unc)
            lst_sts.append(TSObj.dt[i])

        # filter out uncertainties; unless you show scatter no point
        df_mt = pd.DataFrame({'rate': lst_rates, 'unc': lst_uncs}, index=lst_sts)
        df_mt['1sig'] = df_mt['rate'].abs() >= df_mt['unc']
        df_mt['2sig'] = df_mt['rate'].abs() >= 2*df_mt['unc']

        df_plot = df_mt#df_mt[df_mt['2sig']]

        fig, axes = plt.subplots(figsize=(10, 4))
        # axes.scatter(lst_sts, lst_rates, color='k', s=15)
        l, = axes.plot(df_plot.index, df_plot['rate'], color='k', alpha=1.0)
        axes.fill_between(df_plot.index, df_plot['rate']-df_plot['unc'], df_plot['rate']+df_plot['unc'],
                          color=l.get_color(), alpha=0.3)
        axes.set_ylabel('Vup (mm/yr)', fontsize=14)
        axes.set_xlabel('Start Date', fontsize=14)
        axes.set_title(f'{locs[0]}-{locs[1]} Moving Trends', fontsize=16)
        axes.grid(color='k', linestyle='--', alpha=0.1)

        fig.set_label(f'{locs[0]}_{locs[1]}_moving')
        return fig, axes


    def make_transect_numbers(self, npix=0, transect='a'):
        dct_locs_a = {
            'DCDC': (38.9460665587, -77.0979088481),
            'Burleith - Hillandale': (38.91468, -77.07348),
            # 'White House': (38.8977, -77.0365),
            'Lincoln_Memorial': (38.88943, -77.04988),
            'Jefferson_Memorial': (38.881389, -77.036528),
            'DCA': (38.851774, -77.039282),
            'Potomoc_Yard': (38.828186, -77.04630),
            'Alexandria': (38.7979, -77.06332),
        }
        dct_locs_b = {
            'Ivy_City': (38.9148, -76.9859),
            'US_Capital': (38.8899, -77.0091),
            'Marine Barracks': (38.8799, -76.9943),
            'Navy Yard': ( 38.874654, -76.994709),
            'Fair Lawn': (38.86796, -76.985842),
            'SS': (38.851801, -77.015553),
            'Beltway_295': (38.793084, -77.015572),
            'National Harbor': ( 38.783015, -77.01648),
        }

        dct_locs = dct_locs_a if transect.lower() == 'a' else dct_locs_b

        log.warning(f'npix={npix}')
        for lbl, lalo in dct_locs.items():
            self.get_vel_unc_near(lalo, npix=npix, lbl=lbl)

        return


    def _add_ART(self, axes, fontsize=14, linewidth=10, color='deeppink'):
        """ Plot the Anacostia River Transect """
        art_parms = dict(color=color, transform=self.proj, ax=axes)
        gdf_tunnels = gpd.read_file(self.path_shapes / 'Tunnel_Lines.GeoJSON')
        gdf_art = gdf_tunnels.set_index('name').loc[['Anacostia_River_Tunnel']]
        gdf_art.plot(linewidth=linewidth, linestyle='-', **art_parms, edgecolor='k')

        gdf_tpoints = gpd.read_file(self.path_shapes / 'Tunnel_Points.GeoJSON')
        gdf_007 = gdf_tpoints.set_index('name').loc[['CSO 007']]
        gdf_007.plot(marker='o', markersize=45, edgecolor='w', **art_parms)

        art_parms.pop('ax')
        # txt_parms = dict(path_effects=[pe.withStroke(linewidth=1.5, foreground='w')])
        # axes.text(-76.975, 38.875, 'ART', rotation=60, fontsize=fontsize,
        #             **txt_parms, **art_parms)
        return axes


    def _add_hydrogeo(self, axes, inset=False):
        fs = 30 if not inset else 22
        lw = 2.75 if not inset else 2.25
        gdf_contours = gpd.read_file(self.path_shapes / 'Potentiometric_Surfaces.GeoJSON')
        gdf_c0 = gdf_contours[gdf_contours['surface'] == 0].set_index('year').sort_index()
        gdf_c0.index = '2011 2019'.split() # manually overwrite for now; 2007 == 2011
        # xys = [(38.8555, -76.9180), (38.8705, -76.932)]
        xys = [(38.8555, -76.9160), (38.8705, -76.932)]
        rotations = [49, 37]
        for i, color in enumerate('fuchsia m'.split()): # 2011, 2017
            gdfi = gdf_c0.iloc[[i]]
            gdfi.plot(ax=axes, transform=self.proj, color=color, zorder=30, linewidth=lw)
            axes.text(xys[i][1], xys[i][0], gdfi.index.item(),
                    fontsize=fs, color=color, zorder=30,
                    transform=self.proj, rotation=rotations[i],
                    path_effects=[pe.withStroke(linewidth=1, foreground='k')])


        ## add the faults
        gdf_faults = gpd.read_file(self.path_shapes / 'Possible_Faults.GeoJSON')
        # gdf_faults = gdf_faults[gdf_faults['name'].isin('NE_DFZ_Long NE_BFZ_Long'.split())]
        # yxs = [(38.8360, -77.0780), (38.8190, -77.0425)]
        gdf_faults = gdf_faults[gdf_faults['name'].isin('NE_DFZ NE_BFZ'.split())]
        yxs = [(38.8350, -77.0780), (38.8200, -77.0425)]
        rotations = [55, 55]
        for i, color in enumerate('sandybrown darkorange'.split()):
            color = 'darkorange'
            gdfi = gdf_faults.iloc[[i]]
            name = gdfi['name'].item().split('_')[1]
            gdfi.plot(ax=axes, transform=self.proj, color=color, zorder=30,
                      linewidth=lw, linestyle='--')
            axes.text(yxs[i][1], yxs[i][0], name,
                    fontsize=fs, color=color, zorder=30,
                    transform=self.proj, rotation=rotations[i],
                    path_effects=[pe.withStroke(linewidth=1, foreground='k')])
        return


    def plot_seasonal_phase(self, mask_amp=True):
        assert self.neofs == 0, 'Use neofs=0 to make sure phase is in there'
        # self.pparms['cmap'] = 'cmo.phase'
        self.pparms['cmap'] = 'cmo.curl'
        self.pparms['norm'] = mpl.colors.TwoSlopeNorm(0, -3.14, 3.14)
        self.show_gps=False
        arr, meta= readfile.read(self.path_vup_geo, datasetName='annualPhase')
        da_seasonal = self.da.copy()
        da_seasonal.data = np.flipud(arr) if meta['ORBIT_DIRECTION'] == 'ASCENDING' else arr

        arr, meta= readfile.read(self.path_vup_geo, datasetName='annualAmplitude')
        da_amp = self.da.copy()
        da_amp.data = np.flipud(arr) if meta['ORBIT_DIRECTION'] == 'ASCENDING' else arr

        self.da = da_seasonal.where(self.da.notnull(), np.nan)

        # self.da = xr.where(da_amp>=da_amp.mean(), self.da, np.nan) if mask_amp else self.da
        self.da = xr.where(da_amp>=da_amp.quantile(0.25), self.da, np.nan) if mask_amp else self.da

        f, a = self.plot_basic_DC()
        # a.set_extent([-77.075, -76.85, 38.75, 38.920])
        self._add_hydrogeo(a)
        im = a.get_children()[1]
        cbar = im.colorbar
        cbar.set_label('Annual Phase (radians)', rotation=270, labelpad=0, fontsize=18)
        cbar.set_ticks([-3.14, 3.14], labels='$-\\pi$ $\\pi$'.split(), fontsize=20)

        f.set_label('')
        f.set_label(f'SeasonalPhase_{self.reg}')
        return


def plot_gps_insar():
    from VLM.bzFRInGE.analysis import analyzeGPS
    Exp0   = ExpBase(DC_SR, 'ERA5_SET_PM_ex_Fast', 'USN7', neofs=20)
    GPSObj = analyzeGPS.GPS_Rate_Analysis(Exp0)
    GPSObj.plot_gps_vs_insar()
    GPSObj.plot_gps_vs_insar_qq()


def plot_gradient_map_ts(ExpI, buff=0):
    """ Plot the TWO maps and gradient ts and arrange them in keynote

    buff = 0 will plot points with atmosphere removed
    buff >0 will average in box around the point
    """
    kind = 'unc'
    PObjA = PlotDC(ExpI, kind, sten='_20161229-20200822')

    ## so no gridlines
    fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': PObjA.basemap.crs})
    ax0 = PObjA.plot_gradient_map(buff, show_box=False, lbl='A', axes=axes)[1]
    im0 = ax0.get_children()[0]
    im0.colorbar.set_label('mm/yr', labelpad=-40, fontsize=24)
    im0.colorbar.ax.tick_params(labelsize=20)
    # im0.colorbar.remove() ## use this to make the maps the same size

    PObjB = PlotDC(ExpI, kind, sten='_20200822-20240427')
    fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': PObjB.basemap.crs})
    ax1 = PObjB.plot_gradient_map(buff, lbl='B', axes=axes)[1]
    im1 = ax1.get_children()[0]
    im1.colorbar.remove()

    PObjC = PlotDC(ExpI, kind='rate') # full for gradient ts
    if buff:
        ax2 = PObjC.plot_gradient_ts(buff, show_each=True)[1]
    else:
        ax2 = PObjC.plot_gradient_ts_atm(show_each=False)[1]


##------------------------------------------------------------------ DEPRECATED
def plot_gradient_maps3(kind='rate'):
    """ Plot three 3 year gradient maps """
    sts = '20150310 20180318 20210326'.split()
    ens = '20180318 20210326 20240427'.split()
    for i in range(3):
        sten = f'{sts[i]}_{ens[i]}'
        nes = 15 if i < 2 else 10
        sb = False if i < 2 else True
        Expi = ExpBase(DC_SR, f'ERA5_SET_PM_ex_Fast_{sten}', 'USN8', neofs=nes)
        PObj = PlotDC(Expi, kind=kind)
        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': PObj.basemap.crs})
        ax1 = PObj.plot_gradient_map(buff=500, show_box=sb, axes=axes,
                                     lbl=string.ascii_uppercase[i])[1]
        im1 = ax1.get_children()[0]
        im1.colorbar.remove()
        st = pd.to_datetime(sts[i]).strftime('%b, %Y')
        en = pd.to_datetime(ens[i]).strftime('%b, %Y')
        ax1.set_title(f'{st} to {en}', fontsize=17)
    return


def gridspec_gradient_map_ts():
    """ Dont do this to much work; arrange it manually in keynote """
    Exp0 = ExpBase(DC_SR, 'ERA5_SET_PM_ex_Fast', 'USN8', neofs=15)
    PObj0 = PlotDC(Exp0, kind='rate')

    Exp1 = ExpBase(DC_SR, 'ERA5_SET_PM_ex_Fast_2020', 'USN8', neofs=15)
    PObj1 = PlotDC(Exp0, kind='rate')

    fig = plt.figure(figsize=(16, 12))

    gs = GridSpec(2, 2, height_ratios=[1, 0.5], width_ratios=[1, 1])
    axes = []
    axes.append(plt.subplot(gs[0, 0], projection=PObj.basemap.crs))
    axes.append(plt.subplot(gs[0, 1], projection=PObj.basemap.crs))

    ax0 = PObj0.plot_gradient_map(buff=500, axes=axes[0], show_box=False, lbl='A')[1]
    im = ax0.get_children()[0]
    im.colorbar.set_label('')

    ax1 = PObj1.plot_gradient_map(buff=500, axes=axes[1], lbl='B')[1]
    im = ax1.get_children()[0]
    im.colorbar.remove()

    axes.append(plt.subplot(gs[1, :2]))
    ax2 = PObj.plot_gradient_ts(buff=500, show_each=False, axes=axes[-1])[1]
    fig.set_label('DC_Gradient_Map_TS')


# --------------------------------------------------------------------------- #
def plot_gradient_insets(exp0):
    """ Wrapper for many gradient inset maps """
    PObj = PlotDC(exp0, 'rate', show_gps=True)
    gdf_pts = gpd.read_file(PObj.path_wd / 'shapes' / 'MD_Points1.GeoJSON', index_col=0)
    for i in range(1, 1+int(gdf_pts.shape[0]/2)):
        # if i!=3:
        #     continue
        ## add lat lon to points dct
        for loci in 'MD_Up MD_Down'.split():
            loc = f'{loci}{i}'
            pt = gdf_pts.set_index('name').loc[loc].geometry
            DCT_PTS[loc] = [(pt.y, pt.x), 'DC']
        PObj.plot_inset_grad(f'MD_Up{i} MD_Down{i}'.split(), double_diff=True, npix=1)
    return


def make_all_figs(ExpObj):
    PObj = PlotDC(ExpObj, kind='rate', show_gps=True)
    PObj.plot_basic_DC()
    PObj.plot_transect_map()
    plot_gradient_map_ts(ExpObj)
    plot_gradient_map_moving(ExpObj)

    PObj.plot_insets('Jefferson_Memorial Anacostia'.split(), npix=1)
    # PObj.plot_inset('Alexandria', show_stable=False)


def make_crop():
    """ Make the crop of the InSAR """
    cities = ["Washington DC", 'Alexandria, VA', 'Arlington, VA']#, 'Fairfax County, VA']
    lst_polys = [city_shape(city, 0)[0] for city in cities]
    gdf_cities = pd.concat(lst_polys).set_index('name')
    return gdf_cities

if __name__ == '__main__':
    # Exp0 = ExpBase(DC_SR, 'Base_ex_Fast_2017', 'USN8', neofs=0)
    Exp0 = ExpBase(DC_SR, 'ERA5_SET_PM_ex_Fast_2017', 'USN8', neofs=20)
    # make_numbers(Exp0, PlotDC, qH=0.997)
    # make_all_figs(Exp0)

    PObj = PlotDC(Exp0, 'rate', show_gps=True)
    ## main figures ----------------------------------------------------
    PObj.plot_basic_DC()
    # PObj.plot_transect_map()
    # PObj.plot_inset_grad(double_diff=True)
    # plot_gradient_map_ts(Exp0, 0)

    # PObj.plot_seasonal_phase()

    ## supplement figures ----------------------------------------------------
    # plot_gradient_insets(Exp0)

    # PObj = PlotDC(Exp0, kind='unc')
    # PObj.plot_basic_DC(False, False, False, False)
    # PObj.make_transect_numbers(npix=1, transect='A') ## Subsidence
    # PObj.make_transect_numbers(npix=1, transect='B') ## Uplift
    # PObj.plot_insets()
    # PObj.plot_inset('Alexandria', npix=0)
    # PObj.plot_city('Ivy_City', loc=1, alpha=1, crop=True)
    # PObj.plot_city('Capitol Hill, Washington DC', alpha=1, crop=True)

    # da_tay = get_tay22('DC')
    # PObj.plot_basic(da=da_tay)

    # PObj.plot_inset('Jefferson_Memorial', npix=1)

    path_figs = op.join(PATH_RES, f'{Exp0.reg}_2024')
    # bbPlot.savefigs(path_figs, True, True, dpi=300)
    plt.show()

