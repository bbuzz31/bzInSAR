from VLM.bzFRInGE.plotting.plotExp import *
from VLM.bzFRInGE.plotting.plot_TS1 import PlotTS1, DCT_PTS
# from VLM.bzFRInGE.plotting.plotTS import plot_ts_vup1 ## see plot DC
from matplotlib.gridspec import GridSpec#, GridSpecFromSubplotSpec
from mintpy.utils import utils as ut
import geopandas as gpd

DCT_OILN = {}
DCT_OILS = {}


class PlotHouston(PlotExp):
    def __init__(self, exp, kind='rate', show_gps=True, continuous=True):
        super().__init__(exp, kind, show_gps, continuous)
        self._set_plot_parms()
        self.gdf_facs = self.get_oil_locs()
        self.show_shipc = True
        self.show_counties = True


    def _set_plot_parms(self):
        super()._set_plot_parms()
        ## Plot Parms
        self.basemap = cimgt.GoogleTiles(
            # url='https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}')
            # url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
            url='https://server.arcgisonline.com/arcgis/rest/services/Canvas/'\
            'World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}.jpg')

        i = 0 if self.continuous else 1
        if self.kind == 'Rate':
            self.pparms['cmap'] = 'cmo.ice'
        elif self.kind == 'Tvar':
            self.show_gps = False
        # elif self.kind == 'Unc':
        #     self.pparms['norm'] = mpl.colors.Normalize(0.5, 1.5)

        return


    def plot_basic0(self, da=None, axes=None):
        self.show_gps = False
        fig, axes = super().plot_basic(da=da, axes=axes)
        if axes is None: # there may already be one
            im  = axes.get_children()[1]
            im.colorbar.remove()
            bbPlot.cartopy_cbar(im, ylabel=self.cbar_lbl, fontsize=16, pad=0.15)

        self.show_gps = True
        return fig, axes

    def plot_basic(self, da=None, show_pts=True, axes=None):
        self.show_gps_names = False
        if self.kind == 'Tvar':
            self.show_gps = False

        fig, axes = super().plot_basic(da=da, axes=axes)
        if axes is None: # there may already be one
            im  = axes.get_children()[1]
            im.colorbar.remove()
            bbPlot.cartopy_cbar(im, ylabel=self.cbar_lbl, fontsize=16, pad=0.15, extend='min')

        ## add the pts if wanted
        pts = 'BBayou Texas_City'.split()
        pts.append('Nonlinear') if self.kind == 'Tvar' else ''
        pts = [] if not show_pts else pts
        color, edgecolor = 'gold', 'white' #edgecolor=white normally, gold for nonlinear
        from VLM.bzFRInGE.plotting.plot_TS1 import DCT_PTS
        sty_box = dict(zorder=50, alpha=1, facecolor='none', edgecolor=edgecolor,
                    linestyle='--', linewidth=1.5, path_effects=[pe.withStroke(linewidth=3.0, foreground=color)])
        sty_pt = dict(path_effects=[pe.withStroke(linewidth=2.5, foreground=color)],
                        color=edgecolor, transform=ccrs.PlateCarree(), s=30, marker='^', facecolor='none')

        for i, k in enumerate(pts):
            lalo = DCT_PTS[k][0]
            SNWE = buffer_point(*lalo, 2000)
            bbox = bbGIS.bbox2poly(SNWE=SNWE)

            ## add triangle in th middle and box around it
            # axes.scatter(lalo[1], lalo[0], **sty_pt)
            # axes.add_geometries([bbox], crs=self.proj, **sty_box)

        ## add the ship channel
        gdf_shipc = gpd.read_file(self.path_wd / 'Houston_Ship_Channel.GeoJSON').dropna(subset=['geometry'])
        sty_shipc = dict(facecolor='none', edgecolor='white', alpha=0.9, linestyle='--',
                         path_effects=[pe.withStroke(linewidth=2.5, foreground='red')],
                         zorder=40, crs=self.proj)
        # geopandas causes second basemap to not show
        # gdf_shipc.plot(ax=axes, **sty_shipc) if self.show_shipc else ''
        # axes.add_geometries(gdf_shipc.geometry, **sty_shipc)

        # lon , lat rotation
        dct_counties = {'Harris': [-95.84, 29.925, 0], 'Galveston': [-95.06, 29.2375, 0],
                        'Chambers': [-94.8343, 29.826, 310], 'Fort Bend': [-95.8435, 29.60, 0],
                        'Brazoria': [-95.40, 29.3, 0], 'Liberty': [-94.84, 30.0, 0]}
        dct_counties = {} if not self.show_counties else dct_counties
        sty_cnty  = dict(zorder=40, alpha=0.9, facecolor='none', edgecolor='white',
                    linestyle=':', linewidth=1.25)
        sty_name = dict(facecolor='lightgray', alpha=0.5, boxstyle='round')
        for county, lst in dct_counties.items():
            if not county in 'Galveston Liberty'.split() and self.show_counties:
                gser = self.get_county_gser(county)
                # geom = gser.geometry if not county == 'Galveston' else gser.geometry.convex_hull
                axes.add_geometries(gser.geometry, crs=self.proj, **sty_cnty)

            if self.kind == 'Rate' and show_pts:
                county = 'Fort\nBend' if county == 'Fort Bend' else county
                x, y, rot = lst
                # y, halign = 0, 'left'
                geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=y)
                axes.text(x, y, county, transform=text_transform, zorder=20, rotation=rot,
                        verticalalignment='center', horizontalalignment='left',
                            bbox=sty_name, size=11.5, color='black')

        ## inset GPS plot
        add_inset = False
        if add_inset:
            # fig.add_axes arguments:
            # [left, bottom, width, height] - all values between 0-1 representing fraction of figure dimensions
            inset_ax = fig.add_axes([0.2, 0.2, 0.3, 0.2])
            GPSAnalysis = analyzeGPS.GPS_Rate_Analysis(self)
            GPSAnalysis.plot_gps_vs_insar_qq(npix=1, llim=-20, ulim=5, axes=inset_ax)
            fig.set_label(f'{fig.get_label()}_inset')

        # bbPlot._scale_bar(axes, 25, (0.13, 0.015))
        bbPlot._scale_bar(axes, 25, (0.87, 0.015))

        return fig, axes


    def plot_nonlinear_mask(self, thresh=None):
        assert self.kind == 'Tvar', 'Use temporal variabilty'
        if thresh is None:
            thresh = self.da.mean() + 2*self.da.std()
            print (f'Calculating thresh using 2-sigma range: {thresh:.2f}')
        da_mask = xr.where(self.da < thresh, 1, 0) # 1 is linear
        da_mask = self.da.where(self.da.isnull(), da_mask)
        n_valid = da_mask.notnull().sum()
        n_linear = da_mask.sum()
        log.info (f'Percent Linear: {100*n_linear/n_valid:.2f}%')


        cmap = mpl.colors.ListedColormap('maroon white'.split())
        bounds = np.arange(-0.5, 2.5, 1)
        labels = 'Nonlinear Linear'.split()
        centers = bounds[:-1] + (np.diff(bounds)/2)
        norm  = mpl.colors.BoundaryNorm(bounds, cmap.N)

        self.pparms['cmap'] = cmap
        self.pparms['norm'] = norm

        fig, axes = self.plot_basic(da=da_mask)
        cbar = axes.get_children()[1].colorbar
        cbar.set_ticks(centers, labels=labels, rotation=-90)
        cbar.ax.tick_params(length=0)
        cbar.set_label('')
        fig.set_label('Nonlinearity_Mask')
        return fig, axes


    def crop_county(self, county='Harris'):
        pass


    def plot_inset(self, loc='BBayou', npix=1, double_diff=True, show_fit=True):
        """ Plot one point, row wise """
        log.warning('Not double differencing...') if not double_diff else ''
        self.show_gps_names = False
        self.show_gps = False
        fig = plt.figure(figsize=(18, 10), constrained_layout=True)
        gs  = GridSpec(2, 3, width_ratios=[1.3, 1, 1], figure=fig)
        axes = []
        axes.append(plt.subplot(gs[0, 0], projection=self.basemap.crs))
        axes.append(plt.subplot(gs[0, 1:]))

        ## Plot the Map
        lalo = DCT_PTS[loc][0] #  use plot_ts to make a box around the point
        S, N, W, E = buffer_point(*lalo, 2000)
        axes[0].set_extent([W, E, S, N], crs=self.proj)
        self.zoom  = 15

        ax  = self.plot_basic(axes=axes[0], show_pts=False)[1]
        im  = ax.get_children()[0]
        im.colorbar.remove()
        im.set_alpha(0.75)

        ## add actual point
        pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                        color='w', transform=ccrs.PlateCarree(), s=50)
        ax.scatter(lalo[1], lalo[0], **pt_parms)
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
            axeTS, df_ts = TSObj.plot_ts_vup_dd(loc, axes=axes[1], show_fit=show_fit)[1:]
        else:
            axeTS, df_ts = TSObj.plot_ts_vup(loc, axes=axes[1], show_fit=show_fit)[1:]
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

        da_tvar = self.get_da('tvar')
        da_tvar1 = bbGIS.get_surrounding_pixels(da_tvar, lalo[0], lalo[1], npix=npix)
        tvar = da_tvar1.mean().item()
        if show_fit:
            axeTS.fill_between(df_ts.index, df_ts['trend']-tvar, df_ts['trend']+tvar,
                                color='magenta', alpha=0.7, facecolor='none', linestyle='--')
        log.info (f'Temporal Variability of {loc}: {tvar:.2f}')

        fig.set_label(f'Inset_{loc}')
        return


    def plot_insets(self, locs='BBayou Texas_City'.split(), double_diff=True, npix=1):
        """ Plot two points, row wise """
        # dct_ylims = {'Backliff': [-11, 6], 'BBayou': [-1, -20]}
        dct_ylims = {}#'Backliff': [-11, 6], 'BBayou': [-1, -20]}

        log.warning('Not double differencing...') if not double_diff else ''
        self.show_gps = False
        self.zoom  = 15
        fig = plt.figure(figsize=(14, 9), constrained_layout=False)
        gs  = GridSpec(2, 3, width_ratios=[1, 1, 1], hspace=0.33, wspace=0.0, figure=fig)
        axeMaps = []
        axeTSs  = []
        axeMaps.append(plt.subplot(gs[0, 0], projection=self.basemap.crs))
        axeMaps.append(plt.subplot(gs[1, 0], projection=self.basemap.crs))
        axeTSs.append(plt.subplot(gs[0, 1:]))
        axeTSs.append(plt.subplot(gs[1, 1:]))
        for i, loc in enumerate(locs):
            lalo = DCT_PTS[loc][0] #  use plot_ts to make a box around the point
            S, N, W, E = buffer_point(*lalo, 2000)
            axeMaps[i].set_extent([W, E, S, N])

            ax  = self.plot_basic(axes=axeMaps[i], show_pts=False)[1]
            # ax.text(-0.01, 1.01, string.ascii_uppercase[i], transform=ax.transAxes,
            ax.text(1.0, 0.999, string.ascii_uppercase[i+1], transform=ax.transAxes,
                    verticalalignment='top', fontsize=20, color='w',
                    bbox=dict(boxstyle='square,pad=0.05', facecolor='black'))
            im  = ax.get_children()[0]
            im.set_alpha(0.75)
            im.colorbar.remove()
            if loc == 'Texas_City':
                bbPlot._scale_bar(ax, 1, (0.835, 0.015))
            else:
                bbPlot._scale_bar(ax, 1, (0.20, 0.015))

            # put the point in the middle
            pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                            edgecolor='w', facecolor='none', transform=ccrs.PlateCarree(), s=30)
            ax.scatter(lalo[1], lalo[0], **pt_parms)

            # gl = ax.gridlines(draw_labels=True)#, xlocs=np.arange(-180, 190, 0.01),
                        # ylocs=np.arange(-90, 100, 0.01))

            # bbPlot.fmt_gridlines(gl, left=True, right=False, bottom=True, size=GS-1)
            axeMaps[i].add_image(self.basemap, self.zoom) # it's getting removed in the second one
            TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs, npix)
            show_fit = False if loc == 'Nonlinear' else True
            if double_diff and not loc == 'Nonlinear':
                axeTS, df_ts = TSObj.plot_ts_vup_dd(loc, radius=2.5e3, corr_thresh=0.95,
                                                       axes=axeTSs[i], show_fit=show_fit)[1:]
                axeTS.set_ylim(dct_ylims.get(loc, (-30, 30)))
            else:
                axeTS, df_ts = TSObj.plot_ts_vup(loc, axes=axeTSs[i], show_fit=show_fit)[1:]

            # pos = axeTS.get_position().bounds
            # pos[0] = pos[0]-0.1
            # axeTS.set_position(mpl.transforms.Bbox.from_bounds(*pos))

            axeTS.set_title('')
            ## set ytick lims
            axeTS.tick_params(labelbottom=False) if i == 0 else ''
            axeTS.tick_params(axis='both', which='major', labelsize=13)

            if i == 1:
                axeTS.set_ylabel(axeTS.get_ylabel(), labelpad=12)

            axeTS.set_ylabel('Vertical Displacment (mm)', fontsize=13,
                             rotation=270, labelpad=20)
            axeTS.yaxis.set_label_position('right')
            axeTS.yaxis.set_ticks_position('right')

            ## add the temporal variability
            da_tvar = self.get_da('tvar')
            da_tvar1 = bbGIS.get_surrounding_pixels(da_tvar, lalo[0], lalo[1], npix=npix)
            tvar = da_tvar1.mean().item()
            axeTS.fill_between(df_ts.index, df_ts['trend']-tvar, df_ts['trend']+tvar,
                               color='magenta', alpha=0.7, facecolor='none', linestyle='--')
            log.info (f'Temporal Variability of {loc}: {tvar:.2f}')

        # make a single colorbar across
        cbar_ax = fig.add_axes([0.165, 0.475, 0.160, 0.0160]) # left, bot, wid, hgt
        cbar_ax.set_title('V$_{Up}$  (mm yr$^{-1}$)', fontsize=11)#, labelpad=1v#, rotation=270, )
        # cbar_ax.yaxis.set_label_position('right')
        cbar    = plt.colorbar(im, cax=cbar_ax, orientation='horizontal', spacing='uniform', extend='both')
        # gs.update(left=0.07, right=0.23)#, bottom=0.1, top=0.9)
        # gs.update(hspace=-0.5)#, wspace=10.00)
        fig.set_label(f'Insets_{"_".join(locs)}')
        fig.constrained_layout=True
        return


    def plot_transect_map(self):
        self.gps_s = 100 # increase GPS sizes
        self.show_gps_names = False

        transects = 'Shipping_Channel_South NE_Transect'.split()
        fig, axes = super().plot_transect_map(transects, color='orange', fontsize=24)
        axes.set_extent([-95.4626, -94.8630, 29.4932, 29.8890])
        im  = axes.get_children()[1]
        # im.colorbar.remove()
        gl = axes.gridlines(draw_labels=True)#, xlocs=np.arange(-180, 190, 0.01),
                    # ylocs=np.arange(-90, 100, 0.01))

        bbPlot.fmt_gridlines(gl, left=True, right=False, bottom=True, size=GS-1)
        return fig, axes


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


    def get_oil_locs(self):
        df_facs = pd.read_csv(op.join(self.path_wd, 'Houston_Facilities.txt'),
                    sep='\\s+', header=None, index_col=0, names='lat lon'.split())
        gdf_facs = bbGIS.df2gdf(df_facs)

        #get the closest reference GPS station; for stitching not really needed
        ref_stas = 'NASA UHDT CFJV SESG TXTG NASA CSTE DMFB TXB6 TXLQ'.split()
        gdf_refs  = bbGIS.df2gdf(self.df_gps[self.df_gps.sta.isin(ref_stas)].set_index('sta'), 'all')
        gdf_facs = gdf_facs.to_crs(4087).sjoin_nearest(gdf_refs.to_crs(4087), distance_col='dist').to_crs(4326)
        gdf_facs['dist'] /= 1000
        return gdf_facs


    def update_unc(self):
        """ Use Marins Uncertainty """
        from argparse import Namespace
        from mintpy.save_gdal import save_gdal
        from mintpy.utils import writefile
        assert None, 'Would rather not overwrite file.'
        path_vup_geo = Path(self.path_vup_geo)
        path_vup_corr = path_vup_geo.parent / 'geo_velocity_recon10corr.h5'
        unc = readfile.read(path_vup_corr, datasetName='velocityStdCorr')[0]
        inc  = readfile.read(self.path_geom_mp_geo, datasetName='incidenceAngle')[0]

        ser_ref  = self.df_gps[self.df_gps.sta == self.ref_sta]
        unc_ref  = ser_ref.u_sig.item()/1000
        unc  /= np.cos(np.deg2rad(inc)) # bumps it up slightly
        unc   = np.sqrt(unc**2 + unc_ref**2)

        print (f'Tieing to: {ser_ref.sta.item()}')
        writefile.write_hdf5_block(path_vup_geo, unc, 'velocityStd')
        inps= Namespace(file=self.path_vup_geo, dset='velocityStd',
                        outfile=self.path_std_nc, out_format='netCDF')
        save_gdal(inps)

        ## apply the mask to rate and std
        mask = xr.open_dataset(self.path_mask_vup_nc)['Band1']
        mask1 = mask.where(~np.isclose(mask, 0), np.nan)

        unc   = xr.load_dataset(self.path_std_nc)['Band1']

        unc_m   = unc  * mask1

        units    = 'm/y'
        unc_m.attrs = {'units': units, 'region': self.reg}

        # write the masked versions for plotting and transects
        try:
            os.remove(self.path_std_msk_nc)
        except:
            pass

        unc_m.to_netcdf(self.path_std_msk_nc)

        print (f'Wrote corrected Std to {self.path_std_msk_nc}, {datetime.now()}.')
        print ('Finished @:', datetime.now().strftime('%m/%d/%Y %H:%M:%S'))
        # del mask, rate, unc, rate_m, unc_m
        return


    def get_county_gser(self, county):
        """ Search texas counties to find the shapefiles """
        src = self.path_wd / 'Texas_County_Boundaries.GeoJSON'
        assert src.exists(), f'{src} doesnt exist!'
        gdf_counties = gpd.read_file(src)
        assert county in gdf_counties['CNTY_NM'].tolist(), f'{county} not in Texas counties'
        gser_county = gdf_counties[gdf_counties['CNTY_NM'] == county]
        return gser_county


def scatter_facilities_vlm(ExpI, npix=0):
    """ Plot VLM at the facilities on a map """
    import seaborn as sns
    ExpObj = PlotHouston(ExpI)
    gdf_locs = ExpObj.gdf_facs
    df_res = vlm_at_points(ExpObj, gdf_locs, npix).rename_axis('name')
    fig, axes = plt.subplots(figsize=(10, 6))
    axes.errorbar(np.arange(df_res.shape[0]), df_res['ins_vel'],
                            yerr=df_res['ins_unc'], color='k', fmt='o')

    # sns.histplot(df_res, x='name', y='ins_vel', kde=True, color='dimgray', ax=axes)
    axes.set_ylabel('VLM (mm/yr)', fontsize=14)
    axes.set_xlabel('Facility', fontsize=14)
    axes.grid(color='k', linestyle='--', alpha=0.1)
    fig.set_label('vlm_facilities_scatter')
    return df_res, fig, axes


def plot_kde_facilities_vlm(ExpI, npix=0):
    """ Plot total histogram and just the vlm at the facilities on a map """
    import seaborn as sns
    ExpObj = PlotHouston(ExpI)
    df_vlm = ExpObj.da.to_dataframe().reset_index().dropna()

    gdf_locs = ExpObj.gdf_facs
    df_res = vlm_at_points(ExpObj, gdf_locs, npix).rename_axis('name')

    lst_sers = df_vlm['Rate'], df_res['ins_vel']

    fig, axes = plt.subplots(figsize=(10, 10), ncols=2, sharey=True, sharex=True)
    for i, (ax, ser) in enumerate(zip(axes, lst_sers)):
        # m, l, h = ser.median(), ser.quantile(0.003), ser.quantile(0.997)
        m, l, h = ser.median(), ser.quantile(0.1), ser.quantile(0.9)
        sns.kdeplot(ser, ax=ax, color='dimgray', fill=True)
        # ax.axvline(m, color='darkred', label=f'{m:.2f} +/- {(l+h)/2:.2f}')
        # ax.axvline(l, linestyle='--', color='darkred', alpha=0.75)
        # ax.axvline(h, linestyle='--',color='darkred', alpha=0.75)
        # ax.legend()
        ax.grid(color='k', linestyle='--', alpha=0.1)
        ax.tick_params(axis='both', labelsize=12)

    axes[0].set_xlabel('VLM (mm/yr)', labelpad=10, fontsize=14)
    axes[0].set_ylabel('Density', labelpad=10, fontsize=14)
    axes[1].set_xlabel('Facility VLM (mm/yr)', labelpad=10, fontsize=14)
    # axes[1].set_ylabel('')

    fig.set_label('vlm_facilities_kde')
    return fig, axes


def plot_insar_oil_box():
    """ Plot the insar locations, oil facs, and the boxes of interest """
    Exp0 = ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=10)
    PObj = PlotHouston(Exp0, kind='rate')
    gdf_facs = PObj.gdf_facs

    # SNWE
    box1 = [29.69, 29.89, -95.31, -94.855]
    box2 = [29.35, 29.475, -95.075, -94.865]

    gdf_facs_b1 = bbGIS.in_dfll(gdf_facs, SNWE=box1)
    gdf_facs_b2 = bbGIS.in_dfll(gdf_facs, SNWE=box2)

    f, a = PObj.plot_basic(da=None)
    PObj.add_pts(a, gdf_facs)
    # PObj.add_bbox(a, box1)
    PObj.add_bbox(a, box2)
    # bbPlot._scale_bar(a, 25, (0.13, 0.015))

    # bbPlot.savefigs(path_figs, True, True, dpi=300)


def plot_oil_ts():
    Exp0 = ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=10)
    PObj = PlotHouston(Exp0, kind='rate')
    box1 = [29.69, 29.89, -95.31, -94.855]
    box2 = [29.35, 29.475, -95.075, -94.865]
    gdf_facs = PObj.gdf_facs.dropna()

    gdf_facs_b1 = bbGIS.in_dfll(gdf_facs, SNWE=box1)
    gdf_facs_b2 = bbGIS.in_dfll(gdf_facs, SNWE=box2)

    ## plot the rest of the facs
    ix_done    = gdf_facs_b1.index.tolist() + gdf_facs_b2.index.tolist()
    gdf_facs3 = gdf_facs[~gdf_facs.index.isin(ix_done)]

    gdf_facs_plot = gdf_facs_b2

    for i, (sta, ser) in enumerate(gdf_facs_plot.iterrows()):
        plot_ts_vup1(Exp0, sta, (ser.lalo.y, ser.lalo.x), npix=0, save=False)
        print (f'\nFinished {i+1} of {len(gdf_facs_plot)}')


def plot_oil_ts_dd_one(sta, npix=0, thresh=0.5, save=False):
    Exp0  = ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=10)
    PObj = PlotHouston(Exp0, kind='rate')
    gdf_facs0 = PObj.gdf_facs # target location
    ser0       = gdf_facs0.loc[f'R6-TX-{sta}']
    loc_target = [sta, (ser0.lalo.y, ser0.lalo.x)]

    gdf_facs  = gpd.read_file(Path(Exp0.path_wd) / 'stable_oil.GeoJSON') # stable
    ser        = gdf_facs[gdf_facs.name == sta]
    ser_keep   = ser[ser.velocity <= thresh]

    if ser_keep.empty:
        log.error(f'R6-TX-{sta} has no stable points!')
        return

    ser1 = ser_keep.sort_values('velocity').iloc[0]
    loc_stable = [f'{sta}_stable', (ser1.geometry.y, ser1.geometry.x)]

    f, a = plot_ts_vup1_dd(Exp0, loc_target, loc_stable, npix=npix, show_ts=False)

    if save:
        path_figs = op.join(PATH_RES, f'{Exp0.reg}_2024', 'Figures', 'double_diffs')
        bbPlot.savefigs(path_figs, figures=False)
    return f, a


def plot_oil_ts_dd_all(thresh=0.5, npix=0):
    Exp0  = ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=10)
    PObj = PlotHouston(Exp0, kind='rate')
    gdf_facs0 = PObj.gdf_facs # target location
    gdf_facs  = gpd.read_file(Path(Exp0.path_wd) / 'stable_oil.GeoJSON') # stable

    for i, (sta, ser) in enumerate(gdf_facs.groupby('name')):
        ser_keep = ser[ser.velocity <= thresh]
        if ser_keep.empty:
            log.warning(f'R6-TX-{sta} has no stable points.')
            continue
        ser1 = ser_keep.sort_values('velocity').iloc[0]

        ## setup target location
        ser0 = gdf_facs0.loc[f'R6-TX-{sta}']
        dct_pts[sta]    =  [(ser0.lalo.y, ser0.lalo.x), Exp0.reg]

        ## setup stable location
        dct_stable[Exp0.mp_exp0] =  {f'{sta}_stable': [(ser1.geometry.y, ser1.geometry.x), Exp0.reg]}

        try:
            plot_ts_vup1_dd(Exp0, sta, f'{sta}_stable', npix=npix, show_ts=False)
        except Exception as E:
            print (E)
            log.warning(f'{sta} failed, continuing')

    path_figs = op.join(PATH_RES, f'{Exp0.reg}_2024', 'double_diffs')
    bbPlot.savefigs(path_figs)
    return


def plot_gps_insar(ExpI, npix=1):
    from VLM.bzFRInGE.analysis import analyzeGPS
    llim = -18
    ulim = 0

    GPSObj = analyzeGPS.GPS_Rate_Analysis(ExpI)
    fig, axes = GPSObj.plot_gps_vs_insar(npix=npix)
    axes.tick_params(axis='x', rotation=90)
    axes.set_ylim([llim, ulim])
    GPSObj.plot_gps_vs_insar_qq(npix=npix, llim=llim, ulim=ulim)


def plot_scenarios(wd, axes=None):
    """ CSV of IPCC from Marin; adjusted to intermediate low scenario:

    Galveson Pier 21 IPCC 2-4.5 in 2050 is 47.4 cm; Scenario is 50 cm
    Note these are relative to 2000 baseline, the web numbers are relative to today (2020)
    """
    # wd = Path.home() / 'Downloads'
    locs = 'BBAY TEXAS_CITY'.split()
    lbls = 'Buffalo Bayou', 'Texas City'
    colors = 'purple pink'.split()
    if axes is None:
        fig, axes = plt.subplots(figsize=(10, 4))
    else:
        fig = axes.get_figure()

    df0 = pd.read_csv(wd / 'GALVESTON__psml828_tg.csv', index_col=0)
    df0 = df0[((df0.index>2019) & (df0.index<2051))] / 10
    df0 += 2.6
    l, = axes.plot(df0.index, df0[f'IPCC-50'], label='Galveston Pier (Int. Low)', color='k')#, zorder=50)
    axes.plot(df0.index, df0[f'IPCC-17'], color=l.get_color(), linestyle='--', alpha=0.8)#, zorder=50)
    axes.plot(df0.index, df0[f'IPCC-83'], color=l.get_color(), linestyle='--', alpha=0.8)#, zorder=50)

    for i, loc in enumerate(locs):
        df = pd.read_csv(wd / f'{loc}.csv', index_col=0)
        df = df[((df.index>2019) & (df.index<2051))]/10
        df += 2.6 # Projection to scenario
        ext = 'wVLM' #if i > 0 else ''
        l, = axes.plot(df.index, df[f'IPCC-50{ext}'], label=f'{lbls[i]}', color=colors[i])
        axes.fill_between(df.index, df[f'IPCC-17{ext}'], df[f'IPCC-83{ext}'],
                          color=l.get_color(), linestyle='-', alpha=0.3)
    axes.set_xlim([2019, 2050])
    axes.legend(loc='upper left', fontsize=11)
    axes.set_ylabel('RSL (cm)', fontsize=13)
    axes.grid(color='k', linestyle='--', alpha=0.1)
    fig.set_label('Projections')
    return fig, axes, df0


## what is this for?
def find_oil_stable():
    from VLM.bzFRInGE.plotting.plotTS_1D import find_stable
    Exp0  = ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=10)

    df_facs = pd.read_csv(op.join(Exp0.path_wd, 'Houston_facs.txt'),
                sep='\\s+', header=None, index_col=0, names='lat lon'.split())

    lat, lon = df_facs.iloc[0].to_numpy()
    loc = df_facs.iloc[0]

    lst_gdfs = []
    for sta, row in df_facs.iterrows():
        gdf_m = find_stable(Exp0, [lat, lon], rad=10000, loc=sta.replace('R6-TX-', ''), plot=False)
        lst_gdfs.append(gdf_m)

    gdf_ms = pd.concat(lst_gdfs).drop(columns='crs index_right'.split())
    gdf_ms.index.rename('name', inplace=True)

    dst = op.join(Exp0.path_wd, 'stable_oil.GeoJSON')
    gdf_ms.to_file(dst)
    print (f'Wrote: {dst}')


if __name__ == '__main__':
    Exp0 = ExpBase(Houston_SR, 'ERA5_SET_PM_ex_Fast', 'NASA', neofs=5)
    PObj = PlotHouston(Exp0, kind='rate')
    # PObj.plot_basic()
    # PObj.plot_inset('BBayou')
    # PObj.plot_inset('Nonlinear', double_diff=False)
    PObj.plot_insets('BBayou Texas_City'.split(), npix=2)
    # PObj.plot_transect_map()
    # plot_gps_insar(Exp0)
    # scatter_facilities_vlm(Exp0, npix=1)
    # plot_kde_facilities_vlm(Exp0, npix=1)

    # plot_insar_oil_box()

    # plot_oil_ts()

    # plot_oil_ts_dd()

    # plot_oil_ts_dd_all()

    path_figs = op.join(PATH_RES, f'{Exp0.reg}_2024')
    # bbPlot.savefigs(path_figs, True, True, dpi=300)

    plt.show()

