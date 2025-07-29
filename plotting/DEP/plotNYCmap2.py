"""
Plot NYC Maps in Python for Paper
"""
import xarray as xr
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
from cartopy.io import img_tiles as cimgt

from BZ import bbLogger, bbGIS, bbPlot

from VLM.bzFRInGE.contrib import *
from VLM.bzFRInGE.contrib.experiments import *
from VLM.bzFRInGE.contrib.FRINGEBase import ExpBase


## for inset plots
from VLM.bzFRInGE.contrib.plotting.plotTS import *

DCT_TRANSECTS = {'NYC': 'NJ Manhattan Brooklyn'.split()}

## plot parms
class PlotNYC(ExpBase):
    def __init__(self, exself, kind='rate', show_gps=True, use_stitched=False):
        super().__init__(exself.dct_exp, exself.mp_exp0, exself.ref_sta, exself.neofs)
        self.show_gps     = show_gps
        self.use_stitched = use_stitched
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
            self.pparms = {'cmap': 'cmc.vik', 'norm': mpl.colors.TwoSlopeNorm(0, -5, 5)}
            # self.cbar_lbl = 'V$_{Up}$ (mm/yr)'
            self.cbar_lbl = 'Vertical Rate (mm/yr)'
        elif self.kind == 'Unc':
            self.pparms = {'cmap': 'cmc.lajolla', 'norm':mpl.colors.Normalize(0, 1.5)}
            self.cbar_lbl = 'V$_{Up}$ Uncertainty (mm/yr)'

        else:
            raise Exception ('Incorrect kind: {self.kind}. Use rate or unc')

        return


    def get_da(self):
        """ Get the masked netcdf of the data to plot """
        if self.use_stitched:
            path = self.path_rate_msk_nc_stitch if self.kind == 'Rate' else self.path_std_msk_nc_stitch
        else:
            path = self.path_rate_msk_nc if self.kind == 'Rate' else self.path_std_msk_nc

        # get the netcdf
        ds = xr.open_dataset(path)
        da = ds['Band1'].rename(self.kind)
        da = da*1000 if 'm' in da.units else da
        return da


    def plot_basic(self, axes=None):
        """ Plot rate or uncertainty in a single panel w/wout GPS """
            # path = self.path_rate_nc if self.kind == 'Rate' else self.path_std_nc
        da = self.da.copy()

        if axes is None:
            fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
            axes.set_extent(self.WESN, crs=self.proj)
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, bottom=True)
        else:
            fig = axes.get_figure()

        axes.add_image(self.basemap, self.zoom)

        ## add the interpolated GPS rates
        im   = axes.pcolormesh(da.lon, da.lat, da, shading='nearest', alpha=self.alpha,
                        transform=self.proj, **self.pparms)

        bbPlot.cartopy_cbar(im, ylabel=self.cbar_lbl)
        fig.set_label(f'Vup_{self.reg}_{self.kind}')

        ## add the GPS
        if self.show_gps:
            self.plot_gps(axes)


        # print(f'Vup {self.kind}: {da.mean().item():.2f} +/- {da.std().item():.2f} mm/yr')
        l, h, mu = da.min().item(), da.max().item(), da.mean().item()
        print (f'InSAR Vup {self.kind}: {mu:.2f} [{l:.2f}, {h:.2f}]')

        interactive = False
        if interactive:
            import hvplot.xarray
            ## otherwise super slow (LaGuardia)
            da1 = da.sel(lon=slice(-73.88460830568238, -73.8486756943176),
                                lat=slice(40.75743769431761, 40.79337030568239))
            plot = da1.hvplot.quadmesh(*f'lon lat {self.kind}'.split(),
                    coastline='10m', geo=True,  cmap='coolwarm', clim=(-5, 5))
            hvplot.show(plot)


        return fig, axes


    def plot_basic_box(self, pts='LaGuardia Williams'.split()):
        """ Plot basic plus a 2000 m bbox around point to show inset locations """
        fig, axes  = self.plot_basic()

        for pt in pts:
            lalo = dct_pts[pt][0] #  use plot_ts to make a box around the point
            SNWE = buffer_point(*lalo, 2000)

            sty  = dict(zorder=50, alpha=0.9, facecolor='none', edgecolor='white',
                        linestyle='--', linewidth=1)
            bbox = bbGIS.bbox2poly(SNWE=SNWE)
            axes.add_geometries([bbox], crs=self.proj, **sty)
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

        else:
            fig = axes.get_figure()
            s   = axes.scatter(self.df_gps['lon'], self.df_gps['lat'], **msty)

        # Add a square reference pixel
        ref_sta = self.df_gps[self.df_gps.sta == self.ref_sta]
        msty    = {**msty, 'marker':'s', 'c':ref_sta['u_vel'], 's': 75}
        s = axes.scatter(ref_sta['lon'], ref_sta['lat'],  **msty)

        # Add text x pixels of the station for rates only
        if self.kind == 'Rate':
            # sty_name = dict(facecolor='sandybrown', alpha=0.5, boxstyle='round')
            sty_name = dict(facecolor='lightgray', alpha=0.5, boxstyle='round')
            names = self.df_gps['sta'].astype('str').values
            for i in range(len(names)):
                geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=+20)
                axes.text(self.df_gps['lon'].iloc[i], self.df_gps['lat'].iloc[i],
                        names[i].strip(), transform=text_transform, zorder=20,
                        verticalalignment='center', horizontalalignment='right',
                        bbox=sty_name)

        bbLogger.logger.info(
            f'MIDAS mu: {self.df_gps.u_vel.mean():.2f} +/- {self.df_gps.u_sig.mean():.2f} mm/yr')

        axes.set_extent(self.WESN, crs=self.proj)
        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)
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
                            color='k', transform=ccrs.PlateCarree(), fontsize=14)
            if i == 1:
                text_parms['va'] = 'bottom'
                en[1] -= 0.005
                st[0] += 0.004
            letter     = string.ascii_uppercase[i]
            axes.text(st[0], st[1], letter, **text_parms)
            axes.text(en[0], en[1], f"{letter}'", **text_parms)
        fig.set_label(f'{fig.get_label()}_transects')
        return


    def plot_hist(self):
        """ Histogram of rate/unc"""
        # shoudl use xarray, its faster
        bbPlot.plot_hist(self.da.data.reshape(-1))
        return


    def plot_insets_col(self, pts='LaGuardia Williams'.split()):
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
            gparms = dict(left=False, right=True) if i == 1 else {}
            bbPlot.fmt_gridlines(gl, bottom=True, **gparms)

            axeTS = plot_ts_vup1_dd(self, pts[i], f'{pts[i]}_Stable', npix=0, axes=axeTSs[i])[1]
            axeTS.set_title('')
            if i == 1:
                axeTS.set_ylabel('')
            else:
                axeTS.set_ylim([-40, 20])

            # axeTSs[i].plot(np.random.rand(180))
        # axeTS0.plot(np.random.rand(180))
        # axeTS1    = fig.add_subplot(gs[1, 1])#, sharey=axeTS0)
        # axeTS1.plot(np.random.rand(180))

        cbar_ax = fig.add_axes([0.5, 0.377, 0.02, 0.457]) # left, bot, wid, hgt
        cbar_ax.set_title('V$_{Up}$')
        cbar_ax.set_xlabel('(mm yr$^{-1}$)')
        cbar    = plt.colorbar(im, cax=cbar_ax,  spacing='uniform')


    def plot_inset(self, pt='LaGuardia'):
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
        bbPlot.fmt_gridlines(gl, bottom=True)#, **gparms)

        axeTS = plot_ts_vup1_dd(self, pt, f'{pt}_Stable', npix=0, axes=axes[1])[1]
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
        fig = plt.figure(figsize=(14, 9), constrained_layout=True)
        gs  = GridSpec(2, 3, width_ratios=[1, 1, 1], hspace=0.15, wspace=0.01, figure=fig)
        axeMaps = []
        axeTSs  = []
        axeMaps.append(plt.subplot(gs[0, 0], projection=self.basemap.crs))
        axeMaps.append(plt.subplot(gs[1, 0], projection=self.basemap.crs))
        axeTSs.append(plt.subplot(gs[0, 1:]))
        axeTSs.append(plt.subplot(gs[1, 1:]))
        npix = 1
        for i, pt in enumerate(pts):
            lalo = dct_pts[pt][0] #  use plot_ts to make a box around the point
            S, N, W, E = buffer_point(*lalo, 2000)

            self.zoom  = 15

            ax  = self.plot_basic(axes=axeMaps[i])[1]
            im  = ax.get_children()[0]
            im.set_alpha(0.75)
            im.colorbar.remove()
            ax.set_extent([W, E, S, N])

            pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                            color='w', transform=ccrs.PlateCarree(), s=30)
            ax.scatter(lalo[1], lalo[0], **pt_parms)

            gl = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 190, 0.01),
                        ylocs=np.arange(-90, 100, 0.01))

            bbPlot.fmt_gridlines(gl, left=True, right=False, bottom=True)

            axeTS = plot_ts_vup1_dd(self, pt, f'{pt}_Stable', npix=npix, axes=axeTSs[i])[1]
            # axeTS = axeTSs[i].plot(np.random.rand(180))
            axeTS.set_title('')
            axeTS.set_ylim([-40, 10]) if pt == 'LaGuardia' else ''

            axeTS.tick_params(labelbottom=False) if i == 0 else ''

            # ## make colorbars
            # divider = make_axes_locatable(axeMaps[i])
            # cbar_ax = divider.new_horizontal(size="5%", pad=0.55, axes_class=plt.Axes, pack_start=True)
            # fig.add_axes(cbar_ax)


            # cbar    = fig.colorbar(im, cbar_ax, label='Vertical Rate (mm/yr)')
            # cbar.ax.yaxis.set_label_position('left')
            # cbar_ax.yaxis.set_ticks_position('left')

        # make a single colorbar across
        cbar_ax = fig.add_axes([0.05, 0.490, 0.21, 0.0175]) # left, bot, wid, hgt
        cbar_ax.set_title('V$_{Up}$  (mm yr$^{-1}$)', fontsize=11)#, labelpad=1v#, rotation=270, )
        # cbar_ax.yaxis.set_label_position('right')
        cbar    = plt.colorbar(im, cax=cbar_ax, orientation='horizontal', spacing='uniform')
        # gs.update(left=0.07, right=0.23)#, bottom=0.1, top=0.9)
        fig.set_label(f'Insets_{"_".join(pts)}')
        return

    def plot_cities(self, city='NJ', show=True):
        """ Plot the rates with city boundaries and calculate the mean

        Manhattan, Brooklyn, Queens, Staten Island, NJ
        """
        import rioxarray as xrr
        from shapely.geometry import mapping
        from BZ.bbGIS import city_shape
        city_wgs, city_utm = city_shape(city)

        da = self.da.copy()
        da.rio.set_spatial_dims(x_dim='lon', y_dim='lat', inplace=True)
        da.rio.write_crs(4326, inplace=True)

        da_reg = da.rio.clip(city_wgs.geometry.apply(mapping))

        if show:
            fig, axes = plt.subplots(figsize=(10,10),
                                    subplot_kw={'projection': self.basemap.crs})
            axes.set_extent(self.WESN, crs=self.proj)
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, bottom=True)

            axes.add_image(self.basemap, self.zoom)

            ## add the interpolated insar rates
            im   = axes.pcolormesh(da_reg.lon, da_reg.lat, da_reg, **self.pparms)

            bbPlot.cartopy_cbar(im, ylabel=self.cbar_lbl)
            fig.set_label(f'Vup_{self.reg}_{self.kind}')

        l, h, mu = da_reg.min().item(), da_reg.max().item(), da_reg.mean().item()
        bbLogger.logger.info (f'{city} Vup {self.kind} mus: {mu:.1f} [{l:.1f}, {h:.1f}]')
        return



## ----------- Wrappers
def make_main(Exp1):
    PObj = PlotNYC(Exp1, kind='rate')
    PObj.plot_basic_box(pts='LaGuardia Williams'.split())
    PObj.plot_insets(pts='LaGuardia Williams'.split())


def make_supp(Exp1):
    PObj = PlotNYC(Exp1, kind='unc')
    PObj.plot_basic_box(pts='Ashe Woodside'.split())
    PObj.plot_insets(pts='Ashe Woodside'.split())
    PObj.plot_transect()


def make_numbers(Exp1):
    PObj = PlotNYC(Exp1, kind='rate')
    l, h, mu = PObj.da.min().item(), PObj.da.max().item(), PObj.da.mean().item()

    ## overall InSAR mean
    bbLogger.logger.info (f'InSAR Vup {PObj.kind}: {mu:.2f} [{l:.2f}, {h:.2f}]')

    ## overall GPS mean
    print ('Stations:', ', '.join(PObj.df_gps.sta.tolist()))
    bbLogger.logger.info(
        f'MIDAS mu: {PObj.df_gps.u_vel.mean():.2f} +/- {Oobj.df_gps.u_sig.mean():.2f} mm/yr'
    )

    ## cities
    for city in 'NJ Manhatten Brooklyn Queens'.split() + ['Staten Island']:
        PObj.plot_cities(city, show=False)




if __name__ == '__main__':
    Exp       = ExpBase(NYC_SRc, 'ERA5_SET_PM_Fast2017', 'NYBK', neofs=20)
    make_main(Exp)
    make_supp(Exp)

    path_figs = op.join(op.expanduser('~'), 'Desktop', 'VLM', 'NYC_2023')
    bbPlot.savefigs(path_figs, False, True)
    plt.show()