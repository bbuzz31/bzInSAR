"""
Plot NYC Maps in Python for Paper
"""
import string
import xarray as xr
from matplotlib.colors import TwoSlopeNorm, BoundaryNorm, Normalize
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
from cartopy.io import img_tiles as cimgt

from mintpy.utils import readfile
from BZ import bbLogger, bbGIS, bbPlot
from BZ.bbGIS import city_shape

from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from VLM.bzFRInGE.analysis import analyzeGPS

## for inset plots
from VLM.bzFRInGE.plotting.plotExp import *
from VLM.bzFRInGE.plotting.plot_TS1 import *


GS  = 13 # map tick size
FS1 = 9 # GNSS station name size


class PlotNYC(PlotExp):
    def __init__(self, exp, kind='rate', show_gps=True, continuous=True):
        super().__init__(exp, kind, show_gps, continuous)
        # self.use_stitched = False
        self.path_figs = op.join(PATH_RES, f'{exp.reg}_2023')
        self._set_plot_parms()


    def _set_plot_parms(self):
        super()._set_plot_parms()
        if self.kind == 'Rate':
            self.show_gps_names = True

        elif self.kind == 'Unc':
            pass

        else:
            raise Exception (f'Incorrect kind: {self.kind}. Use rate or unc')

        return


    def plot_basic_box_NYC(self, pts='LaGuardia Williams'.split(), letter=False):
        """ Plot basic plus a 2000 m bbox around point to show inset locations """
        fig, axes  = self.plot_basic()

        for i, pt in enumerate(pts):
            try:
                lalo = DCT_PTS[pt][0] #  use plot_ts to make a box around the point
            except:
                continue
            SNWE = buffer_point(*lalo, 2000)

            sty  = dict(zorder=50, alpha=0.9, facecolor='none', edgecolor='white',
                        linestyle='--', linewidth=1.5)
            bbox = bbGIS.bbox2poly(SNWE=SNWE)

            ## add the point in the middle
            pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                            color='w', transform=ccrs.PlateCarree(), s=30, marker='^')

            ## add triangle
            axes.scatter(lalo[1], lalo[0], **pt_parms)

            ## add letter or the box around the pt
            if letter:
                pt_parms.pop('s'); pt_parms.pop('marker')
                lbl = string.ascii_uppercase[i]
                dct_let = {'B':( -0.01, 0), 'C': (-0.007, -0.005), 'D': (-0.009, 0)}
                xoff, yoff = dct_let.get(lbl, (0, 0))
                axes.text(lalo[1]+xoff, lalo[0]+yoff, lbl, fontsize=18, **pt_parms)
            else:
                axes.add_geometries([bbox], crs=self.proj, **sty)

        ## add city borders and names
        # THESE ONLY LOOK WRITE WHEN SAVED WITH 300 DPI AND OPENED IN PREVIEW
        dct_names = {'Bronx': (-73.8145, 40.8060), 'Brooklyn': (-73.9762, 40.6355),
                    # 'Bronx': (-73.8145, 40.7950), 'Brooklyn': (-73.9762, 40.625)
                    # 'Bergen County': (-74.13, 40.79), 'Hudson': (-74.069, 40.650),
                    'Bergen County': (-74.13, 40.805), 'Hudson': (-74.069, 40.6627),
                    'Manhattan': (-73.990, 40.784), 'Queens': (-73.8634, 40.709),
                    # 'Manhattan': (-73.990, 40.770), 'Queens': (-73.8634, 40.695),
                    # 'Staten Island': (-74.111, 40.5155), 'Essex': (-74.19036, 40.780)}
                    'Staten Island': (-74.111, 40.5292), 'Essex': (-74.1902, 40.793)}

        for city in DCT_CITIES[self.reg]:
            continue
            loc = 4 if city == 'Essex' else 0
            city_wgs, city_utm = city_shape(city, loc)
            sty  = dict(zorder=50, alpha=0.9, facecolor='none', edgecolor='white',
                        linestyle=':', linewidth=1.0)
            axes.add_geometries(city_wgs.geometry, crs=self.proj, **sty)
            # if city == 'Bronx':
            #     breakpoint()

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
        pt_parms = dict(edgecolors='white', facecolor='none',
                        transform=ccrs.PlateCarree(), s=70, marker='v')
        axes.scatter(lalo[1], lalo[0], **pt_parms)


        if pts[0] is not None:
            fig.set_label(f'{fig.get_label()}_{"_".join(pts)}')
        else:
            fig.set_label(f'{fig.get_label()}_clean')

        return


    def plot_transects(self):
        transects = 'Transect_NJ1 Transect_Manhattan1 Transect_Brooklyn1'.split()
        fig, axes = super().plot_transects(transects)

        ## add the Battery TG
        lalo = (40.7010692,  -74.0143231)
        pt_parms = dict(edgecolors='white', facecolor='none',
                        transform=ccrs.PlateCarree(), s=70, marker='v')
        axes.scatter(lalo[1], lalo[0], **pt_parms)

        # axes = bbPlot._scale_bar(axes, 10, (.680, 0.9625))
        axes = bbPlot._scale_bar(axes, 10, (0.85, 0.015))
        return


    def plot_insets_col(self, pts='LaGuardia Williams'.split(), npix=3):
        """ Plot a close up of two locations, column wise """
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        self.show_gps = False
        fig  = plt.figure(figsize=(13, 13))
        gs   = GridSpec(nrows=2, ncols=2, figure=fig, wspace=0.25, hspace=0.01, height_ratios=[2.5,1])

        ## plot the map with s6, j3, tg
        axeMap0 = plt.subplot(gs[0, 0], projection=self.basemap.crs)
        axeMap1 = plt.subplot(gs[0, 1], projection=self.basemap.crs)
        axeMaps = [axeMap0, axeMap1]

        axeTS0    = fig.add_subplot(gs[1, 0])
        axeTS1    = fig.add_subplot(gs[1, 1])
        axeTSs    = [axeTS0, axeTS1]

        TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs, npix)
        for i, pt in enumerate(pts):
            lalo = DCT_PTS[pt][0] #  use plot_ts to make a box around the point
            S, N, W, E = buffer_point(*lalo, 2000)

            self.zoom = 15
            axes  = self.plot_basic(axes=axeMaps[i])[1]
            im    = axes.get_children()[0]
            im.colorbar.remove()
            im.set_alpha(0.75)
            axes.set_extent([W, E, S, N])

            gl = axes.gridlines(draw_labels=True)
            gparms = dict(left=False, right=True, size=GS-1) if i == 1 else {}
            bbPlot.fmt_gridlines(gl, bottom=True, **gparms)

            # axeTS = plot_ts_vup1_dd(self, pts[i], f'{pts[i]}_Stable', npix=npix, axes=axeTSs[i])[1]
            axeTS = TSObj.plot_ts_vup_dd(pts[i], axes=axeTSs[i])[1]
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


    def plot_inset(self, pt='LaGuardia', npix=3, show_stable=False):
        """ Plot one point, row wise """
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig = plt.figure(figsize=(18, 10), constrained_layout=True)
        gs  = GridSpec(2, 3, width_ratios=[1.3, 1, 1], figure=fig)
        axes = []
        axes.append(plt.subplot(gs[0, 0], projection=self.basemap.crs))
        axes.append(plt.subplot(gs[0, 1:]))

        lalo = DCT_PTS[pt][0] #  use plot_ts to make a box around the point
        S, N, W, E = buffer_point(*lalo, 2000)

        self.zoom  = 15

        ax  = self.plot_basic(axes=axes[0])[1]
        im  = ax.get_children()[0]
        im.colorbar.remove()
        im.set_alpha(0.75)
        gl = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 190, 0.01),
                    ylocs=np.arange(-90, 100, 0.01))

        gparms = dict(left=False, right=True)# if i == 1 else {}
        bbPlot.fmt_gridlines(gl, bottom=True, size=GS-1)#, **gparms)

        ## add actual point
        pt_parms = dict(path_effects=[pe.withStroke(linewidth=2, foreground='k')],
                        color='w', transform=ccrs.PlateCarree(), s=50)
        ax.scatter(lalo[1], lalo[0], **pt_parms)
        if show_stable:
            lalo1 = DCT_STABLE[self.mp_exp0][f'{pt}_Stable'][0]
            ax.scatter(lalo1[1], lalo1[0], marker='s', **pt_parms)

        else:
            ax.set_extent([W, E, S, N], crs=self.proj)

        TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs, npix)
        # axeTS = plot_ts_vup1_dd(self, pt, f'{pt}_Stable', npix=npix, axes=axes[1])[1]

        if pt == 'Woodside':
            # axeTS = plot_woodside_dd(self, npix=npix, axes=axes[1])[1]
            axeTS = TSObj.plot_ts_vup_piece2('20191231', pt, npix=npix, axes=axes[1])[1]

        else:
            axeTS = TSObj.plot_ts_vup_dd(pt, axes=axes[1])[1]

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
        # i accidentally forgot to do 3 for woodside but diff is <0.1 mm
        dct_npix = {'LaGuardia': 3, 'Williams': 3, 'Ashe': 3, 'Woodside': 0}
        for i, pt in enumerate(pts):
            lalo = DCT_PTS[pt][0] #  use plot_ts to make a box around the point
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

            TSObj = PlotTS1(self.dct_exp, self.mp_exp0, self.ref_sta, self.neofs, npix)

            if pt == 'Woodside':
                axeTS = TSObj.plot_ts_vup_piece2('20191231', pt, npix=npix, axes=axeTSs[1])[1]

            else:
                axeTS = TSObj.plot_ts_vup_dd(pts[i], axes=axeTSs[i])[1]

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
                   linestyle='--', label=f'{mu:.1f} $\\pm$ {lh:.1f}')

        ax.axvspan(mu-lh, mu+lh, ymax=0.65, color='darkred', alpha=0.2)
        letter     = string.ascii_uppercase[i]
        ax.text(0.0, 1.0, letter, transform=ax.transAxes,
                fontsize=15, verticalalignment='top', color='white',
                bbox={'facecolor':'k'})


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
                   linestyle='--', label=f'{mu:.1f} $\\pm$ {lh:.1f}')

        ax.axvspan(0, mu+lh, ymax=0.65, color='darkred', alpha=0.2)

        # ax.axvline(mu-sig, ymax=0.65, color='darkred', linestyle=':')
        # ax.axvline(mu+sig, ymax=0.65, color='darkred', linestyle=':')
        ax.grid(color='gray', linestyle = '--', linewidth=0.1)
        ax.set_xlabel(lbl)

        ax.legend()

    fig.set_label('Corr_Hists')


## ----------- more wrappers
def make_scratch(Exp1):
    from VLM.bzFRInGE.analysis import analyzeGPS
    PObj = PlotNYC(Exp1, kind='rate', show_gps=True)
    # PObj.plot_basic_box_NYC('Williams'.split())
    # PObj.plot_basic_box_NYC('Woodside'.split())
    # PObj.plot_insets('LaGuardia Williams'.split())
    PObj.plot_insets('Ashe Woodside'.split())

    # PObj.plot_basic_box_NYC(pts=[None])
    # PObj.plot_city()
    # PObj.plot_transects()
    # PObj.plot_basic_box_NYC(pts='LaGuardia Ashe Williams Woodside'.split(), letter=True)


    # GPSObj = analyzeGPS.GPS_Rate_Analysis(Exp1)
    # GPSObj.plot_trends_new()
    # GPSObj.plot_gps_vs_insar()
    # GPSObj.plot_gps_vs_insar(unc_lsq=True)


def make_main(Exp1):
    PObj = PlotNYC(Exp1, kind='rate', show_gps=True)
    # da_tay = get_tay22()
    # PObj.plot_basic_NYC(da=None)

    # PObj.plot_insets(pts='LaGuardia Williams'.split())
    # PObj.plot_insets(pts='Ashe Woodside'.split())
    # PObj.plot_city('Essex', 4)


def make_supp(Exp1):
    from VLM.bzFRInGE.analysis import analyzeGPS
    PObjU= PlotNYC(Exp1, kind='unc')
    PObjU.plot_basic()
    # PObjU.plot_gps()
    # PObjU.plot_basic_box_NYC(pts='Woodside'.split())
    # PObjR.plot_rm_background('Caron')
    return

    PObjR = PlotNYC(Exp1, kind='rate')
    # PObjR.plot_insets(pts='Ashe Woodside'.split())
    # PObjR.plot_transect()

    GPSObj = analyzeGPS.GPS_Rate_Analysis(Exp1)
    GPSObj.plot_gps_vs_insar()
    # GPSObj.plot_trends_new()

    # plot_hist_ref_sta()
    # plot_hist_corr()

    # PObjU.plot_polys()

    PObjR1= PlotNYC(Exp1, kind='rate')
    PObjR1.plot_rm_background('Garner')


def diagnostics(Exp1):
    from VLM.bzFRInGE.analysis import analyzeGPS
    PObjR = PlotNYC(Exp1, kind='rate')
    PObjU = PlotNYC(Exp1, kind='unc')
    PObjR.plot_basic()
    # PObjR.plot_hist()

    n_pix  = PObjR.da.where(np.isnan(PObjR.da), 1).sum()

    # get a count of where pixels are >= sigma
    da_sig  = xr.where(np.abs(PObjR.da)>=PObjU.da, 1, 0)
    per_sig = 100*(da_sig.sum() / n_pix)
    log.critical (f'Percent significant at 1sig: {per_sig:.1f}%')

    GPSObj = analyzeGPS.GPS_Rate_Analysis(Exp1)
    GPSObj.plot_gps_vs_insar()
    # GPSObj.plot_trends_new()


## need to refactor into class above
def plot_woodside_dd(Exp,  npix=0, mask=True, axes=None, split='20200101'):
    """ Plot the vertical displacement timeseries average (npix)

    at woodside with two trends
    """
    from mintpy.utils import readfile, utils as ut
    loc_target = 'Woodside'
    loc_stable = 'Woodside_Stable'

    vup_true  = readfile.read(Exp.path_vup_geo)[0]*1000
    unc, meta = readfile.read(Exp.path_vup_geo, datasetName='velocityStd')
    arr_iang  = readfile.read(Exp.path_geom_mp_geo, datasetName='incidenceAngle')[0]
    unc *= 1000
    df_gps    = prep_gps(Exp.path_gps, Exp.reg, units='mm').set_index('sta')
    ref_vel   = df_gps.loc[Exp.ref_sta, 'u_vel']
    log.debug('This may be the incorrect ref station depending on loc of ts point')

    # path_ts_geo = Exp.path_ts_geo if 'recon' in Exp.path_ts_geo else Exp.path_ts_geo_ann
    # log.info('Not using seasonally corrected')
    path_ts_geo = Exp.path_ts_geo if 'recon' in Exp.path_ts_geo else Exp.path_ts_geo
    if 'stitch' in Exp.mp_exp.lower():
        log.warning('Using stitched, recon time series')
        path_ts_geo = Exp.path_ts_geo_stitch


    with h5py.File(path_ts_geo, 'r') as h5:
        arr_ts = h5['timeseries'][:]
        arr_dt = [dt.decode('utf-8') for dt in h5['date'][:]]

    dates = pd.to_datetime(arr_dt)
    decyr0 = bbTS.date2dec(dates)

    lalo_target = DCT_PTS[loc_target][0]
    lalo_stable = DCT_STABLE[Exp.mp_exp0][loc_stable][0]

    ## convert whole timeseries to Vup
    arr   = 1000 * arr_ts / np.cos(np.deg2rad(arr_iang))

    if 'PM' in Exp.mp_exp.upper():
        arr_pm, meta1 = readfile.read(op.join(Exp.path_mp_exp_geo, 'ITRF14.h5'))
        assert 'REF_X' in meta1.keys(), 'Need to reference the plate motion'
        arr_pm0= (1000 * arr_pm / np.cos(np.deg2rad(arr_iang))).reshape(-1)
        inters = np.zeros_like(arr_pm0)
        # arr_pm = np.polyval([arr_pm.reshape(-1), np.zeros_like(arr_pm.reshape(-1))], decyr0)
        arr_pm1 = (arr_pm0[:, np.newaxis] * decyr0) + inters[:, np.newaxis]
        arr_pm  = arr_pm1.transpose(1, 0).reshape(arr_ts.shape)
        arr_pm -= arr_pm[0] # reference to first date
        arr    -= arr_pm

    if mask:
        mask = readfile.read(Exp.path_mask_vup)[0]
        arr *= np.where(np.isclose(mask, 0), np.nan, 1)

    ## average a bit around the target point
    coord    = ut.coordinate(meta, lookup_file=Exp.path_geom_mp_geo)
    y, x     = coord.geo2radar(lalo_target[0], lalo_target[1])[0:2]
    log.info ('Loc y/x: %s/%s', y, x)
    if npix == 0:
        ts_vup_target   = arr[:, y, x] # no real diff if surrounding points similar
        vup_true_target = vup_true[y, x]
        unc_target      = unc[y, x]
    else:
        ts_vup_target   = np.nanmean(arr[:, y-npix:y+npix, x-npix:x+npix], axis=(1, 2))
        vup_true_target = np.nanmean(vup_true[y-npix:y+npix, x-npix:x+npix])
        unc_target      = np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])

    rate_target, bias_target = np.polyfit(decyr0, ts_vup_target, 1)

    # Check rates calculate by hand match those in velocity file
    check_rate(vup_true_target, rate_target+ref_vel, npix)

    ## get coords of (stable) point
    y, x     = coord.geo2radar(lalo_stable[0], lalo_stable[1])[0:2]

    # use the stable point exactly
    ts_vup_stable   = arr[:, y, x] # no real diff if surrounding points similar
    vup_true_stable = vup_true[y, x]
    unc_stable      = unc[y, x]

    rate_stable, bias_stable = np.polyfit(decyr0, ts_vup_stable, 1)
    check_rate(vup_true_stable, rate_stable+ref_vel, npix)

    ## double difference
    ts_dd    = ts_vup_target - ts_vup_stable

    ## get gps for propagating the uncertainty
    if 'stitch' in Exp.mp_exp.lower():
        ## get the nearest gps to the target point
        import geopandas as gpd
        try:
            ref_sta = df_gps[df_gps.index.isin(Exp.ref_sta)]
        except:
            ref_sta = df_gps[df_gps.index == Exp.ref_sta]
        gdf_target = gpd.GeoDataFrame(crs=4326,
            geometry=gpd.points_from_xy([lalo_target[1]], [lalo_target[0]]))
        gdf = bbGIS.df2gdf(ref_sta, 'all').to_crs(4087)
        gdf_near = gdf_target.to_crs(4087).sjoin_nearest(gdf,
                                        distance_col='dist').to_crs(4326)
        gps_vel = ref_sta[ref_sta.index.isin(gdf_near.index_right)].u_vel.item()
        gps_sig = ref_sta[ref_sta.index.isin(gdf_near.index_right)].u_sig.item()

    ## mean uncertainty of the two stat
    elif not isinstance(Exp.ref_sta, str) and len(Exp.ref_sta) > 1:
        ref_sta = df_gps[df_gps.index.isin(Exp.ref_sta)]
        gps_vel = np.mean(ref_sta.u_vel)
        gps_sig = np.sqrt(np.mean(ref_sta.u_sig ** 2))

    elif not isinstance(Exp.ref_sta, str) and len(Exp.ref_sta) == 1:
        ref_sta = df_gps[df_gps.index.isin(Exp.ref_sta)]
        gps_vel = ref_sta.u_vel.item()
        gps_sig = ref_sta.u_sig.item()

    else:
        ref_sta = df_gps[df_gps.index == Exp.ref_sta]
        gps_vel = ref_sta.u_vel.item()
        gps_sig = ref_sta.u_sig.item()

    #  note pm may have been removed, and this is relative to GPS
    log.info (f'Vup Rate Target: {rate_target+gps_vel:.2f}, '\
              f'Vup Rate Stable: {rate_stable+gps_vel:.2f}')

    # for Emi, maybe something else?
    df  = pd.DataFrame({f'{loc_target}_dd': ts_dd, f'{loc_target}_stable': ts_vup_stable,
                       f'{loc_target}_target': ts_vup_target, 'decyr': decyr0, 'date':dates})
    corr  = True if not 'Base' in Exp.mp_exp0 else False
    corrl = '_corr' if corr else '_raw'
    df['corr'] = corr
    dst = f'{PATH_RES}/{loc_target}_dd{corrl}.csv'
    df.to_csv(dst)
    log.info(f'Wrote: %s', dst)

    rate_dd, unc_dd = bbTS.bootstrap_1D(ts_dd, decyr0, n_boots=1500)[:2]
    rss = np.polyfit(decyr0, ts_dd, 1, full=True)[1]
    rmse = np.sqrt(rss/len(ts_dd)).item()

    # sem2       = np.sqrt(cov_vup[0,0]+gps_sig**2)
    sem2       = np.sqrt(unc_dd**2+gps_sig**2)
    coeffs_vup, cov_vup = np.polyfit(decyr0, ts_dd, 1, cov='unscaled')
    trend_vup  = np.polyval(coeffs_vup, decyr0)


    ## plot the Vup
    if axes is None:
        fig, axes = plt.subplots(figsize=(10, 4))
        axes.set_title(f'{loc_target}-{loc_stable}')#-GPS')
    else:
        fig = axes.get_figure()

    ts0, dates0 = ts_dd, dates
    ix1, ix2 = dates0<=split, dates0>split
    # ix1 *= dates0>='20170101'
    colors = 'crimson gray'.split()
    for i, ix in enumerate([ix1, ix2]):
        dates, ts = dates0[ix], ts0[ix]
        decyr     = bbTS.date2dec(dates)
        unc_dd    = bbTS.bootstrap_1D(ts, decyr, n_boots=1500)[1]
        # hack to add the GPS uncertainty
        sem2      = np.sqrt(unc_dd**2 + 0.76**2)
        coeffs    = np.polyfit(decyr, ts, 1)
        rate      = coeffs[0]
        if np.abs(rate) < 1e-2:
            rate = np.abs(rate)
        trend     = np.polyval(coeffs, decyr)
        axes.scatter(dates, ts, c='k', s=15)
        l, = axes.plot(dates, trend, color=colors[i], linestyle='--',
                       label=f'{rate:.1f}$\\pm${unc_dd:.1f} (mm/yr)')
        axes.fill_between(dates, trend-sem2, trend+sem2, color=l.get_color(), alpha=0.3)

        # l1, = axes.fill(np.nan, np.nan, l.get_color(), alpha=0.3)

        axes.legend(prop ={"size": 13}, loc='lower right')
        # axes.legend([(l1, l)], )

    # axes.scatter(dates, ts_dd, c='k', s=15)
    # axes.plot(dates, trend_vup, color='r', linestyle='--', label=f'{rate_dd:.1f}$\\pm${sem2:.1f} (mm/yr)')
    # axes.fill_between(dates, trend_vup-sem2, trend_vup+sem2, color='r', alpha=0.3)
    # axes.set_ylabel('Vertical Displacement (mm)', fontsize=13)
    # axes.legend(prop ={"size": 12})
    # axes.tick_params(axis='x', labelsize=13)
    # fig.set_label(f'{loc_target}_dd')

    log.info(f'{loc_target} (npix={npix})')

    axes.set_title(f'Timeseries Std; RMSE; npix: {np.nanstd(ts_dd):.1f}; {rmse:.1f}; {npix}')
    axes.set_ylabel('Vertical Displacement (mm)', fontsize=13)
    axes.grid(color='gray', linestyle = '--', linewidth=0.1)

    fig.set_label(f'tsdd_{loc_target}_{npix}_split')
    # log.info ('Double Diff Rate, Timeseries Std, RMSE: %.1f, %.1f, %.1f',
    #             coeffs_vup[0], np.nanstd(ts_dd), rmse)

    ## show correlation
    rho0 = np.corrcoef(ts_vup_target, ts_vup_stable)[0,1]
    ts_vup_target_detrend = ts_vup_target - np.polyval([rate_target, bias_target], decyr0)
    ts_vup_stable_detrend = ts_vup_stable - np.polyval([rate_stable, bias_stable], decyr0)
    rho_detrend = np.corrcoef(ts_vup_target_detrend, ts_vup_stable_detrend)[0,1]
    # fig, axes1 = plt.subplots(figsize=(10, 6), nrows=2, sharex=True)
    # axes1[0].scatter(dates0, ts_vup_target, s=15, color='k')
    # axes1[0].scatter(dates0, ts_vup_target_detrend, s=15, color='r')
    # axes1[0].set_title(loc_target)
    # axes1[1].scatter(dates0, ts_vup_stable, s=15, color='k')
    # axes1[1].scatter(dates0, ts_vup_stable_detrend, s=15, color='r')
    # axes1[1].set_title(loc_stable)
    # for ax in axes1:
    #     ax.set_ylabel('Vertical Displacement (mm)')
    #     ax.grid(color='k', alpha=0.1, linestyle='--')

    log.info(f'Correlation between {loc_target} and {loc_stable}: {rho0:.2f}')
    log.info (f'Correlation between DETRENDED {loc_target} and {loc_stable}: {rho_detrend:.2f}')

    return fig, axes


if __name__ == '__main__':
    Exp0      = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20)
    # Exp0      = ExpBase(NYC_ARIA, 'ERA5_SET_PM_ex_Fast', 'NYBK', neofs=0)

    make_scratch(Exp0)
    # plot_together(Exp0)
    # make_main(Exp0)

    # make_supp(Exp0)
    # make_numbers(Exp0, PlotNYC)
    # diagnostics(Exp0)

    # GPSObj = analyzeGPS.GPS_Rate_Analysis(Exp0)
    # GPSObj.plot_gps_vs_insar()

    path_figs = op.join(PATH_RES, f'{Exp0.reg}_2023')
    # bbPlot.savefigs(path_figs, True, True, dpi=300)

    # import subprocess
    # subprocess.call('open /Users/buzzanga/Desktop/VLM/NYC_2023/Vup_NYC_Rate_LaGuardia_Williams.png'.split())


    plt.show()