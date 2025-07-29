""" Plot all kinds of VLM together: UNR GPS, UNR IMAGING (grid/TG), GIA (Caron, Peltier)  """
from BZ import *
from BZ import bbGIS, bbPlot
from VLM.bzFRInGE import *
import matplotlib as mpl
# import cartopy.feature as cfeature
import pyproj
import cartopy.crs as ccrs
from cartopy.io import img_tiles as cimgt
import cmcrameri
import xarray as xr
import geopandas as gpd
from shapely.geometry import box

from VLM.bzFRInGE.FRINGEBase import ExpBase

# """
# This only works in VLM environment
# """
# assert os.environ['CONDA_DEFAULT_ENV'] == 'VLM', 'Must be in VLM environment'

ALPHA    = 0.85

## region to cmap, colorbounds, custom fuction
DCT_PARMS = {'HR': ['Blues_r', np.arange(-3, 0.5, 0.5)],
             'Charleston': ['Blues_r', np.arange(-2, 0.25, 0.25)],
             'Maine': ['cmc.vik', np.arange(-2, 2.5, 0.5)],
             'Savannah': ['cmc.vik', np.arange(-2, 2.5, 0.5)],
             'NYC': ['Blues_r', np.arange(-3, 0.25, 0.25)],
            #  'NYC': ['Blues_r', np.arange(-3, 0.5, 0.5)],
             'NNJ': ['Blues_r', np.arange(-3, 0.5, 0.5)],
             'Kennedy': ['cmc.vik', np.arange(-2, 2.5, 0.5)],
            #  'Houston': ['cmc.vik', np.arange(-5, 6, 1)] # GIA,
             'Houston': ['cmo.ice', np.arange(-10, 1, 1)] ,
             'DC': ['cmc.vik', np.arange(-2, 2.5, 0.5)],
             'Philly': ['cmc.vik', np.arange(-3, 3.5, 0.5)],
            }

# ToDo:
#     Support TG Imaging VLM


def load_gridded_imaging(var='VU'):
    """ Load/reformat the UNR Imaging data at http://geodesy.unr.edu/vlm.php

    VU or SVU (sigma vup, formal uncertainty)
    """
    assert var in 'VU SVU'.split(), 'Incorrect variable'
    Obj  = bzBase()
    src  = op.join(Obj.path_vlm, 'GNSS', 'UNR', 'imaging', 'VLM_Global_Imaged')
    dst  = f'{src}_BB.nc'
    if op.exists(dst):
        return xr.open_dataset(dst)[var]

    ds   = xr.open_dataset(f'{src}.nc')
    lat  = ds.LatGrid.data[0, :]
    lon  = ds.LonGrid.data[:, 0]
    lst_das = []
    for var in 'VU SVU'.split():
        dat  = ds[var].data
        da   = xr.DataArray(dat.T, name=var, coords={'lat': lat, 'lon': lon},
                    dims=['lat', 'lon'])
        lst_das.append(da)
    da_m = xr.merge(lst_das)
    da_m.to_netcdf(dst)
    return da


def load_GIA(kind='ICE6G'):
    cols = ['VLM']
    if kind == 'ICE6G':
        src = op.join(bzBase().path_vlm, 'GIA', f'drad.12mgrid_512.nc')
        assert op.exists(src), f'{src} doesnt exist'
        ds  = xr.open_dataset(src).rename(Lat='lat', Lon='lon', Drad_250='VLM')
        ds  = ds.assign_coords({'lon':
                                (((ds.lon + 180) % 360) - 180)}).sortby('lat lon'.split())

    else:
        src = op.join(bzBase().path_vlm, 'GIA', 'GIA_Caron_05.nc')
        ds  = xr.open_dataset(src).rename(x='lon', y='lat', rad_mean='VLM',
                                         rad_sterr='VLM_unc')
        ds  = ds.assign_coords({'lon':
                    (((ds.lon + 180) % 360) - 180)}).sortby('lat lon'.split())
        cols.append('VLM_unc')

    return ds[cols]


class PlotVLM(ExpBase):
    def __init__(self, dct_exp, mp_exp='Base_ex_Fast', ref_sta=None, neofs=15):
        super().__init__(dct_exp, mp_exp, ref_sta, neofs)
        ## in case you want a bigger area
        # self.SNWE = self.SNWE[0]-1, self.SNWE[1]+1, self.SNWE[2]-1, self.SNWE[3]+1
        self.path_res  = Path.home() / 'Desktop' / 'VLM' / f'{self.reg}_2024'
        self.lst_parms = DCT_PARMS[self.reg]
        # self.basemap = cimgt.GoogleTiles(url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
        self.basemap   = cimgt.GoogleTiles(url='https://server.arcgisonline.com/'\
                        'arcgis/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}.jpg')
        self.norm    = mpl.colors.BoundaryNorm(boundaries=self.lst_parms[1], ncolors=256)
        # self.norm    = mpl.colors.Normalize(self.lst_parms[1].min(), self.lst_parms[1].max())

        self.proj = ccrs.PlateCarree()

        self.WESN = [*self.SNWE[2:], *self.SNWE[:2]]
        self.bounds = box(self.SNWE[2], self.SNWE[0], self.SNWE[3], self.SNWE[1])


        self.df_gps  = prep_gps(self.path_gps, self.reg, units='mm')

        # self.df      = bbGIS.in_dfll_poly(self.df, op.join(self.path_crops, 'NYC_East.GeoJSON'))


    def __call__(self):
        # self.plot_gps_ts() # gps time series, play with jumps
        self.plot_mu()
        self.plot_GPS()
        self.plot_interp_GPS()
        self.plot_imaging_GPS_grid()
        self.plot_GIA('ICE6G')
        self.plot_GIA('Caron')


    def plot_mu(self):
        """ Plot the mean and std of each type of VLM if exists """
        fig, axes = plt.subplots(figsize=(10,4))
        lbls = []
        for i, k in enumerate('InSAR MIDAS ICE6G Caron Imaging'.split()):
            if k in 'ICE6G Caron Imaging'.split():
                if k == 'Imaging':
                    da     = load_gridded_imaging()
                else:
                    ds = load_GIA(k).load()
                    da = ds['VLM'] # VLM_unc

                ## just cropping the plot (for the means)
                # buff = 1e-3
                # da_reg = da.sel(lat=slice(self.SNWE[0]-buff, self.SNWE[1]+buff),
                #                 lon=slice(self.SNWE[2]-buff, self.SNWE[3]+buff))
                da = da.rio.write_crs(4326).rio.set_spatial_dims(x_dim='lon', y_dim='lat')
                da_reg = da.rio.clip([self.bounds], all_touched=True).drop('spatial_ref')

                # adjust lons
                if not da_reg.data.any():
                    log.warning(f'Taking nearest single {k} data point in {self.reg}')
                    # take the average of the coordinates and find the closest
                    lon = np.mean(self.SNWE[2:])
                    lat = np.mean(self.SNWE[:2])
                    da_reg = da.sel(lon=lon, lat=lat, method='nearest')
                    df = pd.Series(da_reg.item(), name=k).to_frame()
                else:
                    df = da_reg.rename(k).to_dataframe().reset_index().drop(
                        'lat lon'.split(), axis=1)

            elif k == 'MIDAS':
                df = self.df_gps['u_vel'].rename(k)

            elif k == 'InSAR':
                from mintpy.utils import readfile
                path_vup = self.path_vup_geo_stitch if 'Stitch' \
                    in self.mp_exp0.lower() else self.path_vup_geo
                if op.exists(path_vup):
                    arr_vel = readfile.read(path_vup)[0]*1000
                    df = pd.Series(arr_vel.ravel())
                else:
                    log.error(f'Couldnt get InSAR mean/std. {path_vup} doesnt exist.')
                    continue

            # dont compute std for very small
            if df.dropna().empty:
                log.info(f'{k} data is all nan!')

            elif df.shape[0] < 2 or np.isnan(df.std().item()):
                axes.scatter(str(i), df[k].dropna(), color='k')
            else:
                mu, std = df.mean(), df.std()
                axes.errorbar(str(i), mu, yerr=std, color='k', fmt='o')

            lbls.append(k)

        # axes.set_xticklabels(lbls)
        axes.set_xticks(axes.get_xticks(), lbls)
        axes.grid(color='gray', linestyle = '--', linewidth=0.1)
        axes.set_ylabel('Vup (mm/yr; $\\pm1\\sigma$ Variability)')
        fig.set_label('VLM_scatter')
        return


    def plot_GPS(self, axes=None, vel_col='u_vel', df=None):
        """ Plot the MIDAS GPS locations (+ rates and names if axes=None) """
        df   = self.df_gps if df is None else df
        lbl = 'Up' if vel_col == 'u_vel' else 'E' if vel_col == 'e_vel' else 'N'
        msty = dict(marker='o', s=50, linewidth=1.5, transform=self.proj,
                      cmap=self.lst_parms[0], norm=self.norm, alpha=1, zorder=30, edgecolors='white')
        # fig  = None
        if axes is None:
            fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
            axes.add_image(self.basemap, 10)
            s = axes.scatter(df['lon'], df['lat'], c=df[vel_col], **msty)

            bbPlot.cartopy_cbar(s, f'V$_{{{lbl}}}$ (mm/yr)')
            axes.set_title(f'MIDAS GPS {self.reg}')
            fig.set_label(f'Midas_{self.reg}')

            ## format the names
            sty_name = dict(facecolor='sandybrown', alpha=0.3, boxstyle='round')
            names = df['sta'].astype('str').values
            # Add text x pixels of the station.
            for i in range(len(names)):
                geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=+5+i)
                axes.text(df['lon'].iloc[i], df['lat'].iloc[i],
                        names[i].strip(), transform=text_transform, zorder=20,
                        verticalalignment='center', horizontalalignment='right',
                        bbox=sty_name)

            mu, sig = df.u_vel.mean(), df.u_sig.mean()
            mi, ma  = df.u_vel.min(), df.u_vel.max()

            logger.info(f'Raw GPS Min/Mean/Max: {mi:.2f} / {mu:.2f} +/- {sig:.2f} / {ma:.2f} (mm/yr)')
            axes.set_extent(self.WESN, crs=self.proj)
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, bottom=True)
        else:
            s = axes.scatter(df['lon'], df['lat'], c=df[vel_col], **msty)
            fig = s.axes.get_figure()

        return fig, axes


    def plot_gps_ts(self):
        """ Plot the GPS time series in the stuy reagion a linear fit, and the MIDAS rate

        Manually play with start/end to see how jumps look
        """
        from VLM.bzGNSS.get_MIDAS import GPS_TS
        from BZ import bbTS
        GPSObj = GPS_TS(self.reg)
        for ix, row in self.df_gps.iterrows():
            sta, uvel, usig = row['sta'], row['u_vel'], row['u_sig']
            if not sta == 'ROG1':
                continue
            fig, axes = GPSObj.plot_ts_trend(sta, st='20150101')# en='20200701')# , )
            x     = axes.get_lines()[0].get_xdata()
            xdec  = bbTS.date2dec(pd.to_datetime(x))
            trend = np.polyval([uvel, 0], xdec)
            axes.plot(x, trend, color='g', linestyle=':',
                      label=f'{uvel:.2f}$\\pm${usig:.2f} (MIDAS) ')
            axes.legend()
        return


    def plot_interp_GPS(self):
        """ Interpolate and plot the MIDAS GPS data """
        ds  = bbGIS.pts2grd(self.df_gps, 'u_vel', '4087', 90, 90, self.SNWE)
        arr = ds.ReadAsArray()

        WESN   = bbGIS.get_extent(ds)
        epsg   = pyproj.CRS(ds.GetProjection()).to_epsg()
        proj   = ccrs.epsg(epsg) if int(epsg) != 4326 else ccrs.PlateCarree()

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent(WESN, crs=proj)
        axes.add_image(self.basemap, 10)

        ## add the interpolated GPS rates
        im   = axes.imshow(arr, cmap=self.lst_parms[0], norm=self.norm,
                        zorder=20, interpolation='nearest', alpha=ALPHA,
                        extent=WESN, transform=proj, origin='upper')
        bbPlot.cartopy_cbar(im, 'V$_{up}$ (mm/yr)')

        axes.set_title(f'Interpolated MIDAS GPS {self.reg}')
        fig.set_label(f'GPSInterp_{self.reg}')

        ## add the GPS
        self.plot_GPS(axes)

        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)
        logger.info(f'Interpolated GPS: {np.nanmean(arr):.2f} +/- {np.nanstd(arr):.2f} mm/yr')
        return fig, axes


    def plot_imaging_GPS_grid(self):
        """ Plot the gridded UNR imagery for study area """
        da = load_gridded_imaging('VU')
        da = da.rio.write_crs(4326).rio.set_spatial_dims(x_dim='lon', y_dim='lat')
        da_reg = da.rio.clip([self.bounds], all_touched=True)

        ## to show the larger region
        # da_reg = da.sel(lat=slice(np.floor(self.SNWE[0]), np.ceil(self.SNWE[1])),
                        # lon=slice(np.floor(self.SNWE[2]), np.ceil(self.SNWE[3])))
        buff = 1e-2
        Sb, Nb, = self.SNWE[0]-buff, self.SNWE[1]+buff
        Wb, Eb, = self.SNWE[2]-buff, self.SNWE[3]+buff

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent(self.WESN, crs=self.proj)
        axes.set_extent([Wb, Eb, Sb, Nb], crs=self.proj)
        # axes.add_image(self.basemap, 10)

        norm = self.norm # mpl.colors.BoundaryNorm(np.arange(-2, 0, 0.1), 256)

        ## add the interpolated GPS rates
        im   = axes.pcolormesh(da_reg.lon, da_reg.lat, da_reg, shading='nearest', alpha=ALPHA,
                        transform=self.proj, cmap=self.lst_parms[0], norm=norm)
        bbPlot.cartopy_cbar(im, ylabel='V$_{up}$ (mm/yr)')
        axes.set_title(f'Imaging (Grid) {self.reg} V$_{{up}}$')
        fig.set_label(f'Imaging_grid_{self.reg}')


        ## add the GPS
        self.plot_GPS(axes)
        logger.info(f'Imaging: {da_reg.mean().item():.2f} +/- {da_reg.std().item():.2f} mm/yr')

        # bbPlot.plot_bbox(SNWE=self.SNWE, axes=axes, info=False, edgecolor='w')

        return fig, axes


    ## this (and imaging) should be a weighted average by the actual scene coverage
    def plot_GIA(self, model='ICE6G', interactive=False):
        """ Plot regional GIA from ICE6G or Caron """
        ds = load_GIA(model).load()
        da = ds['VLM'] # VLM_unc

        ## just cropping the plot (for the means)
        # da_reg = da.sel(lat=slice(self.SNWE[0], self.SNWE[1]),
                        # lon=slice(self.SNWE[2], self.SNWE[3]))
        da = da.rio.write_crs(4326).rio.set_spatial_dims(x_dim='lon', y_dim='lat')
        da_reg = da.rio.clip([self.bounds], all_touched=True)

        sig = da_reg.std().item()
        if not da_reg.data.any():
            # take the average of the coordinates and find the closest
            lon = np.mean(self.SNWE[2:])
            lat = np.mean(self.SNWE[:2])
            da_reg = da.sel(lon=lon, lat=lat, method='nearest')
            sig = np.nan


        # to see the values under cursor
        if interactive:
            import hvplot.xarray
            da_plot = da if da_reg.size == 1 else da_reg
            plot = da_plot.hvplot.quadmesh(*f'lon lat VLM'.split(),
                    coastline='10m', geo=True,  cmap=self.lst_parms[0],
                    width=500)#, clims=(self.norm.vmin, self.norm.vmax))
            hvplot.show(plot)

        ## show a litle bit extra
        buff = 1e-2
        Sb, Nb, = self.SNWE[0]-buff, self.SNWE[1]+buff
        Wb, Eb, = self.SNWE[2]-buff, self.SNWE[3]+buff

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent([Wb, Eb, Sb, Nb], crs=self.proj)
        axes.add_image(self.basemap, 10)

        ## add the interpolated GPS rates
        im   = axes.pcolormesh(da.lon, da.lat, da, shading='nearest', alpha=ALPHA,
                        transform=self.proj, cmap=self.lst_parms[0], norm=self.norm)
        cb = bbPlot.cartopy_cbar(im, ylabel='V$_{up}$ (mm/yr)', pad=0.15)
        cb.set_ticks(self.norm.boundaries)

        axes.set_title(f'GIA: {model}')

        ## add the GPS
        self.plot_GPS(axes)

        bbPlot.plot_bbox(SNWE=self.SNWE, axes=axes, info=False, edgecolor='w')

        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)

        fig.set_label(f'GIA_{model}')
        # print (da_reg.data)
        logger.info(f'{model}: {da_reg.mean().item():.2f} +/- {sig:.2f} mm/yr')
        return fig, axes


    def plot_interp_GPS_horizontal(self):
        """ Interpolate and plot the MIDAS GPS data """
        df_gps = prep_gps(self.path_gps, self.reg, horiz=True, units='mm')

        msty = dict(marker='o', s=50, linewidth=1.5,
                alpha=1, zorder=30, edgecolors='white', cmap=self.lst_parms[0])

        fig, axes = plt.subplots(figsize=(10,10),
                                 subplot_kw={'projection': self.basemap.crs}, sharey=True, ncols=2)

        for i, ax in enumerate(axes):
            vcol = 'e_vel' if i == 0 else 'n_vel'
            lbl = 'E' if i == 0 else 'N'
            ds = bbGIS.pts2grd(df_gps, vcol, '4087', 90, 90, self.SNWE)
            epsg   = pyproj.CRS(ds.GetProjection()).to_epsg()
            proj   = ccrs.epsg(epsg) if int(epsg) != 4326 else ccrs.PlateCarree()

            arr = ds.ReadAsArray()
            WESN = bbGIS.get_extent(ds)
            ax.set_extent(WESN, crs=proj)
            ax.add_image(self.basemap, 10)

            ## add the interpolated GPS rates
            norm = mpl.colors.BoundaryNorm(np.arange(np.ceil(arr.min()), np.floor(arr.max())+0.5, 0.5), 256)
            im   = ax.imshow(arr, cmap=self.lst_parms[0], norm=norm,
                            zorder=20, interpolation='nearest', alpha=ALPHA,
                            extent=WESN, transform=proj, origin='upper')
            bbPlot.cartopy_cbar(im, f'V$_{{{lbl}}}$ (mm/yr)')

            ax.set_title(f'Interpolated MIDAS GPS {self.reg}')
            fig.set_label(f'GPSInterp_{self.reg}_{vcol}')

            ## add the GPS
            bbPlot.plot_dfll(df_gps, c=df_gps[vcol], axes=ax, norm=norm, draw_grid=False, **msty)

            left = not i
            gl = ax.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, left=left, bottom=True)
            logger.info(f'Interpolated GPS: {np.nanmean(arr):.2f} +/- {np.nanstd(arr):.2f} mm/yr')

        return fig, axes


    def plot_GPS_horizonatal(self):
        """ Plot the GPS locations with horizontal uncertainty """
        lst_ds, names = [ds_e, ds_n], 'East North'.split()

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs}, ncols=2)
        for i, ax in enumerate(axes):
            ax.set_extent(self.WESN, crs=self.proj)
            ax.add_image(self.basemap, 10)
            col = 'e_vel' if i == 0 else 'n_vel'
            bbPlot.plot_dfll(self.df_gps, c=self.df_gps[col], s=100, zorder=100, axes=ax)

        return fig, axes


## same as previous except no dct exp
class PlotVLM2(bzBase):
    def __init__(self, region, df_los=None):
        super().__init__()
        self.reg       = region
        self.lst_parms = DCT_PARMS[self.reg]
        # self.basemap = cimgt.GoogleTiles(url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
        # self.basemap   = cimgt.GoogleTiles(url='https://server.arcgisonline.com/'\
                        # 'arcgis/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}.jpg')
        # self.basemap   = cimgt.GoogleTiles(url='https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}')

        self.basemap   = cimgt.GoogleTiles(url='https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}')
        self.norm    = mpl.colors.BoundaryNorm(boundaries=self.lst_parms[1], ncolors=256)
        self.proj    = ccrs.PlateCarree()
        # self.SNWE    = regions[region]
        self.SNWE = [40.530, 40.822960, -74.170, -73.756498] \
                            if self.reg == 'NYC' else regions[region]
        self.WESN    = [*self.SNWE[2:], *self.SNWE[0:2]]

        path_gps     = op.join(self.path_unr, f'Midas_GPS-{self.reg}_cln.csv')
        self.df_gps  = prep_gps(path_gps, self.reg)
        self.df_gps  = bbGIS.in_dfll(self.df_gps, WESN=self.WESN) # some regions smaller
        self.df_gps['u_vel'] *= 1000
        self.df_gps['u_sig'] *= 1000


    def __call__(self):
        self.plot_GPS()
        self.plot_interp_GPS()
        self.plot_imaging_GPS_grid()
        # self.plot_GIA('ICE6G')
        # self.plot_GIA('Caron')


    def plot_GPS(self, axes=None):
        """ Plot the MIDAS GPS locations (+ rates and names if axes=None) """
        msty = dict(marker='o', s=50, linewidth=1.5, transform=self.proj,
                      cmap=self.lst_parms[0], norm=self.norm, alpha=1, zorder=30, edgecolors='white')
        WESN = self.WESN
        if axes is None:
            fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
            axes.set_extent(WESN, crs=self.proj)
            axes.add_image(self.basemap, 13)
            s = axes.scatter(self.df_gps['lon'], self.df_gps['lat'], c=self.df_gps['u_vel'], **msty)

            bbPlot.cartopy_cbar(s, 'V$_{up}$ (mm/yr)')
            axes.set_title(f'MIDAS GPS {self.reg}')
            fig.set_label(f'Midas_{self.reg}')

            ## format the names
            sty_name = dict(facecolor='sandybrown', alpha=0.3, boxstyle='round')
            names = self.df_gps['sta'].astype('str').values
            # Add text x pixels of the station.
            for i in range(len(names)):
                geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=+5+i)
                axes.text(self.df_gps['lon'].iloc[i], self.df_gps['lat'].iloc[i],
                        names[i].strip(), transform=text_transform, zorder=20,
                        verticalalignment='center', horizontalalignment='right',
                        bbox=sty_name)
            self.log.info(f'Raw GPS: {self.df_gps.u_vel.mean():.2f} +/- {self.df_gps.u_sig.mean():.2f} mm/yr')
            gl = axes.gridlines(draw_labels=True)
            bbPlot.fmt_gridlines(gl, bottom=True)
        else:
            s = axes.scatter(self.df_gps['lon'], self.df_gps['lat'], c=self.df_gps['u_vel'], **msty)
            fig = s.axes.get_figure()
        return fig, axes


    def plot_interp_GPS(self):
        """ Interpolate and plot the MIDAS GPS data """
        ds  = bbGIS.pts2grd(self.df_gps, 'u_vel', '4087', 90, 90, self.SNWE)
        arr = ds.ReadAsArray()

        WESN   = bbGIS.get_extent(ds)
        epsg   = pyproj.CRS(ds.GetProjection()).to_epsg()
        proj   = ccrs.epsg(epsg) if int(epsg) != 4326 else ccrs.PlateCarree()

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent(WESN, crs=proj)
        axes.add_image(self.basemap, 10)

        ## add the interpolated GPS rates
        im   = axes.imshow(arr, cmap=self.lst_parms[0], norm=self.norm,
                        zorder=20, interpolation='nearest', alpha=0.5,
                        extent=WESN, transform=proj, origin='upper')
        bbPlot.cartopy_cbar(im, 'V$_{up}$ (mm/yr)')

        axes.set_title(f'Interpolated MIDAS GPS {self.reg}')
        fig.set_label(f'GPSInterp_{self.reg}')

        ## add the GPs
        self.plot_GPS(axes)
        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)
        self.log.info(f'Interpolated GPS: {np.nanmean(arr):.2f} +/- {np.nanstd(arr):.2f} mm/yr')
        return fig, axes


    def plot_imaging_GPS_grid(self):
        """ Plot the gridded UNR imagery for study area """
        da     = load_gridded_imaging()
        ## just cropping the plot
        da_reg = da.sel(lat=slice(np.floor(self.SNWE[0]), np.ceil(self.SNWE[1])),
                        lon=slice(np.floor(self.SNWE[2]), np.ceil(self.SNWE[3])))

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent(self.WESN, crs=self.proj)
        axes.add_image(self.basemap, 10)

        ## add the interpolated GPS rates
        im   = axes.pcolormesh(da.lon, da.lat, da, shading='nearest', alpha=ALPHA,
                        transform=self.proj, cmap=self.lst_parms[0], norm=self.norm)
        bbPlot.cartopy_cbar(im, 'V$_{up}$ (mm/yr)')
        axes.set_title(f'Imaging (Grid) {self.reg} V$_{{up}}$')
        fig.set_label(f'Imaging_grid_{self.reg}')

        ## add the GPs
        self.plot_GPS(axes)

        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)
        self.log.info(f'Imaging: {da_reg.mean().item():.2f} +/- {da_reg.std().item():.2f} mm/yr')
        return fig, axes


    def plot_GIA(self, model='ICE6G'):
        """ Plot regional GIA from ICE6G or Caron """
        da     = load_GIA(model).load()
        da_reg = da.sel(lat=slice(self.SNWE[0], self.SNWE[1]),
                        lon=slice(np.floor(self.SNWE[2]), self.SNWE[3])).sortby('lat lon'.split())

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent(self.WESN, crs=self.proj)
        axes.add_image(self.basemap, 10)

        ## add the interpolated GPS rates
        im   = axes.pcolormesh(da.lon, da.lat, da, shading='nearest', alpha=ALPHA,
                        transform=self.proj, cmap=self.lst_parms[0], norm=self.norm)
        bbPlot.cartopy_cbar(im, 'V$_{up}$ (mm/yr)')
        axes.set_title(f'GIA: {model}')

        ## add the GPs
        self.plot_GPS(axes)

        fig.set_label(f'GIA_{model}')

        ## this doesnt work right
        # fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        # axes.set_extent(WESN, crs=self.proj)
        # axes.add_image(self.basemap, 10)
        # quad = da_reg.plot.pcolormesh(alpha=0.5, subplot_kws={'projection': self.basemap.crs}, zorder=20)
        # ax   = quad.axes
        # ax.add_image(self.basemap, 10)
        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)
        self.log.info(f'{model}: {da_reg.mean().item():.2f} +/- {da_reg.std().item():.2f} mm/yr')
        return fig, axes


class PlotVLM_Manoo(PlotVLM):
    def __init__(self, exp):
        super().__init__(exp.dct_exp, exp.mp_exp0, exp.ref_sta, exp.neofs)
        # assert self.reg in 'NYC MD Houston'.split(), f'Unsupported region: {self.reg}'
        # dct_mregs  = {'Houston': 'gulf_vlm.csv', 'NYC': 'NY_v3_manoo.txt', 'MD':'MD.txt'}
        # dct_mregs  = {'Houston': 'gulf_vlm.csv', 'NYC': 'NY.txt'}
        path_manoo = Path(self.path_vlm) / 'Other' / 'Manoo'
        path_vlm = path_manoo / 'EC_vlm.csv' # Atlantic Coast
        # path_vlm = path_manoo / 'gulf_vlm.csv' # Gulf Coast
        # path_vlm = path_manoo / 'NY.txt'
        sep = '\\s+' if path_vlm.suffix == '.txt' else ','

        df = pd.read_csv(path_manoo / path_vlm, sep=sep, header=None)
        df.columns = 'lon lat rate unc'.split()
        df['rate'] *= 10
        df['unc']  *= 10
        W,E,S,N=-74.5, -73.0, 40.0, 41.5

        # self.df_reg = bbGIS.in_dfll(df, SNWE=self.SNWE)
        self.df_reg = bbGIS.in_dfll(df, SNWE=[S,N,W,E])

        m, l, h = self.df_reg.rate.median(), self.df_reg.rate.quantile(0.003), self.df_reg.rate.quantile(0.997)
        print (f'{m.item():.2f} [{l.item():.2f}, {h.item():.2f}]')
        print (f'99% Range: {h-l:.2f} mm/yr')
        print (f'Full Range: {self.df_reg.rate.min():.2f} to {self.df_reg.rate.max():.2f} mm/yr')
        self.plot_df()
        self.plot_gps_vs_insar_qq()
        # self.plot_city('New York')


    def plot_city(self, city='Manhattan', crop=True):
        """ Plot Manoo VLM wihtin a city; only supports close up (crop=True)"""
        from BZ.bbGIS import city_shape
        city_wgs, city_utm = city_shape(city)
        gdf_reg = bbGIS.df2gdf(self.df_reg, 'all')
        gdf_city = gdf_reg[gdf_reg.geometry.within(city_wgs.geometry.iloc[0])]

        cmap = 'cmc.vik'
        norm = mpl.colors.Normalize(-4, 4)
        proj = ccrs.PlateCarree()

        ## add the interpolated insar rates
        fig, axes = bbPlot.plot_dfll(gdf_city, c=gdf_city['rate'], zoom=13,
                        use_dark=True, adj_fact=0.1, cmap= cmap, norm=norm)
                        #  cmap='cmc.oslo', s=1,
        # axes.set_extent(W,E,S,N, crs=proj)
        # cb.axes.set_yscale('linear')

        # add the outline of the city
        sty  = dict(zorder=50, alpha=0.9, facecolor='none', edgecolor='white',
        linestyle=':', linewidth=1.0)
        axes.add_geometries(city_wgs.geometry, crs=proj, **sty)
        fig.set_label(f'Manoo_Vup_{city}')


    def plot_df(self):
        """ Plot the dataframe as points """
        parms = DCT_PARMS.get(self.reg)
        if self.reg == 'Houston':
            # norm = mpl.colors.Normalize(-10, 0)
            norm = mpl.colors.TwoSlopeNorm(0, -15, 5)
            # cmap = 'cmo.ice'
            cmap = 'cmc.vik'
        elif self.reg == 'NYC':
            # cmap = 'cmc.vik'
            cmap = 'jet_r'
            # norm = mpl.colors.BoundaryNorm(np.arange(-2, 3, 0.2), 256)
            # norm = mpl.colors.Normalize(-5, 5)
            norm = mpl.colors.BoundaryNorm(np.arange(-3.5, 4.0, 0.5), 256)
        else:
            cmap = 'cmc.vik'
            norm = mpl.colors.Normalize(-3, 3)

        f, a = bbPlot.plot_dfll(self.df_reg, c=self.df_reg['rate'], zoom=10, use_dark=False,
                        #  cmap='cmc.vik', s=1, norm=mpl.colors.TwoSlopeNorm(0, -5, 5))
                        #  cmap=cmap, norm=mpl.colors.TwoSlopeNorm(0, -3, 3))
                        #  cmap='cmc.oslo', norm=mpl.colors.Normalize(-2.25, -0.25))
                        cmap= cmap, norm=norm)
        cbar = a.get_children()[0].colorbar
        # cbar.set_label('Vertical Rate (mm/yr)', fontsize=13)
        # cb.axes.set_yscale('linear')
        WESN = [*self.SNWE[2:], *self.SNWE[:2]]
        # a.set_extent(WESN)
        a.set_extent([-74.5, -73.0, 40.0, 41.5])
        f.set_label(f'Manoo_{self.reg}')

        self.lst_parms[0] = cmap
        self.norm = norm
        # self.plot_GPS(axes=a)

        # import seaborn as sns
        # fig, axes = plt.subplots(figsize=(12, 8))
        # sns.kdeplot(self.df_reg['rate'], fill=True, color='gray', ax=axes)
        # bbPlot.plot_hist(self.df_reg['rate'])


    def plot_gps_vs_insar_qq(self, radius=30, llim=-20, ulim=10):
        gdf_gps = bbGIS.df2gdf(self.df_gps, 'all').to_crs(4087)
        gdf_vlm = bbGIS.df2gdf(self.df_reg, 'all').to_crs(4087)
        dct_cmp = {}
        for ix, row in gdf_gps.iterrows():
            gdf_poly = gpd.GeoDataFrame(geometry=gpd.GeoSeries(row['lalo'].buffer(radius), crs=4087))
            gdf_vlm_at_gps = gdf_vlm.sjoin(gdf_poly, predicate='intersects')
            if gdf_vlm_at_gps.empty:
                print (f'No VLM data at {row.sta}')
                arr_vel, arr_unc = np.nan, np.nan
            arr_vel, arr_unc = np.nanmean(gdf_vlm_at_gps['rate']), np.nanmean(gdf_vlm_at_gps['unc'])
            dct_cmp[row.sta] = [row.u_vel, row.u_sig, arr_vel, arr_unc]

        df = pd.DataFrame(dct_cmp).T
        df.columns='gps_vel gps_unc ins_vel ins_unc'.split()
        df['resid'] = df['gps_vel'] - df['ins_vel']
        n0 = df.shape[0]
        df = df.dropna(subset=['resid'])
        log.warning (f'Dropped {n0-df.shape[0]} of {n0} GPS stations without nearby InSAR.')

        ## this is copied from analyGPS
        df_cmp = df
        ser_resid = df_cmp['resid']

        nn = self.df_gps.shape[0]
        log.info(f'Considering {ser_resid.shape[0]} stations of {nn} stations')
        rmse = ser_resid.abs().mean()
        r2   = calc_r2(df_cmp['ins_vel'], df_cmp['gps_vel'])

        fs = 14
        fig, axes = plt.subplots(figsize=(10,6))
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

        return df


def plot_davy(reg='DC'):
    """ Script for plotting the global dataset from machine learning model

    https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GL104497
    """
    import rioxarray as xrr
    from shapely.geometry import box
    import geopandas as gpd
    src = Path(os.getenv('dataroot')) / 'VLM' / 'ds03.tif'
    da = -1*xrr.open_rasterio(src).sel(band=1) # i think they mask uplift
    s, n, w, e = regions[reg]
    gser_reg = gpd.GeoSeries(box(w, s, e, n))
    da_reg = da.rio.clip(gser_reg.geometry, crs=4326, all_touched=True)
    da_reg = da_reg.where(da_reg<1e3) #

    proj    = ccrs.PlateCarree()
    basemap = cimgt.GoogleTiles(
            url='https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}')
            # url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
            # url='https://server.arcgisonline.com/ArcGIS/rest/services/' \
            # 'Elevation/World_Hillshade_Dark/MapServer/tile/{z}/{y}/{x}.jpg')


    cmap = 'cmc.vik'
    # norm = mpl.colors.TwoSlopeNorm(0, -4, 4)
    norm = mpl.colors.Normalize(-4, 4)

    fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': basemap.crs})
    axes.set_extent([w,e,s,n], crs=proj)
    gl = axes.gridlines(draw_labels=True)
    bbPlot.fmt_gridlines(gl, bottom=True, size=9)
    axes.add_image(basemap, 13)
    im   = axes.pcolormesh(da_reg.x, da_reg.y, da_reg, shading='nearest',
                    transform=proj, norm=norm, cmap=cmap)

    bbPlot.cartopy_cbar(im, ylabel='Vertical Rate (mm/yr)', fontsize=10, pad=0.15)
    fig.set_label(f'Vup_{reg}_davydzenska')

    ## plot the histogram
    plt.figure()
    da_reg.plot.hist()
    bbPlot.da_mmm(da_reg)
    print (f'Rounded unique values:', np.unique(np.round(da_reg.data, 2)))
    # breakpoint()
    return


def rate_oelsman_tg(reg='DC', psmsl_id=360):
    """ Print the VLM rate at a PSMSL tide gauge """
    f = 'Vertical_land_motion_preliminary_point_wise_trends.nc'
    src = Path(os.getenv('dataroot')) / 'VLM' / 'Other' / f
    assert src.exists(), 'Cant find Oelsmann data'
    ds = xr.open_dataset(src)
    ds_tg = ds.set_index(x='ID').sel(x=str(psmsl_id))
    vlm_rate_tg = ds_tg['trend'].item() * 1000
    vlm_unc_tg = ds_tg['trend_un'].item() * 1000
    print (f'Rate at TG: {psmsl_id}= {vlm_rate_tg:.2f} +/- {vlm_unc_tg:.2f} (mm/yr)')
    return


def plot_oelsman_reg(reg='DC'):
    """ Print the VLM rate at a PSMSL tide gauge """
    f = 'VLM_reconstruction.nc'
    src = Path(os.getenv('dataroot')) / 'VLM' / 'Other' / f
    assert src.exists(), 'Cant find Oelsmann data'
    ds = xr.open_dataset(src).rename(VLM_trend_coefficient_mean='VLM_rate')
    da = ds['Full_VLM_reconstruction_mean']
    df_rates = ds['VLM_rate'].to_dataframe()
    # df_rates_reg = bbGIS.in_dfll(df_rates, SNWE=regions[reg])
    ## plotting here US to check if there is data available in DC
    # there is not.
    df_rates_reg = bbGIS.in_dfll(df_rates, SNWE=regions['US'])
    df_rates_reg['VLM_rate'] *= 1000
    norm = mpl.colors.Normalize(-5, 5)
    f, a = bbPlot.plot_dfll(df_rates_reg, c=df_rates_reg['VLM_rate'], zoom=5,
                    adj_fact=1, use_dark=True, cmap='cmc.vik', norm=norm)
    # check if the data hits somewhere
    a.scatter(-77.021667, 38.873333, marker='*', transform=ccrs.PlateCarree(), s=100)
    bbPlot.plot_bbox(regions[reg], axes=a)
    return


def plot_gia_horiz(reg, buff=0.05, show_reg=True):
    src = bzBase().path_vlm / 'GIA' / 'Hvel.12mgrid_512.nc'
    ds = xr.open_dataset(src).rename(Lat='lat', Lon='lon',
                                           East_250='east', North_250='north')
    ds  = ds.assign_coords({'lon':
                    (((ds.lon + 180) % 360) - 180)}).sortby('lat lon'.split())
    S, N, W, E = regions[reg]
    S-=buff; N+=buff; W-=buff; E+=buff
    ds_crop = ds.sel(lat=slice(S, N), lon=slice(W, E))
    da_e, da_n = ds_crop['east'], ds_crop['north']
    units = da_e.units
    da_mmm(da_e, prefix='East: ', units=units)
    da_mmm(da_n, prefix='North: ', units=units)

    ## plot params
    basemap = cimgt.GoogleTiles(
            url='https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}')
            # url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
            # url='https://server.arcgisonline.com/ArcGIS/rest/services/' \
            # 'Elevation/World_Hillshade_Dark/MapServer/tile/{z}/{y}/{x}.jpg')

    norm = mpl.colors.TwoSlopeNorm(0, -1, 1)
    cmap, proj, alpha = 'cmc.roma_r', ccrs.PlateCarree(), 0.7

    fig, axes = plt.subplots(figsize=(12, 9), ncols=2, sharey=True,
                             subplot_kw={'projection': proj})
    for i, (ax, da) in enumerate(zip(axes, [da_e, da_n])):
        im = ax.pcolormesh(da.lon, da.lat, da, transform=proj,
                    norm=norm, cmap=cmap, alpha=alpha)
        ax.add_image(basemap, 10)
        ax.set_extent([W, E, S, N])
        cbar = bbPlot.cartopy_cbar(im, xlabel=units)
        if i == 0:
            cbar.remove()
            left, ti = True, 'East Velocity'
        else:
            left, ti = False, 'North Velocity'
        ax.set_title(ti, fontsize=13)
        gl = ax.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, left=left, bottom=True, size=11)
        if show_reg:
            ax.add_geometries(box(W,S,E,N), crs=proj, edgecolor='k', linestyle='--', facecolor='none')

    fig.subplots_adjust(wspace=0.01)
    fig.set_label(f'GIA_Horiz_{reg}')
    return fig, axes



if __name__ == '__main__':
    # PlotVLM(HR_SR)()
    # PlotVLM(Charleston_SR)()

    # PlotVLM(HR_SR)()
    # PlotVLM(NNJ_SR)()
    # Obj = PlotVLM(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20)
    # Obj = PlotVLM(DC_SR, 'ERA5_SET_PM_ex_Fast_2017', neofs=15)
    # Obj = PlotVLM(Houston_SR, 'ERA5_SET_PM_ex_Fast', neofs=5)

    # Obj.plot_mu()
    # Obj.plot_GIA('ICE6G')
    # Obj.plot_GIA('Caron')
    # Obj.plot_imaging_GPS_grid()
    # Obj.plot_interp_GPS_horizontal()
    # Obj()

    # Exp = ExpBase(DC_SR, 'ERA5_SET_PM_Fast', 'USN8', neofs=20)
    Exp = ExpBase(Houston_SR, 'ERA5_SET_PM_Fast', 'NASA', neofs=5)
    # Exp = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_Fast', 'NYBK', neofs=20)
    PlotVLM_Manoo(Exp)

    # plot_davy()
    # rate_oelsman_tg()
    # plot_oelsman_reg()

    # bbPlot.savefigs(Obj.path_res, True, True)
    plt.show()
