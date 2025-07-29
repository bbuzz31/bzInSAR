from VLM.bzFRInGE.FRINGEBase import *
from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from BZ import *
from BZ import bbGIS, bbPlot

from mintpy.objects.gnss import GNSS
from mintpy.utils import readfile
from mintpy.cli import save_gdal
import xarray as xr

import geopandas as gpd
import pyproj
import cartopy.crs as ccrs
from cartopy.io import img_tiles as cimgt

PATH_RES = op.join(op.expanduser('~'), 'Desktop', 'VLM', 'NYC_2023')
ALPHA    = 0.8

## region to cmap, colorbounds, custom function
DCT_PARMS = {#'HR': ['Blues_r', np.arange(-3, 0.5, 0.5)],
             #'Charleston': ['Blues_r', np.arange(-2, 0.25, 0.25)],
             #'Maine': ['cmc.vik', np.arange(-2, 2, 0.5)],
             #'Savannah': ['cmc.vik', np.arange(-2, 2, 0.5)],
             'NYC': ['Reds', np.arange(6, 8.25, 0.25)],
             'DC': ['cmc.vik', np.arange(-3, 4, 1)],
            }

def get_geom(Exp):
    """ Get the incidence and azimuth angles into xarrays """
    path_geom = Exp.path_geom_mp_geo
    path_inp  = op.dirname(path_geom)
    path_inc  = op.join(path_inp, 'incidenceAngle.nc')
    path_az   = op.join(path_inp, 'azimuthAngle.nc')
    if op.exists(path_inc):
        ds_inc = xr.open_dataset(path_inc)
        ds_az  = xr.open_dataset(path_az)
    else:
        cmd = f'{path_geom} -d incidenceAngle -o {path_inc} --of netCDF'
        save_gdal.main(cmd.split())

        cmd = f'{path_geom} -d azimuthAngle -o {path_az} --of netCDF'
        save_gdal.main(cmd.split())

        ds_inc = xr.open_dataset(path_inc)
        ds_az  = xr.open_dataset(path_az)

    return ds_inc['Band1'], ds_az['Band1']


def load_ICE6G():
    """ Get vertical and horizontal (dont have latter for Caron) of ICE6G """
    src = op.join(bzBase().path_vlm, 'GIA', f'drad.12mgrid_512.nc')
    assert op.exists(src), f'{src} doesnt exist'
    ds   = xr.open_dataset(src).rename(Lat='lat', Lon='lon', Drad_250='u_vel')
    ds_u = ds.assign_coords({'lon':
                            (((ds.lon + 180) % 360) - 180)}).sortby('lat lon'.split())

    src = op.join(bzBase().path_vlm, 'GIA', f'Hvel.12mgrid_512.nc')
    assert op.exists(src), f'{src} doesnt exist'
    ds    = xr.open_dataset(src).rename(Lat='lat', Lon='lon',
                                        East_250='e_vel', North_250='n_vel')
    ds_en = ds.assign_coords({'lon':
                            (((ds.lon + 180) % 360) - 180)}).sortby('lat lon'.split())


    return xr.merge([ds_en, ds_u])


## from mintpy utils.utils0
def enu2los(v_e, v_n, v_u, inc_angle, head_angle=None, az_angle=None):
    """Project east/north/up motion into the line-of-sight (LOS) direction defined by incidence/azimuth angle.
    Parameters: v_e        - np.ndarray or float, displacement in east-west   direction, east  as positive
                v_n        - np.ndarray or float, displacement in north-south direction, north as positive
                v_u        - np.ndarray or float, displacement in vertical    direction, up    as positive
                inc_angle  - np.ndarray or float, incidence angle from vertical, in the unit of degrees
                head_angle - np.ndarray or float, azimuth angle of the SAR platform along track direction
                             measured from the north with clockwise direction as positive, in the unit of degrees
                az_angle   - np.ndarray or float, azimuth angle of the LOS vector from the ground to the SAR platform
                             measured from the north with anti-clockwise direction as positive, in the unit of degrees
                             head_angle = 90 - az_angle
    Returns:    v_los      - np.ndarray or float, displacement in LOS direction, motion toward satellite as positive
    """
    # unite (los_)head/az_angle into (los_)az_angle
    if az_angle is None:
        if head_angle is not None:
            az_angle = heading2azimuth_angle(head_angle)
        else:
            raise ValueError(f'invalid az_angle: {az_angle}!')

    # project ENU onto LOS
    v_los = (  v_e * np.sin(np.deg2rad(inc_angle)) * np.sin(np.deg2rad(az_angle)) * -1
             + v_n * np.sin(np.deg2rad(inc_angle)) * np.cos(np.deg2rad(az_angle))
             + v_u * np.cos(np.deg2rad(inc_angle)))

    return v_los


def midas2los(Exp, en_only=False, check=False):
    df_gps   = prep_gps(Exp.path_gps, Exp.reg, horiz=True)
    da_inc, da_az = get_geom(Exp)

    rows = []
    for ix, row in df_gps.iterrows():
        e, n, u  = row.e_vel, row.n_vel, row.u_vel
        e_sig, n_sig, u_sig  = row.e_sig, row.n_sig, row.u_sig
        if en_only:
            u     = 0
            u_sig = 0

        lat, lon  = row.lat, row.lon
        tol = 0.5 # degrees to search
        # fails if outside of study area
        try:
            inc = da_inc.sel(lat=lat, lon=lon, method='nearest', tolerance=tol)
            az  = da_az.sel(lat=lat, lon=lon, method='nearest', tolerance=tol)
            rate_los = enu2los(e, n, u, inc, az_angle=az).item()
            unc_los  = enu2los(e_sig, n_sig, u_sig, inc, az_angle=az).item()
        except:
            rate_los, unc_los = np.nan, np.nan
            inc = np.nan
            ## optionally check to find out why the point is nan
            if check:
                from shapely.geometry import Polygon, Point
                S, N, W, E = bbGIS.get_extent_nc(da_inc, msg=False)
                poly = Polygon.from_bounds(W, S, E, N)
                ## probably its outside of the bounding box entirely
                if poly.contains(Point(lon, lat)):
                    # otherwise, it might just be outside of the data area
                    breakpoint()

        row['los_vel'] = rate_los
        row['los_sig'] = unc_los
        rows.append(row)
    df_gps_los = pd.DataFrame(rows).dropna()
    return df_gps_los


def GIA_to_LOS(Exp, ew_only=False):
    """
    Convert GIA to LOS using nearest non/nan angles

    Verified order is correct using u_vel
    """
    da_inc, da_az = get_geom(Exp)
    S, N, W, E = Exp.SNWE
    ## round to the nearest int?
    SN = np.floor(S-0.05), np.ceil(N+0.05)
    WE = np.floor(W-0.05), np.ceil(E+0.05)
    ds_GIA = load_ICE6G().sel(lat=slice(*SN), lon=slice(*WE))

    ## use geopandas to convert the GIA or you'll have nans because of the cropping
    df_GIA  = ds_GIA.to_dataframe().reset_index()
    gdf_GIA = bbGIS.df2gdf(df_GIA, 'all')
    gdf_inc = bbGIS.df2gdf(da_inc.to_dataframe().reset_index().dropna(), 'all')
    gdf_m  = gpd.sjoin_nearest(gdf_GIA.to_crs(4087), gdf_inc.to_crs(4087),
                            distance_col='distance').rename(columns={'Band1':'incAngle'})

    gdf_az = bbGIS.df2gdf(da_az.to_dataframe().reset_index().dropna(), 'all')
    gdf_m1 = gpd.sjoin_nearest(gdf_GIA.to_crs(4087), gdf_az.to_crs(4087)).to_crs(4326
                                            ).rename(columns={'Band1':'azAngle'})

    assert (gdf_m.index_right == gdf_m1.index_right).all(), \
        'Az/Inc locations dont match?'

    gdf_m['azAngle'] = gdf_m1['azAngle']

    arr_los = []
    for ix, row in gdf_m.iterrows():
        e = row['e_vel']
        n = row['n_vel']
        u = row['u_vel']
        if ew_only:
            # print ('Setting vup to 0')
            u = 0
        arr_los.append(enu2los(e, n, u, row['incAngle'], az_angle=row['azAngle']))

    gdf_m['u_vel1'] = arr_los
    gdf_m.to_crs(4326, inplace=True)

    ## need to pivot to reshape correctly
    gdf_m['lat'] = gdf_m.geometry.y
    gdf_m['lon'] = gdf_m.geometry.x
    df_m = pd.DataFrame(gdf_m).drop(columns='lalo').pivot_table(
                index='lat', columns='lon', values='u_vel1')

    da_GIA  = ds_GIA['u_vel'].copy()
    da_GIA_los = da_GIA.copy()
    da_GIA_los.data = df_m.to_numpy()

    return da_GIA_los


def sanity_check(geom, sta='NJHT'):
    """ Use MintPy to project the GPS time series and rate to LOS

    sanity_check(Exp.path_geom_mp_geo, Exp.ref_sta)
    """

    Obj = GPS(site=sta)
    # dates, timeseries displacement, timeseries und uncertainty in (m), loc of gps and ref if given
    ts_gps_los = Obj.read_gps_los_displacement(geom)

    # rate los (not MIDAS; for sanity check)
    vel_gps_los = Obj.get_gps_los_velocity(geom)[0]
    print (f'Velocity at {sta}: {vel_gps_los:.4f} (m/y)')
    return


def plot_ICE6G_enu(SNWE=None, reg=None):
    """ Plot the E/N/U Components of ICE6G """
    assert SNWE or reg, 'must pass some location info'
    S, N, W, E = regions[reg] if SNWE is None else SNWE
    WESN = W,E,S,N
    ds   = load_ICE6G()

    basemap   = cimgt.GoogleTiles(url='https://server.arcgisonline.com/'\
                   'arcgis/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}.jpg')
    proj    = ccrs.PlateCarree()
    norms   = (mpl.colors.BoundaryNorm(boundaries=np.arange(-0.1, 0.1, 0.02), ncolors=256),
               mpl.colors.BoundaryNorm(boundaries=np.arange(0, 0.2, 0.02), ncolors=256),
               mpl.colors.BoundaryNorm(boundaries=np.arange(-2.2, -1.2, 0.1), ncolors=256))
    cmaps   = ['coolwarm', 'Reds', 'Blues_r']
    lst_das = []
    for i, c in enumerate('e n u'.split()):
        da = ds[f'{c}_vel']
        da = da.sel(lat=slice(np.floor(S), np.ceil(N)), lon=slice(np.floor(W), np.ceil(E)))
        # print (da)

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': basemap.crs})
        axes.set_extent(WESN, crs=proj)
        axes.add_image(basemap, 10)

        im   = axes.pcolormesh(da.lon, da.lat, da, shading='nearest', alpha=ALPHA,
                        transform=proj, cmap=cmaps[i], norm=norms[i])

        bbPlot.cartopy_cbar(im, f'{c.title()} Vel (mm/yr)')
        axes.set_title(f'ICE6G_{c.upper()}')
        fig.set_label(f'ICE6G_{c.upper()}')
        print (f'Mean {da.name}: {da.mean().item():.2f}')
        if c in 'e n':
            lst_das.append(da)

    ## plot the horizontal motion
    da_horiz = np.sqrt(lst_das[0]**2 + lst_das[1]**2)
    fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': basemap.crs})
    axes.set_extent(WESN, crs=proj)
    axes.add_image(basemap, 10)

    norm = mpl.colors.BoundaryNorm(boundaries=np.arange(0, 0.25, 0.05), ncolors=256)
    im   = axes.pcolormesh(da_horiz.lon, da_horiz.lat, da_horiz, alpha=ALPHA,
                    shading='nearest', transform=proj, cmap='Reds', norm=norm)

    bbPlot.cartopy_cbar(im, f'Horizontal (mm/yr)')
    axes.set_title(f'ICE6G$_{{EN}}$')
    fig.set_label(f'ICE6G_EN')
    print (f'Mean EW: {da_horiz.mean().item():.2f}')


    return


def proj_gps_en(Exp, en_only=True, plot=False):
    """ Project just the horizontal componets of GPS to midas """

    df_gps_los = midas2los(Exp, en_only=en_only)
    xy      = 0.000277778 # ~ 30 m
    ds_geo  = bbGIS.pts2grd(df_gps_los, 'los_vel', '4326', xy, xy, Exp.SNWE)
    path_gps_los = f'{op.splitext(Exp.path_gps)[0]}_los.nc'
    gdal.Translate(path_gps_los, ds_geo,  format='netcdf')
    da = xr.open_dataset(path_gps_los)['Band1']
    os.remove(path_gps_los)

    ## now interp this ds to the correct posting
    da_inc = get_geom(Exp)[0]
    da_gps = da.interp_like(da_inc).rename('')

    if plot:
        plt.figure()
        (1000*da_gps).rename('Horizontal vLOS').plot()

    return da_gps


def fit_plane(da):
    """ Fit a plane to a 2D array """
    from mintpy.objects import ramp
    dat  = np.flipud(da.data)
    # meta = da.attrs if not da.attrs.empty else None
    meta = None
    dat_new, ramp = ramp.deramp(dat, ramp_type='linear', metadata=meta)

    da_out = da.copy()
    da_out.data = np.flipud(da.data)
    return da_out


class PlotLOS(ExpBase):
    def __init__(self, exp):
        super().__init__(exp.dct_exp, exp.mp_exp0, exp.ref_sta, exp.neofs)
        self.lst_parms = DCT_PARMS[self.reg]
        # self.basemap = cimgt.GoogleTiles(url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
        self.basemap   = cimgt.GoogleTiles(url='https://server.arcgisonline.com/'\
                        'arcgis/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}.jpg')
        self.norm    = mpl.colors.BoundaryNorm(boundaries=self.lst_parms[1], ncolors=256)

        self.proj    = ccrs.PlateCarree()
        self.WESN    = [*self.SNWE[2:], *self.SNWE[0:2]]

        self.df      = midas2los(self)
        # self.df      = bbGIS.in_dfll(self.df, WESN=self.WESN) # some regions smaller
        self.df['los_vel'] *= 1000
        self.df['los_sig'] *= 1000


    def __call__(self):
        # self.plot_GPS()
        # self.plot_interp_GPS()
        # self.plot_interp_GPS_EW()
        # self.plot_GIA()
        self.plot_GIA_plane()


    def plot_GPS(self, axes=None):
        """ Plot the MIDAS GPS locations (+ rates and names if axes=None) """
        basemap = cimgt.GoogleTiles(url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
        msty = dict(marker='o', s=50, linewidth=1.5, transform=self.proj,
                      cmap=self.lst_parms[0], norm=self.norm, alpha=1, zorder=30, edgecolors='white')
        fig  = None
        if axes is None:
            fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
            axes.add_image(self.basemap, 10)
            s = axes.scatter(self.df['lon'], self.df['lat'], c=self.df['los_vel'], **msty)

            bbPlot.cartopy_cbar(s, 'LOS (mm/yr)')
            axes.set_title(f'MIDAS GPS {self.reg}')
            fig.set_label(f'MidasLOS_{self.reg}')

            ## format the names
            sty_name = dict(facecolor='sandybrown', alpha=0.3, boxstyle='round')
            names = self.df['sta'].astype('str').values
            # Add text x pixels of the station.
            for i in range(len(names)):
                geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(axes)
                text_transform = mpl.transforms.offset_copy(geodetic_transform, units='dots', y=+5+i)
                axes.text(self.df['lon'].iloc[i], self.df['lat'].iloc[i],
                        names[i].strip(), transform=text_transform, zorder=20,
                        verticalalignment='center', horizontalalignment='right',
                        bbox=sty_name)
            print (self.df.dropna())
        else:
            s = axes.scatter(self.df['lon'], self.df['lat'], c=self.df['los_vel'], **msty)

        axes.set_extent(self.WESN, crs=self.proj)
        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)
        return fig, axes


    def plot_interp_GPS(self):
        """ Interpolate and plot the MIDAS GPS data """
        ds  = bbGIS.pts2grd(self.df, 'los_vel', '4087', 90, 90, self.SNWE)
        arr = ds.ReadAsArray()

        # basemap = cimgt.Stamen('terrain-background')
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
        fig.set_label(f'GPSInterpLOS_{self.reg}')

        ## add the GPs
        self.plot_GPS(axes)
        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)
        return fig, axes


    def plot_interp_GPS_EW(self):
        """ Plot horizontal MIDAS rates interpolated to InSASR WGS84 grid """
        da  = proj_gps_en(self)*1000
        norm = mpl.colors.BoundaryNorm(np.arange(7, 9.25, 0.25), ncolors=256)

        # basemap = cimgt.Stamen('terrain-background')
        S,N,W,E = bbGIS.get_extent_nc(da, False)
        epsg    = 4326
        proj    = ccrs.epsg(epsg) if int(epsg) != 4326 else ccrs.PlateCarree()

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent([W,E,S,N], crs=proj)
        axes.add_image(self.basemap, 10)

        im   = axes.pcolormesh(da.lon, da.lat, da, shading='nearest', alpha=ALPHA,
                        transform=proj, cmap=self.lst_parms[0])#, norm=self.norm)
        bbPlot.cartopy_cbar(im, 'MIDAS LOS (mm/yr)')

        axes.set_title(f'Interpolated MIDAS GPS (EW) {self.reg}')
        fig.set_label(f'GPSInterpLOS_{self.reg}_WGS84')

        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)

        ## add the GPS
        # self.plot_GPS(axes)
        return fig, axes


    def plot_GIA(self, ew_only=True):
        """ Plot regional GIA from ICE6G or Caron """
        da_reg = GIA_to_LOS(self, ew_only)
        if ew_only:
            # norm = mpl.colors.TwoSlopeNorm(0, -0.02, 0.02)
            norm = mpl.colors.TwoSlopeNorm(0, -0.1, 0.1)
            cmap = 'coolwarm'
        else:
            norm = mpl.colors.BoundaryNorm(np.arange(-1.6, -1.4, 0.02), ncolors=256)
            cmap = 'Blues_r'

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent(self.WESN, crs=self.proj)
        axes.add_image(self.basemap, 10)

        ## add the GIA
        im   = axes.pcolormesh(da_reg.lon, da_reg.lat, da_reg, shading='nearest', alpha=ALPHA,
                        transform=self.proj, cmap=cmap, norm=norm)
        bbPlot.cartopy_cbar(im, 'LOS (mm/yr)')
        axes.set_title(f'GIA: ICE6G')

        ## add the GPs
        # self.plot_GPS(axes)

        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)

        fig.set_label(f'GIA_ICE6G_LOS')
        return fig, axes


    def plot_GIA_plane(self, ew_only=True):
        """ Plot regional GIA from ICE6G or Caron """
        da_reg = GIA_to_LOS(self, ew_only)
        if ew_only:
            # norm = mpl.colors.TwoSlopeNorm(0, -0.02, 0.02)
            norm = mpl.colors.TwoSlopeNorm(0, -0.1, 0.1)
            cmap = 'coolwarm'
        else:
            norm = mpl.colors.BoundaryNorm(np.arange(-1.6, -1.4, 0.02), ncolors=256)
            cmap = 'Blues_r'

        da_plane = fit_plane(da_reg)

        fig, axes = plt.subplots(figsize=(10,10), subplot_kw={'projection': self.basemap.crs})
        axes.set_extent(self.WESN, crs=self.proj)
        axes.add_image(self.basemap, 10)

        ## add the GIA
        im   = axes.pcolormesh(da_reg.lon, da_reg.lat, da_reg, shading='nearest', alpha=ALPHA,
                        transform=self.proj, cmap=cmap, norm=norm)

        interactive = False
        if interactive:
            import hvplot.xarray
            plot = da_reg.hvplot(cmap='coolwarm')
            hvplot.show(plot)


        bbPlot.cartopy_cbar(im, 'LOS (mm/yr)')
        axes.set_title(f'GIA: ICE6G')

        ## add the GPs
        # self.plot_GPS(axes)

        gl = axes.gridlines(draw_labels=True)
        bbPlot.fmt_gridlines(gl, bottom=True)

        fig.set_label(f'GIA_ICE6G_LOS')
        return fig, axes


if __name__ == '__main__':
    # Exp = ExpBase(NYC_SRc, 'ERA5_SET_PM_Fast', 'NYBK', neofs=20)
    Exp = ExpBase(DC_SR, 'ERA5_SET_PM_Fast', 'USN8', neofs=20)

    # df_gps_los = midas2los(Exp)

    # plot_ICE6G_enu(Exp.SNWE)

    PlotLOS(Exp)()
    # proj_gps_en(Exp, plot=True)

    # bbPlot.savefigs(PATH_RES, False, True)
    plt.show()