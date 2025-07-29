from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase

import xarray as xr
import seaborn as sns
from collections import OrderedDict
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict

from BZ import bbPlot
from BZ import *


def close_up_EC():
    """ Zoomed in map with NYC, HR, Charleston """
    # sty = dict(marker=7, s=750, c='white', cmap='coolwarm', zorder=10)
    # sty = dict(marker='s', s=750, facecolor='none', edgecolor='white', cmap='coolwarm', zorder=10)
    sty = dict(marker='.', s=750, c='white', cmap='coolwarm', zorder=10)


    SNWE       = regions['EC']
    S, N, W, E = [float(i) for i in SNWE];
    WESN       = bbPlot._adjust_llbounds([W,E,S,N], 0)
    proj0      = ccrs.LambertCylindrical()
    proj1      = ccrs.PlateCarree()
    proj2      = ccrs.Sinusoidal()

    proj      = proj1
    url       = 'https://server.arcgisonline.com/arcgis/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}.jpg'

    fig, axes = plt.subplots(figsize=(14,8), subplot_kw={'projection': proj})
    basemap   = cimgt.GoogleTiles(url=url)
    axes.add_image(basemap, 8)


    cities = 'Charleston', 'Hampton Roads', 'NYC'
    lats   = [32.7765, 36.9339, 40.7128]
    lons   = [-79.9311, -76.3637, -74.006]


    s   = axes.scatter(lons, lats, transform=ccrs.PlateCarree(), **sty)

    sty_name = dict(facecolor='white', alpha=0.5, boxstyle='round')
    yxs   = [(+5, +100), (-20, +60), (-20, -10)]
    names = cities
    # Add text x pixels of the station.
    for i in range(len(names)):
        geoT = ccrs.PlateCarree()._as_mpl_transform(axes)
        txtT = mpl.transforms.offset_copy(geoT, units='dots', y=yxs[i][0], x=yxs[i][1])
        axes.text(lons[i], lats[i], names[i].strip(), transform=txtT,
                verticalalignment='center', horizontalalignment='right', zorder=20,
                bbox=sty_name, fontsize=15)

    gl = axes.gridlines(draw_labels=True)
    bbPlot.fmt_gridlines(gl, bottom=True, size=15)
    fig.patch.set_alpha(0)


    fmt      = 'pdf'
    path_fig = op.join(PATH_RES, 'Figures', f'Overview_Map.{fmt}')
    plot_kws = dict(dpi=600, bbox_inches='tight', pad_inches=0.025, transparent=False,
                edgecolor=fig.get_edgecolor(), format=fmt)

    fig.savefig(path_fig, **plot_kws)
    print ('Saved:', path_fig)
    return


def plot_globe(clon=-74, clat=38):
    """ Plot an inset globe (for e.g., showing where EC is) """
    url  = 'https://server.arcgisonline.com/arcgis/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}.jpg'

    proj = ccrs.NearsidePerspective(clon, clat)
    fig, axes = plt.subplots(figsize=(14,8), subplot_kw={'projection': proj})
    basemap   = cimgt.GoogleTiles(url=url)
    # basemap   = cimgt.Stamen('terrain-background')
    axes.add_image(basemap, 2)
    axes.coastlines()

    cities = 'Charleston', 'Hampton Roads', 'NYC'
    lats   = [32.7765, 36.9339, 40.7128]
    lons   = [-79.9311, -76.3637, -74.006]
    sty = dict(marker='.', s=500, c='white', edgecolor='r', zorder=10)

    s   = axes.scatter(lons, lats, transform=ccrs.PlateCarree(), **sty)
    axes.set_global()

    fig.set_label('Overview_EC')
    return


def plot_houston(clon=-94.18-3, clat=31.7953):
    """ Plot an inset globe (for e.g., showing where EC is) """
    from matplotlib.patches import Rectangle

    url='https://www.google.cn/maps/vt?lyrs=s@189&gl=cn&x={x}&y={y}&z={z}'
    # url='https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png')
    # url='https://server.arcgisonline.com/ArcGIS/rest/services/' \
    # 'Elevation/World_Hillshade_Dark/MapServer/tile/{z}/{y}/{x}.jpg')
    # url  = 'https://server.arcgisonline.com/arcgis/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}.jpg'

    # proj = ccrs.NearsidePerspective(clon, clat, satellite_height=1.95e6)
    proj = ccrs.PlateCarree()

    fig, axes = plt.subplots(figsize=(14,8), subplot_kw={'projection': proj})
    basemap   = cimgt.GoogleTiles(url=url)
    axes.add_image(basemap, 5)
    # axes.coastlines()

    cities = ['Houston']
    lats = [29.7550]
    lons = [-95.3621]
    sty = dict(marker='*', s=500, c='white', edgecolor='r', zorder=10)

    s   = axes.scatter(lons, lats, transform=ccrs.PlateCarree(), **sty)
    # axes.set_global()
    axes.set_extent([clon-20, clon+20, clat-20, clat+20], crs=ccrs.PlateCarree())
    border = Rectangle((0, 0), 1, 1, transform=axes.transAxes, color='w',
                       linewidth=3, fill=False)
    axes.add_patch(border)
    # fig.set_label('Overview_Houston')
    return


if __name__ == '__main__':
    # close_up_EC()
    # plot_globe()
    plot_houston()
    bbPlot.savefigs(PATH_RES, True, True)
    plt.show()