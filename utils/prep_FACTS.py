""" Take a raster and convert it to points lat/lon vlm vlm_unc """
from BZ import *
from BZ import bbPlot
from VLM import *
from BZ.bbGIS import city_shape
import cmcrameri
import rioxarray as xrr
import geopandas as gpd

PATH_WD = Path(os.getenv('dataroot')) / 'VLM' / 'Other' / 'FACTS'
CMAP = 'cmc.vik'
NORM = mpl.colors.TwoSlopeNorm(0, -5, 5)

def NYC(show=True):
    lst_dfs = []
    for i, k in enumerate('rate std'.split()):
        name = 'rate' if i == 0 else 'std'
        ds = xr.open_dataset(PATH_WD / f'NYC_{k}.nc')
        da = ds['Band1']*1000
        df = da.rename(name).to_series().reset_index().dropna().reset_index(drop=True)
        lst_dfs.append(df)
        if i == 0 and show:
            fig = plt.figure()
            da.plot(cmap=CMAP, norm=NORM)
    df = lst_dfs[0]
    df['std'] = lst_dfs[1]['std']
    df = df.rename(columns={'x': 'lon', 'y': 'lat'})
    df = df['lon lat rate std'.split()]
    dst = PATH_WD / 'NYC_rate_unc.csv'
    df.to_csv(dst)
    if show:
        bbPlot.plot_dfll(df, c=df['rate'], cmap=CMAP, norm=NORM)
    print (f'Wrote: {dst}')
    return df


def CA(show=True):
    gdf_crop = gpd.read_file(PATH_WD / 'Socal_Coast.GeoJSON')
    lst_dfs = []
    for i, k in enumerate(['', '_std']):
        name = 'rate' if i == 0 else 'std'
        src = PATH_WD / f'CA_VLM{k}.tif'
        da = xrr.open_rasterio(src).squeeze()
        da.rio.write_crs(4326, inplace=True)
        da_crop = da.rio.clip(gdf_crop.geometry, all_touched=True, drop=True)
        df_crop = da_crop.rename(name).to_series().reset_index().dropna()
        lst_dfs.append(df_crop)
        if i == 0 and show:
            plt.figure()
            da_crop.plot(cmap=CMAP, norm=NORM)

    df_crop = lst_dfs[0]
    df_crop['std'] = lst_dfs[1]['std']
    df_crop = df_crop.rename(columns={'x': 'lon', 'y': 'lat'})
    df_crop = df_crop['lon lat rate std'.split()]
    dst = PATH_WD / 'Socal_rate_unc.csv'
    df_crop.reset_index(drop=True).to_csv(dst)
    print (f'Wrote: {dst}')
    if show:
        bbPlot.plot_dfll(df_crop, c=df_crop['rate'], cmap=CMAP, norm=NORM)
    return df_crop


def DC():
    from VLM.bzFRInGE.plotting.plotDC import make_crop
    gdf_crop = make_crop()
    lst_dfs = []
    for i, k in enumerate('rate std'.split()):
        ds = xr.open_dataset(PATH_WD / f'DC_{k}.nc')
        da = ds['Band1']*1000
        da.rio.write_crs(4326, inplace=True)
        da_crop = da.rio.clip(gdf_crop.geometry, all_touched=True, drop=True)
        df_crop = da_crop.rename(k).to_series().reset_index().dropna()
        lst_dfs.append(df_crop)
        if i == 0 and show:
            plt.figure()
            da_crop.plot(norm=NORM, cmap=CMAP)

    df_crop = lst_dfs[0]
    df_crop['std'] = lst_dfs[1]['std']
    df_crop = df_crop['lon lat rate std'.split()]
    dst = PATH_WD / 'DC_rate_unc.csv'
    df_crop.reset_index(drop=True).to_csv(dst)
    print (f'Wrote: {dst}')

    if show:
        bbPlot.plot_dfll(df_crop, c=df_crop['rate'], cmap=CMAP, norm=NORM)
    return df_crop


def Houston(show=True):
    from VLM.bzFRInGE.plotting.plotDC import make_crop
    lst_dfs = []
    for i, k in enumerate('rate std'.split()):
        ds = xr.open_dataset(PATH_WD / f'Houston_{k}.nc')
        da = ds['Band1']*1000
        da.rio.write_crs(4326, inplace=True)
        df = da.rename(k).to_series().reset_index().dropna()
        lst_dfs.append(df)
        if i == 0 and show:
            plt.figure()
            da.plot(norm=NORM, cmap=CMAP)

    df = lst_dfs[0]
    df['std'] = lst_dfs[1]['std']
    df = df['lon lat rate std'.split()]
    dst = PATH_WD / 'Houston_rate_unc.csv'
    df.reset_index(drop=True).to_csv(dst)
    print (f'Wrote: {dst}')

    if show:
        bbPlot.plot_dfll(df, c=df['rate'], cmap=CMAP, norm=NORM)
    return df



if __name__ == '__main__':
    # NYC()
    # CA()
    Houston()
    # DC()
    plt.show()