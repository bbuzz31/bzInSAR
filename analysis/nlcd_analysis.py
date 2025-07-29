from contrib import *
from contrib.experiments import *
from contrib.FRINGEBase import ExpBase

import xarray as xr
import seaborn as sns
from collections import OrderedDict
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
from matplotlib.backends.backend_pdf import PdfPages

from BZ import bbPlot

# import log

# NLCD code: my code, my label name
DCT_LULC = {11: [0, 'Water'],
            21: [1, 'Developed'],
            22: [1, 'Developed'],
            23: [1, 'Developed)'],
            24: [1, 'Developed)'],
            # 21: [7, 'Dev. (Low)'],
            # 22: [7, 'Dev. (Low)'],
            # 23: [9, 'Dev. (High)'],
            # 24: [9, 'Dev. (High)'],
            41: [2, 'Forest'],
            42: [2, 'Forest'],
            43: [2, 'Forest'],
            51: [3, 'Shrub'],
            52: [3, 'Shrub'],
            71: [4, 'Ag'],
            81: [4, 'Ag'],
            82: [4, 'Ag'],
            90: [5, 'Wetlands'],
            95: [5, 'Wetlands'],
            31: [6, 'Barren'],
            0:  [-1, 'nan']
            }

DCT_LULC_LBL = {0:'Water', 1: 'Developed', 7:'Dev. (Low)', 9:'Dev. (High)', 2:'Forest', 3: 'Shrub',
                4:'Ag', 5:'Wetlands', 6:'Barren'}

# DCT_LULC_COL = {'Water': 'Darkblue', 'Dev.': 'Darkred', 'Forest': 'DarkGreen',
#                 'Shrub': 'Green', 'Ag': 'Gold', 'Wetlands': 'Blue', 'Barren': 'Brown'}
# DCT_LULC_COL = {k: mpl.colors.to_rgba(v) for k, v in DCT_LULC_COL.items()}

DCT_LULC_COL = {'Water': (0, 0, 255, 1), 'Developed': (255, 153, 0, 1),
               'Dev. (Low)': (255, 204, 0, 1), 'Dev. (High)': (255, 0, 0, 1),
               'Forest': (0, 102, 0, 1), 'Shrub': (178, 178, 0, 1), 'Ag': (204, 77, 128, 1),
               'Wetlands': (0, 255, 255, 1), 'Barren': (229, 229, 204, 1)}


DCT_VUP = OrderedDict([('fast_s', ('Fast Subsidence', '[-5, -3)')), ('mod_s', ('Subsidence', '[-3, -1)')),
            ('no_s', ('Stable', '[-1,1)')), ('mod_u', ('Uplift', '[1,3]')),
            ('fast_u', ('Fast Uplift', '[3,5) mm/yr'))])


for k, v in DCT_LULC_COL.items():
    try: # so strings are maintained
        vnew = [vv/255.0 for vv in v]
        vnew[-1] = 1
    except:
        vnew = v
    DCT_LULC_COL[k] = vnew

WCOL = 'dimgray'
PLFS = 20 # piechart fontsize


def cmp_NLCD_VLM(exp):
    da_nlcd0        = load_NLCD(exp)['NLCD']
    da_rate, da_std = load_rate_std(exp)
    da_nlcd1        = match_grid(da_rate, da_nlcd0)
    da_nlcd         = combine_lulc(da_nlcd1)

    DCT_SUB         = [-3, -1, 1, 1, 3]

    da_fast_s = xr.where(da_rate<=-3, 1, np.nan).rename('fast_s')
    da_mod_s  = xr.where(((da_rate>-3) & (da_rate<=-1)), 1, np.nan).rename('mod_s')
    da_null   = xr.where(((da_rate>-1) & (da_rate<=1)), 1, np.nan).rename('no_s')
    da_mod_u  = xr.where(((da_rate>1) & (da_rate<=3)), 1, np.nan).rename('mod_u')
    da_fast_u = xr.where(da_rate>3, 1, np.nan).rename('fast_u')

    ds_m = da_nlcd * xr.merge([da_fast_s, da_mod_s, da_null, da_mod_u, da_fast_u])
    return ds_m


def plot_NLCD_VLM_pie(exp):
    """ Plot a pie of each VLM type split up by landcover """
    ds_m = cmp_NLCD_VLM(exp)
    ds_m = ds_m.where(ds_m>0, np.nan) # get rid of water
    keys = list(ds_m.keys())
    dst  = op.join(PATH_RES, 'Figures', f'{Exp.reg}_VLM_Landcover_pie.pdf')
    pdf  = PdfPages(dst)

    ## plot a pie chart of the categories for each type
    for i, key in enumerate(keys):
        fig, axes = plt.subplots(figsize=(10, 10))
        df_vlm_lc  = ds_m[key].to_dataframe().reset_index()
        ser_vlm_lc = df_vlm_lc[key].dropna().astype(int)
        # skip the pie chart if there's basically no data in rate of VLM
        # this is < 0.5%
        npix       = ser_vlm_lc.shape[0]
        ## somewhat arbitrary cutoff
        if npix < 1000:
            print (f'Skipping {key} with only '\
                  f'({npix} pixels)')
            plt.close(fig)
            continue

        uniq, counts = np.unique(ser_vlm_lc, return_counts=True)
        percents     = 100*(counts/len(ser_vlm_lc))
        percents1    = percents[percents>1]
        uniq1        = uniq[percents>1]

        lbls         = [DCT_LULC_LBL[u] for u in uniq1]
        cols         = [DCT_LULC_COL[l] for l in lbls]
        sep          = [0.05]*len(lbls) # break pie

        axes.pie(percents1, explode=sep, labels=lbls,
                 colors=cols, autopct='%1.1f%%', #shadow=True,
                 startangle=45, textprops={'fontsize': PLFS})
        axes.axis('equal')
        axes.set_title(DCT_VUP[key][0])
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    pdf.close()
    print ('Wrote:', dst)
    return


def plot_NLCD_VLM_pie2(exp):
    """ Plot a pie of each landcover type split up by VLM """
    ds_m = cmp_NLCD_VLM(exp)
    ds_m = ds_m.where(ds_m>0, np.nan) # get rid of water
    df_m = ds_m.to_dataframe().reset_index().drop(columns='lat lon'.split())
    keys = list(ds_m.keys())
    classes   = sorted(list(DCT_LULC_LBL.keys()))
    colors    = plt.get_cmap('cmc.roma_r', len(DCT_VUP.keys())).colors
    dct_color = {v[0]: colors[i] for i, v in enumerate(DCT_VUP.values())}
    dst        = op.join(PATH_RES, 'Figures', f'{Exp.reg}_VLM_Landcover_pie2.pdf')
    pdf        = PdfPages(dst)

    for lc in classes:
        fig, axes = plt.subplots(figsize=(10, 10))
        df_m1 = df_m.dropna(how='all')[df_m == float(lc)].melt(var_name='lc').dropna()
        npix  = df_m1.shape[0]
        if npix < 1000:
            print (f'Skipping {DCT_LULC_LBL[lc]} with only '\
                  f'({npix} pixels)')
            plt.close(fig)
            continue

        # count sub/uplift
        uniq, counts = np.unique(df_m1['lc'], return_counts=True)
        percents     = 100*(counts/df_m1.shape[0])
        lbls         = [DCT_VUP[i][0] for i in uniq]
        sep          = [0.1]*len(lbls) # break pie
        axes.pie(percents, explode=sep, labels=lbls, #shadow=True
                colors=[dct_color[lbl] for lbl in lbls],
                 startangle=45, autopct='%1.1f%%', labeldistance=0.95,
                 textprops={'fontsize': PLFS})
        axes.axis('equal') # not allowed for 2 shared axis...
        axes.set_title(DCT_LULC_LBL[lc], fontsize=24)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()
    pdf.close()
    print ('Wrote:', dst)
    return


def plot_NLCD_VLM_hist(exp):
    ds_m = cmp_NLCD_VLM(exp)
    ds_m = ds_m.where(ds_m>0, np.nan) # get rid of water
    keys = list(ds_m.keys())
    ## plot a histogram
    df_m0 = ds_m.to_dataframe().reset_index().drop(columns='lat lon'.split())
    df_m1 = df_m0.where(~np.isnan(df_m0), -1)

    ## nan out where the percent type is too low
    for col in df_m0.columns:
        ser = df_m0[col].dropna()
        uniq, counts = np.unique(ser, return_counts=True)
        percents     = 100*(counts/len(ser))
        percents1    = percents[percents>1]
        uniq1        = uniq[percents>1]
        df_m1[col]   = df_m1[col].where(df_m1[col].isin(uniq1), -1)

    DCT_LULC_LBL[-1] = 'nan'
    df_m  = df_m1.applymap(lambda x: DCT_LULC_LBL[x]).astype('category').melt(var_name='key')
    df_m  = df_m[~(df_m['value'] == 'nan')]
    # cats  = df_m['value'].unique()

    # dct_col = {k: v for k, v in DCT_LULC_COL.items() if k in cats}

    df_m['key'] = df_m['key'].apply(lambda x: f'{DCT_VUP[x][0]} {DCT_VUP[x][1]}')

    fig, axes = plt.subplots(figsize=(12, 8))
    sns.countplot(x='key', data=df_m , hue='value', palette=DCT_LULC_COL)
    axes.grid(color='gray', linestyle = '--', linewidth=0.1)
    # axes.set_yscale('log')
    axes.legend(loc='upper right')
    axes.set_title(f'{Exp.reg} Landcover', fontsize=16)
    axes.set_xlabel('Vertical Land Motion Rate', fontsize=14)
    axes.set_ylabel('Pixel Counts', fontsize=14)

    fig.set_label(f'{Exp.reg}_VLM_Landcover_hist')
    return


def plot_NLCD(da_nlcd0):
    """ Plot the NLCD by itself (although combining classes) """
    da_nlcd  = combine_lulc(da_nlcd0)

    df_nlcd  = da_nlcd.to_dataframe().reset_index()
    ser_nlcd = df_nlcd['NLCD'].dropna().astype(int)
    ser_nlcd1 = ser_nlcd.apply(lambda x: DCT_LULC_LBL[x]).astype('category')

    cmap      = mpl.colors.ListedColormap(DCT_LULC_COL.values())
    basemap   = cimgt.Stamen('terrain-background')
    fig, axes = plt.subplots(figsize=(10, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    axes.add_image(basemap, 7, zorder=1)
    axes.add_feature(cfeature.OCEAN, color=WCOL, edgecolor=WCOL, zorder=5)
    axes.add_feature(cfeature.LAKES, color=WCOL, edgecolor=WCOL, zorder=5)
    im    = axes.pcolormesh(da_nlcd.lon, da_nlcd.lat, da_nlcd, shading='auto', zorder=10,
                        transform=ccrs.PlateCarree(), cmap=cmap)


    bbPlot.cartopy_cbar(im, xlabel='Landcover', pad=0.1)
    axes.set_title(f'{Exp.reg} Landcover')
    gl = axes.gridlines(draw_labels=True)
    bbPlot.fmt_gridlines(gl, bottom=True)
    fig.set_label(f'{Exp.reg}_Landcover_map')

    # plot histogram
    fig, axes = plt.subplots(figsize=(10, 10))
    sns.countplot(x=ser_nlcd1, palette=DCT_LULC_COL)
    axes.set_title(f'{Exp.reg} Landcover')
    axes.set_xlabel('')
    fig.set_label(f'{Exp.reg}_Landcover_hist')


    ### plot piechart
    # convert to percents
    uniq, counts = np.unique(ser_nlcd, return_counts=True)
    percents     = 100*(counts/len(ser_nlcd))
    lbls         = [DCT_LULC_LBL[u] for u in uniq]
    cols         = [DCT_LULC_COL[l] for l in lbls]
    sep          = [0.1]*len(lbls) # break pie

    fig, axes = plt.subplots(figsize=(10, 10))
    axes.pie(percents, explode=sep, labels=lbls,
             colors=cols, autopct='%1.1f%%', #shadow=True,
             startangle=45, textprops={'fontsize': PLFS})
    axes.axis('equal')
    fig.set_label(f'{Exp.reg}_Landcover_pie')


## UTILS -------------------------------------------------------------------- ##
def match_grid(da_rate, da_nlcd):
    """ Regrid the NLCD. (Maybe shifted by a pixel relative to InSAR) """
    import xesmf as xe
    ## this is identical (as it should be) to iterating over lat/lon in rate
    ## and taking nearest and its much faster
    regridder  = xe.Regridder(da_nlcd, da_rate, 'nearest_s2d')
    da_nlcd_re = regridder(da_nlcd)
    return da_nlcd_re


## NLCD ##
def load_NLCD(exp):
    """ Download by hand and unzip. After projecting, remove the 25GB file """
    da_insar  = load_rate_std(Exp)[0] # for more accurate geom
    S, N      = da_insar.lat.min(), da_insar.lat.max()
    W, E      = da_insar.lon.min(), da_insar.lon.max()
    path_root = op.join(os.getenv('dataroot'), 'GIS', 'nlcd_2019_land_cover_l48_20210604')
    path_raw  = op.join(path_root, 'nlcd_2019_land_cover_l48_20210604.img')
    path_crop = op.join(path_root, f'nlcd_2019_{exp.reg}.nc')

    if op.exists(path_crop):
        print ('Loading:', op.basename(path_crop))

    elif op.exists(path_raw):
        ## project to WGS84 and crop
        dst_srs = 'EPSG:4326'#'+proj=longlat +ellps=WGS84'
        WSEN    = [W, S, E, N]
        ds_proj = gdal.Warp(path_crop, path_raw, dstSRS=dst_srs, outputBounds=WSEN,
                                xRes=0.00027778, yRes=0.00027778, format='netcdf')
        ds_proj.FlushCache(); del ds_proj

    else:
        url_nlcd = 'https://s3-us-west-2.amazonaws.com/mrlc/nlcd_2019_land_cover_l48_20210604.zip'
        raise Exception('Download and extact the zip file by hand from:', url_nlcd)

    return xr.open_dataset(path_crop).rename(Band1='NLCD')


def combine_lulc(da_nlcd0):
    """ Combine the landcover classes into a few

    https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
    https://stackoverflow.com/questions/55949809/efficiently-replace-elements-in-array-based-on-dictionary-numpy-python
    """
    # may want to inspect
    uniq_counts = np.unique(da_nlcd0, return_counts=True)
    da_nlcd     = da_nlcd0.copy()

    k = np.array(list(DCT_LULC.keys()))
    v = np.array([int(vv[0]) for vv in list(DCT_LULC.values())])

    arr_map    = np.zeros(k.max()+1, dtype=v.dtype)
    arr_map[k] = v
    res        = arr_map[da_nlcd]

    da_nlcd.data = res

    return da_nlcd.where(da_nlcd>=0, np.nan)


## InSAR ##
def load_exp(fr_exp, mp_exp='Base', gps_exp='Default', neofs='15'):
    if gps_exp == 'Default':
        ref_sta = DCT_REG[fr_exp['root'].split('_')[0]][3][0]
    else:
        ref_sta = gps_exp

    exp     = ExpBase(fr_exp, mp_exp, gps_exp, neofs)
    return exp


def load_rate_std(exp, show=False):
    da_rate = xr.open_dataset(Exp.path_rate_msk_nc)['Band1']*1000
    da_unc  = xr.open_dataset(Exp.path_std_msk_nc)['Band1']*1000
    # da_resid  = xr.open_dataset(Exp.path_resid_msk_nc)['Band1']*1000
    if show:
        das[0].plot(norm=mpl.TwoSlopeNorm(0, -5, 5), cmap='cmc.roma_r')
        das[1].plot(norm=mpl.BoundaryNorm(np.linspace(0, 3, 4), 256), cmap='cmo.amp')

    return da_rate, da_unc


if __name__ == '__main__':
    # Exp     = load_exp(Charleston_SR, 'Base', 'SCHA', neofs=15)
    # Exp     = load_exp(NYC_SR, 'Base', 'NYBP', neofs=15)
    Exp     = load_exp(HR_SR, 'Base', 'LOY2', neofs=15)
    # Exp     = load_exp(Savannah_SR, 'Base', 'SAVA', neofs=15)

    # das   = load_rate_std(Exp, show=False)

    da_nlcd = load_NLCD(Exp)['NLCD']
    # plot_NLCD(da_nlcd)

    plot_NLCD_VLM_pie(Exp)
    plot_NLCD_VLM_pie2(Exp)
    # plot_NLCD_VLM_hist(Exp)


    bbPlot.savefigs(PATH_RES, True, True)
    plt.show()
