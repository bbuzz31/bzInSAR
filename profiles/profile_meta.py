#!/usr/bin/env python
"""
Metadata for profile transects.
** This is called by matlab **
"""
from scipy.io import savemat

from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from VLM.bzFRInGE.ARIABase import ARIABase

# must have at least two and different shapes for matlab to read right
## think none of this is necessary; just need a blank one to initialize region
# DCT_TRANSECTS = {'Charleston' : [
#         # WNES: -80.048750138, 32.8929168940, -79.732916931, 32.641250328
#         [[-80.0473513336293,32.8586780035838],
#         [-79.9651212479214,32.8584202290205],
#         # [-79.9380549187699,32.7831500565229], drop this one
#         [-79.9249084160391,32.7795412126361],
#         [-79.8421627812045,32.7573726001882]],

#         [[-79.8992706969624,32.892916893999995],
#         [-80.0052770268589,32.7957671017285],
#         [-79.9712749965147,32.7329633515634],
#         # [-79.9004707686216,32.7501643786787]],
#         [-79.89991, 32.75112]],

#         [[-79.7836637937923,32.892916893999995],
#         [-79.9992766685629,32.7401637815186]]],
#     # WNES -81.197083813, 32.152083637, -80.98375056500001, 32.006250361999996
#     'Savannah': [
#         [[-81.1165005262072,32.0067339047308],
#         [-81.0831712834927,32.095904247081]],

#         [[-81.173803434734,32.024860334979], # move the first point left to make cleaner
#         [-81.103003434734,32.024860334979], # add a second point to shift midpoint
#         # [[-81.184803434734,32.024860334979],
#         [-81.053935105673,32.0236908878662],
#         [-81.0077419447178,32.0497110861258]],

#         [[-81.1007129901846,32.0073186282872],
#         [-81.1776041378505,32.1505758996039]]],

#     'HR': [
#             # 1 = A, 2 = B, 3 = C
#             [[-75.9410, 36.7870],[-76.2587, 36.7603],[-76.3000, 36.8040],
#             [-76.3485, 36.8816], [-76.36, 36.9044], [-76.4043, 37.0586],
#             [-76.5546, 37.1605], [-76.57, 37.17]],

#             [[-75.9896, 36.9374],[-76.2393, 36.7628],[-76.4091, 36.6633]],

#             # first two are masked in FRInGE; near the VAHP seems to work
#             [[ -76.401, 36.84], [-76.4043, 37.0586], #[-76.575, 36.851],[-76.5546, 36.9689],
#             [-76.3897, 37.0659], [-76.356, 37.085], [-76.45, 37.155]]],
#     }

DCT_TRANSECTS = {k: [] for k in 'HR Charleston Savannah NYC DC Houston'.split()}

## stations to mark as reference
DCT_REF = {'HR': 'LOY2 VAHP LOYZ'.split(),
            'Charleston': ['SCHA'],
            'Savannah': ['SAVA'],
            'NYC'     : ['NYBK'],
            'DC'     : ['USN8'],
            'Houston': ['NASA']
           }

## exlcude these stations
DCT_EX = {'HR': [''], #'WLP2 LNG4'.split(),
           }

## x lims, y lims, radius; x lims only for third row? its hardcoded in ..._REG
DCT_LIMS = {'HR': [(0, 100), (-14, 14), 500],
            'Charleston': [(0, 30), (-10, 5), 250],
            'Savannah': [(0, 18), (-10, 5), 250],
            'NYC'     : [(0, 80), (-5, 4), 100],
            'DC'     : [(0, 28), (-4, 4), 100],
            'Houston'     : [(0, 70), (-15, 1.5), 50],
           }


def line_to_points(src, dst=None):
    """ Convert a linestring to points (for transects)

    convert a google earth file to geojson:
        ogr2ogr transect.GeoJSON transect.kmz
    """
    import geopandas as gpd
    from BZ.bbPlot import plot_bbox
    from shapely.geometry import Polygon, LineString
    from subprocess import call
    if src.suffix == '.kmz':
        f_new = f'{op.splitext(src)[0]}.GeoJSON'
        dst = op.join(dst, op.basename(f_new)) if dst else f_new
        cmd = f'ogr2ogr {dst} {src}'
        call(cmd.split())
        src = dst
        print ('Wrote:', dst)
        if 'Houston' in str(src):
            print ('Reversing line string')
            gdf = gpd.read_file(src)
            gdf['geometry'] = gdf['geometry'].apply(
                lambda geom: LineString(list(geom.coords)[::-1]) if geom.geom_type == 'LineString' else geom
            )
            gdf.to_file(src)

    gdf = gpd.read_file(src)
    gdf1 = gdf.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    df_pts = pd.DataFrame(gdf1.to_numpy()[0])
    df_pts.columns = 'lon lat ?'.split()
    # print (df_pts)

    # for copying into the dict in profiles_meta.py
    lst_pts = [(np.round(x, 6), np.round(y, 6)) for x, y in zip(df_pts.lon, df_pts.lat)]
    # print (lst_pts)
    # breakpoint()
    return lst_pts


def createParser():
    parser = argparse.ArgumentParser(description='Get metadata for matlab profiles.',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('transects', help='Names of the transects')
    parser.add_argument('dct_exp_n', type=str, help='Name of experiment dictionary')
    parser.add_argument('mp_exp', default='Base', type=str, help='mintpy experiment name (default: Base)')
    parser.add_argument('ref_sta', default=None, type=str, help='reference GPS Station (defaults to experiment)')
    parser.add_argument('neofs', default='', type=str, help='number of eofs in reconstruction to use (default: '')')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 3:
        parser.print_help()
        os.sys.exit(0)

    inps = parser.parse_args(args=iargs)
    return inps


def main(inps):
    all_vars = globals()
    dct_exp  = all_vars[inps.dct_exp_n]

    if 'ARIA' in dct_exp['root']:
        Obj = ARIABase(dct_exp, inps.mp_exp)
    else:
        Obj = ExpBase(dct_exp, inps.mp_exp, ref_sta=inps.ref_sta, neofs=inps.neofs)

    stitched = True if 'stitch' in inps.mp_exp.lower() else False

    reg      = Obj.reg

    ## try to add the NYC transects
    # for i, transect in enumerate('NJ1 Manhattan1 Brooklyn1'.split()):
    # for i, transect in enumerate('NJ_NYOB Manhattan Brooklyn'.split()):
    path_transects = Path(Obj.path_wd) / 'Transects'
    transects = inps.transects.split(',') if isinstance(inps.transects, str) else inps.transects
    for i, transect in enumerate(transects):
        path = path_transects / f'{transect}.GeoJSON'
        if not path.exists():
            path = path.with_suffix('.kmz')
        assert path.exists(), f'Could not get: {path}'
        DCT_TRANSECTS[reg].append(line_to_points(path))

    dct_data = {}
    dct_data['REG']      = reg
    dct_data['COORDS']   = DCT_TRANSECTS[reg]

    dct_data['REF_STAS'] = DCT_REF[reg]
    dct_data['EX_STAS']  = DCT_EX.get(reg, [''])

    xlims, ylims, radius = DCT_LIMS[reg]
    dct_data['XLIMS']    = xlims
    dct_data['YLIMS']    = ylims
    dct_data['RADIUS']   = float(radius)
    dct_data['PATH_GPS'] = Obj.path_gps
    dct_data['PATH_RATE']= Obj.path_rate_msk_nc_stitch if stitched else Obj.path_rate_msk_nc
    dct_data['PATH_STD'] = Obj.path_std_msk_nc if stitched else Obj.path_std_msk_nc
    # path_figs            = op.join(op.dirname(dct_data['PATH_RATE']), 'Figures')
    # path_figs            = Path.home() / 'Desktop', 'VLM', f'{reg}_2024'
    dct_data['PATH_FIGS']= str(PATH_RES)

    ## save variable in a matlab file in same directory so matlab doesnt need to find it
    dst = Path.cwd() / 'profile_parms.mat'

    savemat(dst, dct_data, oned_as='column')
    print (f'Wrote {reg} parameters to:', dst)
    return


if __name__ == '__main__':
    inps = cmdLineParse()
    # inps = argparse.Namespace(dct_exp_n='HR_SR', mp_exp='Base', neofs=15)
    main(inps)
