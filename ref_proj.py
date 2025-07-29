#!/usr/bin/env python3
import shutil
from collections import OrderedDict
import h5py
from shapely.geometry import Point, Polygon

# fails depending on directory
from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase

from mintpy.utils import readfile, writefile, utils0, utils as ut

from BZ.bbLogger import logger
gdal.UseExceptions()

## ------------------------------------------------------------------ Parsers ##
def createParser():
    parser = argparse.ArgumentParser(description='Program to project referenced separately velocities\n'\
        'References HR to GPS LOY2 in the south and VAHP in the north.\n\t',
             epilog='Examples of use: ProjHR.py ./MintPy_2alks_5rlks_33_15',
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path_mp', type=str, help="Path to MintPy dir")
    parser.add_argument('method', type=str, choices=['simple', 'stitch', 'stitch+'],
        help="Reference to to stations with 'simple' or stich a few with 'stitch(+)'")
    parser.add_argument('-m', '--mask', type=str, default=None, help=
        "Optional path to external mask. Requires that geo_waterMask.nc exists at 'path_mp'")
    parser.add_argument('-c', '--corr', type=str, default='',
                    help="Optionally use unwrapping error corrections")
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    inps   = parser.parse_args(args=iargs)
    return inps


## --------------------------------------------------------------- Geom Utils ##
def make_circle_rast(lalo, radius, path_ds, dst=None):
    """
    Create a circle of radius 'radius' (m) around 'lalo'
    The output will match the size of the input 'rast' but have values of 1/9999
        1 == in circle, 9999 not within circle
    The geometry will be saved to 'dst'.shp and the resulting raster to 'dst'
        if dst is None the data will be created and erased
    Need to pass a path to a geocoded dataframe, such as the dem
    """
    from osgeo import gdal, ogr, osr
    dst = op.join(os.getcwd(), 'temp') if dst is None else dst
    os.remove(f'{dst}.shp') if op.exists(f'{dst}.shp') else '' # GeoJSON cant overwrite


    ds2copy = gdal.Open(path_ds, gdal.GA_ReadOnly) if \
                                    isinstance(path_ds, str) else path_ds
    nodata  = ds2copy.GetRasterBand(1).GetNoDataValue()

    ds_out = gdal.GetDriverByName('MEM').CreateCopy('', ds2copy)
    ds_out.GetRasterBand(1).WriteArray(np.ones(ds2copy.ReadAsArray().shape))

    ds = ogr.GetDriverByName('GeoJSON').CreateDataSource(f'{dst}.shp')

    # create layer
    layer = ds.CreateLayer('', None, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger)) # Add 1 attribute
    feat = ogr.Feature(layer.GetLayerDefn()) # Create a new feature (attribute and geometry)
    feat.SetField('id', 0)

    # Make a geometry, from input Shapely object
    poly = circle(lalo, radius)
    geom = ogr.CreateGeometryFromWkb(poly.wkb)
    feat.SetGeometry(geom)
    layer.CreateFeature(feat)
    ds = layer = feat = geom = None
    ds_new  = gdal.Warp(dst, ds_out, width=ds_out.RasterXSize, height=ds_out.RasterYSize,
                        format='GTiff', cutlineDSName=f'{dst}.shp', dstNodata=nodata)
    # dict_ds[rad] = ds_new
    ds_new.FlushCache()
    arr = ds_new.ReadAsArray()
    if op.basename(dst) == 'temp': [os.remove(i) for i in [dst, f'{dst}.shp']]
    return ds_new


def circle(lalo, radius, epsg=4087):
    """ Make a circle with radius (meters) around lalo """
    ## ALL UTM FUNCTIONALITY IS DEPRECATED PENDING HR INSAR TEST
    # import utm
    if isinstance(lalo, Point):
        lalo = [lalo.y, lalo.x]
    x, y  = transform_coords(4326, epsg, lalo[1], lalo[0])
    # ext_utm = utm.from_latlon(*lalo) # returns x, y, zone number
    # lat, lon = utm.to_latlon(ext_utm[0], ext_utm[1], ext_utm[2], ext_utm[3])
    # pt_utm  = Point(ext_utm[:2])
    buf = Point([x, y]).buffer(radius, cap_style=1)
    buf_ll = []
    for x, y in zip(buf.exterior.xy[0], buf.exterior.xy[1]):
        # lat, lon = utm.to_latlon(x, y, ext_utm[2], ext_utm[3])
        lon, lat = transform_coords(epsg, 4326, x, y)
        buf_ll.append([lon, lat])

    poly = Polygon(buf_ll)
    return poly


def transform_coords(proj1, proj2, x, y):
    """
    Transform coordinates from proj1 to proj2 (EPSG num).
    e.g. x, y = transform_coords('4326','4087', lon, lat)
    """
    from pyproj import Transformer

    transformer = Transformer.from_crs(int(proj1), int(proj2), always_xy=True)

    return transformer.transform(x, y)


## --------------------------------------------------------------- Read Utils ##
def read_inc(path_mp, geo='geo_', xy=None):
    assert geo in ['geo_', ''], f'Incorrect geo prefix:{geo} passed to read_inc'
    if geo:
        path_pre = ''
    else:
        path_pre = 'inputs'

    path_geom = op.join(path_mp, path_pre, f'{geo}geometryGeo.h5')


    if not op.exists(path_geom):
        path_geom = path_geom.replace(f'{geo}geometryGeo.h5', f'{geo}geometryRadar.h5')

    with h5py.File(path_geom, 'r') as h5:
        inc = h5['incidenceAngle'][:]
        inc_cos  = np.cos(np.deg2rad(inc))
        inc_cos  = inc_cox[xy[1], xy[0]] if xy is not None else inc_cos
    return inc_cos


def read_vel(path_vel, xy=None):
    with h5py.File(path_vel, 'r') as h5:
        rate  = h5['velocity'][:]
        std   = h5['velocityStd'][:]
        res   = h5['residue'][:] if 'residue' in list(h5.keys()) else np.zeros_like(rate)
        units = h5.attrs['UNIT'][0]
        # because GPS are in mm/yr
        if units == 'mm':
            rate /= 1000
            std  /= 1000
            res  /= 1000

    if xy is not None:
        rate = rate[xy[1], xy[0]]
        std  = std [xy[1], xy[0]]
        res  = res [xy[1], xy[0]]

    return rate, std, res


## --------------------------------------------------------------- Mask Utils ##
def apply_emask(path_mask, arr):
    """ Apply an external (h5) mask """
    if path_mask is None:
       emask  = np.ones_like(arr)
    else:
        with h5py.File(path_mask, 'r') as h5:
            try:
                emask = h5['mask'][:]
            except:
                emask = h5['waterMask'][:]
        emask = np.where(np.isclose(emask, 0), np.nan, 1)
    return arr*emask


def make_crop_mask(path_shp, path_wmask, show=False):
    dstb  = f'{op.splitext(op.basename(path_shp))[0]}.vrt'
    dst   = op.join(op.dirname(path_wmask), dstb)
    w, h  = gdal.Info(path_wmask, format='json')['size']
    # print (path_wmask)
    # print (gdal.Info(path_wmask, format='json'))
    # return
    ds    = gdal.Warp(dst, path_wmask, format='VRT', width=w, height=h,
                    cutlineDSName=path_shp, dstNodata=0)

    if show:
        arr   = ds.ReadAsArray()
        plt.imshow(arr, interpolation='nearest', cmap='binary')
        plt.show()
    del ds
    print ('Made water + crop mask', op.basename(dst))# to:', dst)
    return dst


class RefProj(object):
    """ Path emask is usually temporal coherence """
    def __init__(self, path_mp, region, path_wmask, path_emask=None, corr='SR', show_cmds=False):
        self.path_mp = path_mp
        self.reg     = region
        self.corr    = corr # potential unwrapping corrections
        self.root    = os.getenv('dataroot')
        self.path_wmask = path_wmask # geo_waterMask.nc
        self.path_emask = path_emask # additional optional mask (hdf5)
        self.path_gps   = op.join(self.root, 'VLM', 'GNSS', 'UNR',
                                        f'Midas_GPS-{self.reg}_cln.csv')
        self.gps        = pd.read_csv(self.path_gps, index_col=0)
        self.show_cmds  = show_cmds
        self.geo_pre    = '' if corr == 'ARIA' else 'geo_'


    def apply_mask(self, path_mask):
        """ Get LOS vel and apply a regional mask; for LOY2 in stitching """
        path_vel   = op.join(self.path_mp, f'{self.geo_pre}velocity.h5')
        assert op.exists(path_vel), f'Create {path_vel} by geocoding vel (or running MintPy)'
        rate, std, res = read_vel(path_vel)
        inc_cos    = read_inc(self.path_mp, self.geo_pre)
        rate = apply_emask(self.path_emask, rate)
        std  = apply_emask(self.path_emask, std)
        res  = apply_emask(self.path_emask, res)

        ins_uvel = rate / inc_cos

        ins_ustd = std / inc_cos

        mask = gdal.Open(path_mask, gdal.GA_ReadOnly).ReadAsArray()
        mask = np.where(mask == 0, np.nan, 1)

        return mask*ins_uvel, mask*ins_ustd, mask*res


    ## should really be taking a circle around the gps station
    def proj_one(self, sta='VAHP', path_vel=None, path_mask=None):
        """ Tie to a single GPS OR the average of all/ a list

        pass path_mp for stitching (i.e., a different referenced run)
        path_mask is typically the north/south mask
        """
        # path_mp    = self.path_mp if path_mp is None else path_mp
        # ## if the station you pass isn't in the mintpy file directory, switch
        # if not sta in path_mp:
        #     if op.basename(path_mp) == 'geo':
        #         sta_wrong = op.basename(op.dirname(path_mp)).split('_')[0]
        #     else:
        #         sta_wrong = op.basename(path_mp).split('_')[0]
        #
        #     path_mp   = self.path_mp.replace(sta_wrong, sta)


        path_vel   = op.join(self.path_mp, f'{self.geo_pre}velocity.h5') \
                                            if path_vel is None else path_vel
        print ('Using:', path_vel)
        assert op.exists(path_vel), f'Create {path_vel} by geocoding vel'
        # assert sta in path_vel, f'Incorrect mp experiment for this station ({sta})'

        rate, std, res = read_vel(path_vel)

        inc_cos    = read_inc(self.path_mp, self.geo_pre)
        rate = apply_emask(self.path_emask, rate)
        std  = apply_emask(self.path_emask, std)
        res  = apply_emask(self.path_emask, res)

        # regional mask
        if path_mask is None:
            mask = np.ones_like(rate)
        else:
            mask = gdal.Open(path_mask, gdal.GA_ReadOnly).ReadAsArray()
            mask = np.where(mask == 0, np.nan, 1)

        rate*= mask
        std *= mask
        res *= mask

        if isinstance(sta, str):
            ref_gps = self.gps if sta == 'all' else self.gps[self.gps.sta==sta]

        elif isinstance(sta, list):
            ref_gps = self.gps[self.gps.sta.isin(sta)]

        ref_gps_uvel = ref_gps.u_vel.mean()
        ref_gps_ustd = ref_gps.u_sig.mean()

        print (f'GPS Offset: {ref_gps_uvel*1000:.2f} +/- {ref_gps_ustd*1000:.2f} mm/yr')

        ins_gps_uvel = (rate / inc_cos) + ref_gps_uvel

        ins_gps_vstd = np.sqrt( ((std / inc_cos)**2) + (ref_gps_ustd**2) )

        ## write array as mintpy readable object
        if self.geo_pre:
            base = '_'.join(op.basename(op.dirname(self.path_mp)).split('_')[1:])
        else:
            base = '_'.join(op.basename(self.path_mp).split('_')[1:])

        path_mp = self.path_mp if not self.corr == 'ARIA' else op.join(self.path_mp, 'Vup')
        os.makedirs(path_mp, exist_ok=True)

        dst  = op.join(path_mp, f'{self.geo_pre}Vup_{base}_{self.corr}.h5')
        shutil.copy(path_vel, dst)
        ## dont keep GPS name, just the Base/atm etc...
        with h5py.File(dst, 'r+') as h5_new:
            data    = h5_new['velocity']
            data[:] = mask*ins_gps_uvel

            data    = h5_new['velocityStd']
            data[:] = mask * ins_gps_vstd

            if 'residue' in h5_new.keys():
                data    = h5_new['residue']
                data[:] = res
            else:
                h5_new.create_dataset('residue', data=res)

            h5_new.attrs['FILE_PATH'] = dst
        print (f'\nWrote: {dst} (projected {sta}) at {datetime.now()}')

        return dst


    def proj_two(self, path_vel0, path_vel1, path_vel2=None):
        """ Tie crops regions (made manually) to diff GPS stations """
        paths_vels = [path_vel0, path_vel1]

        vels, stds, ress = [], [], []
        for i, path_vel in enumerate(paths_vels):
            sta = op.basename(op.dirname(op.dirname(path_vel))).split('_')[0]
            v   = CROP_MAP[sta]

            ## crop the watermask to the right size
            path_crop   = op.join(self.root, 'Crops', self.reg, f'{v}.shp')
            path_cmask  = make_crop_mask(path_crop, self.path_wmask)

            ## mask, project it to vertical and propagate unc
            path_result = self.proj_one(sta, path_vel, path_mask=path_cmask)
            with h5py.File(path_result, 'r') as h5:
                vels.append(h5['velocity'][:])
                stds.append(h5['velocityStd'][:])
                ress.append(h5['residue'][:])
            # breakpoint()

        # add the offset for HR2020
        # if 'GRL' in self.path_mp:
            # stds[1]  = np.sqrt(stds[1]**2 + (0.51**2))

        vel = np.where(np.isnan(vels[0]), vels[1], vels[0])
        std = np.where(np.isnan(stds[0]), stds[1], stds[0])
        res = np.where(np.isnan(ress[0]), ress[1], ress[0])

        if path_vel2 is not None:
            sta = op.basename(op.dirname(op.dirname(path_vel2))).split('_')[0]
            v   = CROP_MAP[sta]

            ## crop the watermask to the right size
            path_crop   = op.join(self.root, 'Crops', self.reg, f'{v}.shp')
            path_cmask  = make_crop_mask(path_crop, self.path_wmask)

            ## mask, project it to vertical and propagate unc
            path_result = self.proj_one(sta, path_vel2, path_mask=path_cmask)
            with h5py.File(path_result, 'r') as h5:
                vel_new = h5['velocity'][:]
                std_new = h5['velocityStd'][:]
                res_new = h5['residue'][:]
            print ('Masking out LOYZ region')
            vel     = np.where(np.isnan(vel_new), vel, np.nan)
            std     = np.where(np.isnan(std_new), std, np.nan)
            res     = np.where(np.isnan(res_new), res, np.nan)


        ## make it a mintpy object; could possibly hit with watermask
        src  = op.join(self.path_mp, f'{self.geo_pre}velocity.h5')
        ## dont keep GPS name, just the Base/atm etc...
        if self.geo_pre:
            base = '_'.join(op.basename(op.dirname(self.path_mp)).split('_')[1:])
        else:
            base = '_'.join(op.basename(self.path_mp).split('_')[1:])

        path_mp = self.path_mp if not self.corr == 'ARIA' else op.join(self.path_mp, 'Vup')
        dst     = op.join(path_mp, f'{self.geo_pre}Vup_{base}_{self.corr}.h5')
        print (dst)

        with h5py.File(dst, 'r+') as h5_new:
            data    = h5_new['velocity']
            data[:] = vel

            data    = h5_new['velocityStd']
            data[:] = std

            if 'residue' in h5_new.keys():
                data    = h5_new['residue']
                data[:] = res
            else:
                h5_new.create_dataset('residue', data=res)

            h5_new.attrs['FILE_PATH'] = dst

        with h5py.File(dst, 'r') as h51:
            print (h51['residue'][500, 700])

        print (f'\nWrote Vup to: {dst}\n')
        if self.show_cmds:
            print (f'view.py {dst} velocity -u mm -v -10 10 -c roma')
            print (f'view.py {dst} velocityStd -u mm -v 0 5 -c afmhot_r')
        return


    def proj_HR(self, radius, ref_sta='LOY2', ref_staN='VAHP', cut_monroe=True):
        """ Project the south half based on several GPS; radius in m around sta

        Reference the HR Southeast/west separately
        """
        assert ref_sta in self.path_mp, 'You arent using LOY2 as ref sta?'
        stas_south = 'LOY2 LS03 SPVA'.split()

        ## project southern reference sta to vertical and propagate unc
        path_crop   = op.join(self.root, 'Crops', 'HR', f'HR_South.shp')
        path_cmask  = make_crop_mask(path_crop, self.path_wmask)
        # nan velocity pixels in the path_cmask (North)
        rate0, unc0, res0 = self.apply_mask(path_cmask)

        shp0 = rate0.shape
        rate = rate0.reshape(-1)
        unc  = unc0.reshape(-1)

        print ('Radius around GPS:', radius)

        ## ~~~~~~~~~~~~~~~ CALCULATE OFFSET ~~~~~~~~~~~~~~~ ##
        offsets, weights = [], []
        for j, tie_sta in enumerate(stas_south):
            gps  = self.gps[self.gps.sta == tie_sta]

            lalo = gps[['lat', 'lon']].values

            ## get the insar pixel locations around GPS
            ds_circle   = make_circle_rast(*lalo, radius, self.path_wmask)
            nodata      = ds_circle.GetRasterBand(1).GetNoDataValue()
            nodata      = 0 if nodata is None else nodata # might break eventually
            arr         = ds_circle.ReadAsArray().reshape(-1)
            idx         = np.where(~np.isclose(arr, nodata))[0]

            ## get insar rates at the gps
            rate_at_gps = np.ma.masked_invalid(rate[idx]).compressed()
            unc_at_gps  = np.ma.masked_invalid(unc[idx]).compressed()

            if len(rate_at_gps) == 0 or len(unc_at_gps) == 0:
                print (f'No nonnan pixels around {tie_sta}, skipping...')

            else:
                print (f'{len(rate_at_gps)} pixels around {tie_sta}')
                ## note that at the ref sta, because the insar pixels are 0,
                ## this will match the GPS rate
                offsets.append(gps.u_vel.item() - np.nanmean(rate_at_gps))
                weights.append(gps.u_sig.item()**2)


        ## ~~~~~ PERFORM WEIGHTED MEAN OF OFFSETS USING GNSS ONLY ~~~~~ ##
        wgts         = np.sqrt(weights)
        _,     cov0  = np.polyfit(np.ones(3), offsets, 0, w=1/wgts, cov=True)
        ## seems higher for HR, and seems to be more correct based on docs
        offset, cov1 = np.polyfit(np.ones(j+1), offsets, 0, w=1/wgts, cov='unscaled')

        offset_unc   = np.sqrt(np.max([cov0, cov1])) # use the bigger one

        print(f'Stitch Offset: {offset.item():.2e} +/- {offset_unc.item():.2e}')
        # print(f'Mean GPS variance to add: {wgts.mean()**2:.2f}')

        ## add the GPS uncertainty along with the offset?
        gps_ref      = self.gps[self.gps.sta==ref_sta]
        gps_ref_usig = gps_ref.u_sig.item()

        rate_S = rate.reshape(shp0) + offset
        unc_S  = np.sqrt(unc.reshape(shp0)**2 + offset_unc**2 + gps_ref_usig**2)

        ## combine with Northern portion
        path_crop   = op.join(self.root, 'Crops', 'HR', f'HR_North.shp')
        path_cmask  = make_crop_mask(path_crop, self.path_wmask)
        ## project it to vertical and propagate unc
        path_VAHP   = op.join(self.path_mp.replace('LOY2', 'VAHP'), f'{self.geo_pre}velocity.h5')
        path_resN   = self.proj_one(ref_staN, path_VAHP, path_cmask)

        ## this has nans at water and 0 elsewhere
        with h5py.File(path_resN, 'r') as h5:
            rate_N = h5['velocity'][:]
            unc_N  = h5['velocityStd'][:]
            res_N  = h5['residue'][:]

        # where the South values are nan, put the North values
        rate_SN = np.where(~np.isnan(rate_S), rate_S, rate_N)
        unc_SN  = np.where(~np.isnan(unc_S), unc_S, unc_N)
        res_SN  = np.where(~np.isnan(res0), res0, res_N)

        ## combine with HR SouthWest
        path_crop   = op.join(self.root, 'Crops', 'HR', f'HR_Southwest.shp')
        path_cmask  = make_crop_mask(path_crop, self.path_wmask)
        ## project it to vertical and propagate unc
        path_LOYZ   = op.join(self.path_mp.replace('LOY2', 'LOYZ'), f'{self.geo_pre}velocity.h5')
        path_resW   = self.proj_one('LOYZ', path_LOYZ, path_cmask)

        ## water is nans
        with h5py.File(path_resW, 'r') as h5:
            rate_SW = h5['velocity'][:]
            unc_SW  = h5['velocityStd'][:]
            res_SW  = h5['residue'][:]

        # SW is nan outside of crop; keep nonnans else take South/north
        rate = np.where(np.isnan(rate_SW), rate_SN, rate_SW)
        unc  = np.where(np.isnan(unc_SW), unc_SN, unc_SW)
        res  = np.where(np.isnan(res_SW), res_SN, res_SW)

        ## optionally mask out Fort Monroe
        path_crop   = op.join(self.root, 'Crops', 'HR', f'Fort_Monroe.shp')
        path_cmask  = make_crop_mask(path_crop, self.path_wmask)
        cmask = gdal.Open(path_cmask, gdal.GA_ReadOnly).ReadAsArray()
        cmask = np.where(np.isclose(cmask, 1), np.nan, 0)
        if cut_monroe:
            rate += cmask
            unc  += cmask
            res  += cmask

        # mmm(rate*1000) # diagnostic

        src  = op.join(self.path_mp, f'{self.geo_pre}velocity.h5')
        ## dont keep GPS name, just the Base/atm etc...
        base = '_'.join(op.basename(op.dirname(self.path_mp)).split('_')[1:])
        dst  = op.join(self.path_mp, f'{self.geo_pre}Vup_{base}_{self.corr}.h5')
        shutil.copy(src, dst)

        with h5py.File(dst, 'r+') as h5_new:
            data    = h5_new['velocity']
            data[:] = rate

            data    = h5_new['velocityStd']
            data[:] = unc

            data    = h5_new['residue']
            data[:] = res
            h5_new.attrs['FILE_PATH'] = op.join(dst)

        print (f'\nWrote Vup to:', dst)
        if self.show_cmds:
            print (f'\nview.py {dst} velocity -u mm -v -5 5 -c roma')
            print (f'\nview.py {dst} velocityStd -u mm -v 0 5 -c afmhot_r')
        return dst


def proj_all(Exp, path_vel_los, radius=50, exclude=[]):
    """ Use all stations to compute a single offset """
    from BZ import bbGIS
    import geopandas as gpd

    # get the GPS
    df_gps  = prep_gps(Exp.path_gps, Exp.reg)
    ser_ref = df_gps[df_gps.sta == Exp.ref_sta]
    # df_gps.drop(ser_ref.index, inplace=True)

    df_gps = df_gps[~df_gps.sta.isin(exclude)]
    df_gps = bbGIS.in_dfll(df_gps, SNWE=Exp.SNWE).set_index('sta')

    # get the insar rates
    epsg    = 'epsg:4326'
    da_rate = h5_to_xr(path_vel_los, 'velocity')
    da_rate.rio.write_crs(epsg, inplace=True)

    da_unc  = h5_to_xr(path_vel_los, 'velocityStd')
    da_unc.rio.write_crs(epsg, inplace=True)

    da_res  = h5_to_xr(path_vel_los, 'residue')
    da_res.rio.write_crs(epsg, inplace=True)

    inc = readfile.read(Exp.path_geom_mp_geo, datasetName='incidenceAngle')[0]

    ## convert to vertical
    da_rate *= np.flipud((1/np.cos(inc)))
    da_unc  *= np.flipud((1/np.cos(inc)))


    offsets = []
    weights = []
    for sta, row in df_gps.iterrows():
        lat, lon = row['lat'], row['lon']
        gser = gpd.GeoSeries(gpd.points_from_xy([lon], [lat]), crs=epsg)
        poly = gser.to_crs(4087).buffer(radius, cap_style=1).to_crs(epsg).geometry
        rate_at_gps = da_rate.rio.clip(poly, epsg)
        # unc_at_gps  = da_unc.rio.clip(poly, epsg) # not used
        n_nan = np.isnan(np.where(np.isclose(rate_at_gps, 0), np.nan, rate_at_gps)).sum()

        ## wont be 0 due to gdal interpolation
        if sta == ser_ref.sta.item():
            offsets.append(row['u_vel'].item())
            weights.append(row['u_sig'].item()**2)

        if rate_at_gps.size == n_nan:
            print (f'No nonnan/0 pixels around {sta}, skipping...')

        else:
            print (f'{rate_at_gps.size-n_nan} pixels around {sta}')
            ## note that at the ref sta, because the insar pixels are 0,
            ## this will match the GPS rate
            offsets.append(row['u_vel'].item() - np.nanmean(rate_at_gps))
            weights.append(row['u_sig'].item()**2)

    wgts     = np.sqrt(weights)
    # weighted average; take higher
    _, cov0  = np.polyfit(np.ones(len(offsets)), offsets, 0, w=1/wgts, cov=True)
    offset, cov1 = np.polyfit(np.ones(len(offsets)), offsets, 0, w=1/wgts, cov='unscaled')

    offset_unc   = np.sqrt(np.max([cov0, cov1])) # use the bigger one

    print(f'Stitch Offset: {offset.item():.2e} +/- {offset_unc.item():.2e}')
    # print(f'Mean GPS variance to add: {wgts.mean()**2:.2f}')

    ## add the GPS uncertainty along with the offset?
    # no; i already put it into the weights
    # gps_ref_usig = ser_ref.u_sig.item()

    da_rate += offset
    da_unc  += np.sqrt(da_unc**2 + offset_unc**2)

    rate = np.flipud(da_rate.data)
    unc  = np.flipud(da_unc.data)

    shutil.copy(path_vel_los, Exp.path_vup_geo)
    writefile.write_hdf5_block(Exp.path_vup_geo, rate, 'velocity')
    writefile.write_hdf5_block(Exp.path_vup_geo, unc, 'velocityStd')
    return


def proj_NYC(ExpA, ExpB, npix=0):
    """ Project to East and West and then stitch together """
    from BZ import bbGIS
    import geopandas as gpd
    epsg     = 'epsg:4326'

    lst_rates, lst_uncs, lst_ress = [], [], []
    for Expi in [ExpA, ExpB]:
        path_vup = proj_one(Expi, npix=npix)
        da_rate = h5_to_xr(path_vup, 'velocity')
        da_rate.rio.write_crs(epsg, inplace=True)

        da_unc  = h5_to_xr(path_vup, 'velocityStd')
        da_unc.rio.write_crs(epsg, inplace=True)

        # da_res  = h5_to_xr(path_vup, 'residue')*1000
        # da_res.rio.write_crs(epsg, inplace=True)

        ## clip the rates to polygon
        path_crop    = op.join(Expi.path_crops, f'{CROP_MAP[Expi.ref_sta]}.GeoJSON')
        gdf_crop     = gpd.read_file(path_crop)
        da_rate_crop = da_rate.rio.clip(gdf_crop.geometry, epsg, drop=False)
        da_unc_crop  = da_unc.rio.clip(gdf_crop.geometry, epsg, drop=False)
        # da_res_crop  = da_res.rio.clip(gdf_crop.geometry, epsg, drop=False)

        lst_rates.append(da_rate_crop)
        lst_uncs.append(da_unc_crop)

        logger.info (f'Mean Velocity in {Expi.ref_sta} {CROP_MAP[Expi.ref_sta]} '\
               f'{da_rate_crop.mean().item()*1000:.2f} mm/yr')

        # lst_ress.append(da_res_crop)


    rate_stitch = np.nanmean(lst_rates, 0)
    unc_stitch  = np.nanmean(lst_uncs, 0)

    # overlap = np.sum(lst_rates, 0)
    ix = np.logical_and(~np.isnan(lst_rates[0]), ~np.isnan(lst_rates[1]))
    overlap = np.where(ix, 1, np.nan)

    # do a weighted mean in the overlap region
    # technically should propagate variance of GPS but same for NYC (0.8 for NJHT, NYBK)
    if not np.isnan(overlap).all():
        log.warning('There is overlap between the polygons, averaging overlap.')
        unc_ovl = np.stack([lst_uncs[0].data*overlap, lst_uncs[1].data*overlap])
        wgts     = 1/(unc_ovl**2)
        rate_ovl = np.stack([lst_rates[0].data*overlap, lst_rates[1].data*overlap])
        rate_ovl = np.average(rate_ovl, 0, wgts)
        unc_ovl  = np.nanmean(unc_ovl, 0)

        rate_stitch= np.where(np.isnan(rate_ovl), rate_stitch, rate_ovl)
        unc_stitch= np.where(np.isnan(unc_ovl), unc_stitch, unc_ovl)


    ## just make both; they'll be the same
    ref_stas = [readfile.read_attribute(Expi.path_vup_geo
                                ).get('ref_stas') for Expi in [ExpA, ExpB]]
    for Expi in [ExpA, ExpB]:
        shutil.copy(Expi.path_vup_geo, Expi.path_vup_geo_stitch)
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, np.flipud(rate_stitch), 'velocity')
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, np.flipud(unc_stitch), 'velocityStd')
        # writefile.write_hdf5_block(Expi.path_vup_geo_stitch, res_stitch, 'residue')
        with h5py.File(Expi.path_vup_geo_stitch, 'r+') as h5:
            h5.attrs['ref_stas'] = ref_stas
        logger.info ('Wrote: %s', Expi.path_vup_geo_stitch)
    return


## this is wrong; need to get coordinate of new ref sta
def proj_NYC2(ExpA, ExpB, npix=3):
    """ Project to East and West and then stitch together

    East:
        NJYM as base (-1.47 +/- 0.92; 2012 to 2018)
        NJBK  (-1.21  +/- 0.75; 2012 to 2021)
        NYBR  (-0.84 +/- 0.54; 2011 to present); jump
    West:
        NJHT as base (-1.69 +/- 0.73; 2015 to present)
        ROG1 (ok) (-1.57 +/- 1.25; 2011 to 2019)
        NYPR (ok) (-1.39 +/- 0.86; 2014 to 2021)

    """
    from BZ import bbGIS
    import geopandas as gpd
    # get the crops
    epsg     = 'epsg:4326'

    df_gps   = prep_gps(ExpA.path_gps, ExpA.reg, horiz=True)

    lst_rates, lst_uncs = [], []
    for i, Expi in enumerate([ExpA, ExpB]):
        if i == 0:
            assert Expi.ref_sta == 'NJHT'
            sta2 = 'ROG1'
        elif i == 1:
            assert Expi.ref_sta == 'NYJM'
            sta2 = 'NYBK'

        ## load the West/East data
        ## eventually write these to the individual MintPy directory
        path_vel = op.join(Expi.path_mp_exp_geo, f'geo_vlos_recon{Expi.neofs}.h5') \
            if Expi.neofs else Expi.path_vlos_geo

        rate, meta = readfile.read(path_vel, datasetName='velocity')
        unc  = readfile.read(path_vel, datasetName='velocityStd')[0]
        inc  = readfile.read(Expi.path_geom_mp_geo, datasetName='incidenceAngle')[0]
        y, x     = int(meta['REF_Y']), int(meta['REF_X'])

        ## convert it to LOS
        rate /= np.cos(np.deg2rad(inc))
        unc  /= np.cos(np.deg2rad(inc)) # bumps it up slightly

        # now iterate over the ref and extra station to get an offset
        offsets, weights = [], []
        for sta in [Expi.ref_sta, sta2]:
            rate_ref = df_gps[df_gps.sta == sta].u_vel.item()
            unc_ref  = df_gps[df_gps.sta == sta].u_sig.item()

            ## use the difference between the GPS rate and the pixels as offset
            offset = rate_ref - np.nanmean(rate[y-npix:y+npix, x-npix:x+npix])
            ## use the abs difference between GPS and pixels as weight
            weight = unc_ref**2 - np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])**2
            offsets.append(offset)
            weights.append(np.abs(weight))

        ## ~~~~~ PERFORM WEIGHTED MEAN OF OFFSETS USING GNSS ONLY ~~~~~ ##
        wgts       = 1/np.array(weights) # variance
        offset     = np.average(offsets, weights=wgts) # equivalent to polyfit
        offset_unc = np.average(np.sqrt(weights), weights=wgts)

        # _,     cov0  = np.polyfit(np.ones_like(offsets), offsets, 0, w=np.sqrt(wgts), cov=True)
        ## much higher cov2
        # offset1, cov1 = np.polyfit(np.ones_like(offsets), offsets, 0, w=np.sqrt(wgts), cov='unscaled')
        # offset_unc   = np.max([cov0, cov1])
        # offset_unc   = np.sqrt(offset_unc)
        ## equivalent to polyfit one with cov=True
        # cov2   = np.average((offsets-offset1)**2, weights=wgts)


        logger.info (f'Fit Offset: {offset*1000:.2f} +/- {offset_unc*1000:.2f} mm/yr')

        rate += offset
        unc   = np.sqrt(unc**2 + offset_unc**2)


        shutil.copy(path_vel, Expi.path_vup_geo)
        writefile.write_hdf5_block(Expi.path_vup_geo, rate, 'velocity')
        writefile.write_hdf5_block(Expi.path_vup_geo, unc, 'velocityStd')
        logger.info ('Wrote: %s', Expi.path_vup_geo)

        ## now do the netcdf stuff for stitching
        da_rate = h5_to_xr(Expi.path_vup_geo, 'velocity')
        da_rate.rio.write_crs(epsg, inplace=True)

        da_unc  = h5_to_xr(Expi.path_vup_geo, 'velocityStd')
        da_unc.rio.write_crs(epsg, inplace=True)

        ## clip the rates to polygon
        path_crop    = op.join(Expi.path_crops, f'{CROP_MAP[Expi.ref_sta]}.GeoJSON')
        gdf_crop     = gpd.read_file(path_crop)
        da_rate_crop = da_rate.rio.clip(gdf_crop.geometry, epsg, drop=False)
        da_unc_crop  = da_unc.rio.clip(gdf_crop.geometry, epsg, drop=False)

        lst_rates.append(da_rate_crop)
        lst_uncs.append(da_unc_crop)

    overlap = np.sum(lst_rates, 0)
    # nan out the overlap values in the western portion
    if overlap.any():
        log.warning('Using eastern values in overlap region (bad crops)')
        overlap_nan = np.where(np.isnan(overlap), 0, np.nan)
        lst_rates[0] += overlap_nan
        lst_uncs[0]  += overlap_nan

    rate_stitch = np.flipud(np.where(np.isnan(lst_rates[0].data),
                                     lst_rates[1].data, lst_rates[0].data))
    unc_stitch  = np.flipud(np.where(np.isnan(lst_uncs[0].data),
                                     lst_uncs[1].data, lst_uncs[0].data))

    ## just make both; they'll be the same
    for Expi in [ExpA, ExpB]:
        shutil.copy(Expi.path_vup_geo, Expi.path_vup_geo_stitch)
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, rate_stitch, 'velocity')
        writefile.write_hdf5_block(Expi.path_vup_geo_stitch, unc_stitch, 'velocityStd')
        # writefile.write_hdf5_block(Expi.path_vup_geo_stitch, res_stitch, 'residue')
        logger.info ('Wrote: %s', Expi.path_vup_geo_stitch)
    return


def proj_NYC3(ExpA, sta2='NJHT', npix=0):
    """ Project to one stations with slight offset using other GPS

    ## these to are very close
    NYBR  (-0.84 +/- 0.54; 2011 to present; there's a jump in 2021)
    NYBK  (-1.21  +/- 0.75; 2012 to 2021)
    """
    from BZ import bbGIS
    import geopandas as gpd
    # get the crops

    df_gps   = prep_gps(ExpA.path_gps, ExpA.reg, horiz=True)

    lst_rates, lst_uncs = [], []
    # assert ExpA.ref_sta == 'NYBK', 'Use reference NYBK'
    # sta2 = 'NYBR'

    ## eventually write these to the individual MintPy directory
    path_vel = op.join(ExpA.path_mp_exp_geo, f'geo_vlos_recon{ExpA.neofs}.h5') \
        if ExpA.neofs else ExpA.path_vlos_geo

    rate, meta = readfile.read(path_vel, datasetName='velocity')
    unc  = readfile.read(path_vel, datasetName='velocityStd')[0]
    inc  = readfile.read(ExpA.path_geom_mp_geo, datasetName='incidenceAngle')[0]
    y, x  = int(meta['REF_Y']), int(meta['REF_X'])
    coord = ut.coordinate(meta, lookup_file=ExpA.path_geom_mp_geo)

    ## convert it to LOS
    rate /= np.cos(np.deg2rad(inc))
    unc  /= np.cos(np.deg2rad(inc)) # bumps it up slightly

    # now iterate over the ref and extra station to get an offset
    offsets, weights = [], []
    for i, sta in enumerate([ExpA.ref_sta, sta2]):
        ser_sta  = df_gps[df_gps.sta == sta]
        rate_ref = ser_sta.u_vel.item()
        unc_ref  = ser_sta.u_sig.item()
        # unc_ref will be 0 at the real reference station, so then it's propagated later

        ## get the other GPS coordinate
        if i == 1:
            y, x = coord.geo2radar(ser_sta.lat, ser_sta.lon)[0:2]

        ## use the difference between the GPS rate and the pixels as offset
        if npix:
            offset = rate_ref - np.nanmean(rate[y-npix:y+npix, x-npix:x+npix])
            weight = unc_ref**2 - np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])**2
        else:
            offset = rate_ref - rate[y, x]
            weight = unc_ref**2 + unc[y,x]**2

        ## use the abs difference between GPS and pixels as weight
        offsets.append(offset)
        weights.append(np.abs(weight))

    ## ~~~~~ PERFORM WEIGHTED MEAN OF OFFSETS USING GNSS ONLY ~~~~~ ##
    wgts       = 1/np.array(weights) # variance
    offset     = np.average(offsets, weights=wgts) # equivalent to polyfit
    offset_unc = np.average(np.sqrt(weights), weights=wgts)

    logger.info (f'Fit Offset: {offset*1000:.2f} +/- {offset_unc*1000:.2f} mm/yr')

    rate += offset
    unc   = np.sqrt(unc**2 + offset_unc**2)

    shutil.copy(path_vel, ExpA.path_vup_geo)
    writefile.write_hdf5_block(ExpA.path_vup_geo, rate, 'velocity')
    writefile.write_hdf5_block(ExpA.path_vup_geo, unc, 'velocityStd')
    logger.info ('Wrote: %s', ExpA.path_vup_geo)

    epsg    = 'epsg:4326'
    da_rate = h5_to_xr(ExpA.path_vup_geo, 'velocity')
    da_rate.rio.write_crs(epsg, inplace=True)

    da_unc  = h5_to_xr(ExpA.path_vup_geo, 'velocityStd')
    da_unc.rio.write_crs(epsg, inplace=True)
    logger.info ('Wrote: %s', ExpA.path_vup_geo)

    return


def proj_one(ExpA, path_vel=None, npix=0, na12=False, force_ref=()):
    """ Tie to a single GPS station (optionally npix surrounding pixels)

    use force_ref (tuple:rate, unc; m) to manually put in a reference (Kiribati, TG-ALT)
    """
    ## identical to the other class version
    if path_vel is None:
        path_vel = op.join(ExpA.path_mp_exp_geo,
                           f'geo_velocity_recon{ExpA.neofs}.h5') if \
                            ExpA.neofs else ExpA.path_vlos_geo

    logger.info ('Projection velocity at: %s', path_vel)

    ## remove the station from the exclude when using e.g., stitching
    df_gps   = prep_gps(ExpA.path_gps, ExpA.reg, horiz=True)
    ser_ref  = df_gps[df_gps.sta == ExpA.ref_sta]

    assert not ser_ref.empty, f'No station matching {ExpA.ref_sta} in {ExpA.path_gps}'

    logger.info ('Tieing to: %s', ser_ref.sta.item())

    rate0, meta = readfile.read(path_vel, datasetName='velocity')
    unc  = readfile.read(path_vel, datasetName='velocityStd')[0]
    inc  = readfile.read(ExpA.path_geom_mp_geo, datasetName='incidenceAngle')[0]
    az   = readfile.read(ExpA.path_geom_mp_geo, datasetName='azimuthAngle')[0]

    # 73, 470
    y, x     = int(meta['REF_Y']), int(meta['REF_X'])
    if na12:
        rate_ref = utils0.enu2los(ser_ref.e_vel, ser_ref.n_vel, ser_ref.u_vel, inc[y, x], az_angle=az[y,x]).item()
        unc_ref  = utils0.enu2los(ser_ref.e_sig, ser_ref.n_sig, ser_ref.u_sig, inc[y, x], az_angle=az[y,x]).item()
        path_vup_geo = ExpA.path_vup_geo_na12

    ## residue will be copied over if it exists
    # res  = readfile.read(path_vel, datasetName='residue')[0]
    else:
        rate = rate0 / np.cos(np.deg2rad(inc))
        unc  /= np.cos(np.deg2rad(inc)) # bumps it up slightly
        rate_ref = ser_ref.u_vel.item()
        unc_ref  = ser_ref.u_sig.item()
        path_vup_geo = ExpA.path_vup_geo

    ### average a few pixels around the GPS station (at GPS, insar=0)
    if npix:
        rate_ref = np.nanmean(rate[y-npix:y+npix, x-npix:x+npix]) + rate_ref
        unc_ref  = np.nanmean(unc[y-npix:y+npix, x-npix:x+npix])  + unc_ref

    if force_ref:
        rate_ref, unc_ref = force_ref
        log.warning(f'Forcing reference rate/unc to {rate_ref:.2f}, {unc_ref:.2f}')

    logger.info (f'GPS Offset: {rate_ref*1000:.2f} +/- {unc_ref*1000:.2f} mm/yr')
    rate += rate_ref
    unc   = np.sqrt(unc**2 + unc_ref**2)

    # use na12 gps
    if na12:
        rate /= np.cos(np.deg2rad(inc))
        unc  /= np.cos(np.deg2rad(inc)) # bumps it up slightly

    # copy and update layers
    shutil.copy(path_vel, path_vup_geo)
    writefile.write_hdf5_block(path_vup_geo, rate, 'velocity')
    writefile.write_hdf5_block(path_vup_geo, unc, 'velocityStd')
    # writefile.write_hdf5_block(ExpA.path_vup_geo, res, 'residue')

    # add reference station to attributes
    with h5py.File(path_vup_geo, 'r+') as h5:
        h5.attrs['ref_stas'] = ser_ref.sta.item()
        h5.attrs['n_pix'] = npix
        h5.attrs['REF_RATE'] = rate_ref
        h5.attrs['REF_UNC'] = unc_ref

    logger.info ('Wrote: %s', path_vup_geo)
    return path_vup_geo


## not implemented
def proj_plane(ExpA, path_vel_los):
    from BZ import bbGIS
    from BZ.Tools.planeFit import PlaneFit
    import geopandas as gpd
    from mintpy.utils import readfile, writefile

    # get the GPS
    df_gps  = prep_gps(ExpA.path_gps, ExpA.reg)
    # print (GNSS[['sta', 'u_vel']])
    plane = FitPlane(df_gps.lon,df_gps.lat,df_gps.u_vel, df_gps.u_sig)(plot=False)[1]

    # get the insar rates
    epsg    = 'epsg:4326'
    da_rate = h5_to_xr(path_vel_los, 'velocity')
    da_rate.rio.write_crs(epsg, inplace=True)

    da_unc  = h5_to_xr(path_vel_los, 'velocityStd')
    da_unc.rio.write_crs(epsg, inplace=True)

    inc = readfile.read(ExpA.path_geom_mp_geo, 'incidenceAngle')[0]

    ## convert to vertical
    da_rate = (da_rate * (1/np.cos(inc))) - plane
    da_unc  = da_unc  * (1/np.cos(inc))

    rate = np.flipud(da_rate.data)
    unc  = np.flipud(da_unc.data)

    shutil.copy(path_vel_los, ExpA.path_vup_geo)
    writefile.write_hdf5_block(ExpA.path_vup_geo, rate, 'velocity')
    writefile.write_hdf5_block(ExpA.path_vup_geo, unc, 'velocityStd')
    return


if __name__ == '__main__':
    # Exp       = ExpBase(NYC_SR, 'Base', 'NJHT', neofs=15)
    # Exp       = ExpBase(NYC_SR, 'PM_Fast', 'NYBK', neofs=20)
    # proj_one(Exp, npix=0)

    # Exp       = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NJHT', neofs=20)
    # Exp1      = ExpBase(NYC_SRc, 'ERA5_SET_PM_Stitch_ex_Fast', 'NYBK', neofs=20)
    # Exp0      = ExpBase(NYC_SR, 'PM_Fast', 'NYBK', neofs=20)
    # proj_NYC(Exp, Exp1)
    # proj_NYC2(Exp, Exp1)

    Exp1 = ExpBase(DC_SR, 'ERA5_SET_PM_ex_Fast_2017', 'USN8', neofs=0)
    proj_one(Exp1, npix=0)


    # proj_all(Exp, path_vlos)

    ## mostly deprecated below here
    # path_vlos = op.join(Exp.path_mp_exp_geo, f'geo_vlos_recon{Exp.neofs}.h5')
    # Obj = RefProj(Exp.path_mp_exp_geo, Exp.reg, Exp.path_mask_mp_nc, corr=Exp.lbl)


    # inps = cmdLineParse()
    # ProjHR(inps.path_mp, inps.mask, inps.corr, show_cmds=True)(inps.method)

    # Exp = ExpBase(Savannah_SR)
    # Obj = RefProj(Exp.path_mp_exp_geo, Exp.reg, Exp.path_mask_mp_nc, corr=Exp.lbl)
    # Obj.proj_one()
    # Obj.proj_HR(350)
