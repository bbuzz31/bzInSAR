import re
from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *


def split_masks_HR(path_aria_wmask, path_cmask, dst_dir=None, show=True):
    """ Combines the water mask and crop mask for prep_aria

    Must use ARIA water mask that's correctly cropped and shaped
    just cropping with raw mask with gdal may not give right shape,
    as bounds adjusted by ARIA
    """
    path_crops = op.dirname(path_cmask)
    fname = f'{op.splitext(op.basename(path_aria_wmask))[0]}_'\
            f'{op.splitext(op.basename(path_cmask))[0]}.envi'

    dst_dir = op.dirname(path_aria_wmask) if dst_dir is None else dst_dir
    dst     = op.join(dst_dir, fname)

    wmask = gdal.Open(path_aria_wmask, gdal.GA_ReadOnly)
    ds    = gdal.Warp('', wmask, format='VRT', cutlineDSName=path_cmask, dstNodata=0)
    ds.FlushCache()

    # unset the nodata so that 0 is read as water and 1 as land by mintpy
    ds    = gdal.Translate(dst, ds, format='ENVI', noData=-9999)
    arr   = ds.ReadAsArray()

    if show:
        fig, axes = plt.subplots()
        axes.imshow(arr, interpolation='nearest', cmap='binary')

    print (f'Wrote: {dst} {arr.shape}')
    del ds, wmask

    return dst


def read_spa_coh(path_txt):
    pass


def get_weather_model(mp_exp):
    """ Get the raider weather model if it exists"""
    models = 'GMAO ERA5 ERA-5 HRES HRRR'.split()
    wm  = ''
    for model in models:
        if model.lower() in mp_exp.lower():
            wm = model
            break
    return wm


class ARIABase(object):
    """ base class for setting up experiments """
    def __init__(self, dct_exp, mp_exp='Base', ref_sta=None, v=True):
        self.dct_exp = dct_exp
        self.reg     = re.match('[A-Za-z]+', self.dct_exp['root']).group()
        self.mp_exp  = mp_exp # e.g., Base
        self.ref_sta = DCT_REG[self.reg][3][0] if ref_sta is None else ref_sta
        self.corr    = 'ARIA'
        self.wm      = get_weather_model(mp_exp)
        self.neofs   = '' # placeholder
        self._make_meta()
        self._make_struct()
        self._set_paths_base()
        self._set_paths_vup()


    def _make_meta(self, dct_exp=None):
        dct_exp    = self.dct_exp if dct_exp is None else dct_exp
        reg        = re.match('[A-Za-z]+', dct_exp['root']).group()

        self.SNWE  = DCT_REG[reg][0]
        self.SNWEs = ' '.join([str(i) for i in self.SNWE])
        self.track = str(DCT_REG[reg][1])
        self.frame = '_'.join([str(i) for i in DCT_REG[reg][2]])
        self.post  = '-0.00083333 0.00083333'
        return


    def _make_struct(self, dct_exp=None):
        dct_exp         = self.dct_exp if dct_exp is None else dct_exp
        self.path_vlm   = op.join(os.getenv('dataroot'), 'VLM')
        self.path_wd    = op.join(self.path_vlm, 'Sentinel1', dct_exp['root'])
        self.path_crop  = op.join(op.dirname(self.path_vlm), 'Crops', self.reg)
        self.path_masks = op.join(self.path_wd, 'mask')
        self.path_mask  = op.join(self.path_masks, 'OSM_wmask.msk')
        self.path_prods = op.join(self.path_wd, 'products')

        for path in [self.path_wd, self.path_masks, self.path_prods]:
            os.makedirs(path, exist_ok=True)
        return


    def _set_paths_base(self, dct_exp=None):
        dct_exp       = self.dct_exp if dct_exp is None else dct_exp
        self.path_vlm = op.join(os.getenv('dataroot'), 'VLM')
        self.path_gps = op.join(self.path_vlm, 'GNSS', 'UNR', f'Midas_GPS-{self.reg}_cln.csv')

        return


    def _set_paths_vup(self, dct_exp=None):
        dct_exp           = self.dct_exp if dct_exp is None else dct_exp
        self.path_mp_exp      = op.join(self.path_wd, f'{self.ref_sta}_{self.mp_exp}')
        self.path_mask_mp_geo = op.join(self.path_mp_exp, 'waterMask.h5')

        self.path_ts          = op.join(self.path_mp_exp, 'timeseries.h5')
        self.path_ts_dem      = op.join(self.path_mp_exp, 'timeseries_demErr.h5')
        self.path_ts_geo      = self.path_ts_dem # for consistency with FRInGE
        self.path_ts_tropo    = op.join(self.path_mp_exp, 'timeseries_demErr_HRRR.h5')
        self.path_ifgramStack = op.join(self.path_mp_exp, 'inputs', 'ifgramStack.h5')
        self.path_geometryGeo = op.join(self.path_mp_exp, 'inputs', 'geometryGeo.h5')
        self.path_geom_mp_geo = self.path_geometryGeo
        self.path_tropoGeo    = op.join(self.path_mp_exp, 'inputs', f'{self.wm}.h5')
        self.path_geom_mp_geo = self.path_geometryGeo

        self.path_tc          = op.join(self.path_mp_exp, 'temporalCoherence.h5')
        self.path_mask_tc     = op.join(self.path_mp_exp, 'maskTempCoh.h5')
        self.path_mask_cc     = op.join(self.path_mp_exp, 'maskConnComp.h5')
        self.path_vel         = op.join(self.path_mp_exp, 'velocity.h5')
        self.path_vlos        = self.path_vel
        self.path_vlos_geo    = self.path_vlos # consistency w/ FRInGE

        ## paths for referencing, projecting, plotting
        # self.path_mp_vup      = op.join(op.dirname(self.path_mp_exp), f'Vup_{op.basename(self.path_mp_exp)}')
        self.path_mp_vup      = op.join(self.path_mp_exp, 'Vup')

        self.path_mask_vup    = op.join(self.path_mp_vup, 'waterMask.h5')
        self.path_mask_vup_nc = op.join(self.path_mp_vup, 'geo_waterMask.nc')

        self.path_vup_geo     = op.join(self.path_mp_vup, f'Vup_{self.mp_exp}_ARIA.h5')
        self.path_rate_nc     = op.join(self.path_mp_vup, f'geo_rate_ARIA.nc')
        self.path_std_nc      = op.join(self.path_mp_vup, f'geo_std_ARIA.nc')

        # should be automatically masked
        self.path_rate_msk_nc = op.join(self.path_mp_vup, f'geo_rate_ARIA.nc')
        self.path_std_msk_nc  = op.join(self.path_mp_vup, f'geo_std_ARIA.nc')

        self.path_vel_kmz     = op.splitext(self.path_rate_nc)[0]
        self.path_std_kmz     = op.splitext(self.path_std_nc)[0]
        # self.path_resid_kmz   = op.splitext(self.path_resid_nc)[0]

        # for differences made for comparing experiments
        self.path_cmp = op.join(self.path_wd, 'comparisons')

        os.makedirs(self.path_mp_vup, exist_ok=True)
        os.makedirs(self.path_cmp, exist_ok=True)


if __name__ == '__main__':
    ARIABase(HR_ARIA)
