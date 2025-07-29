""" Quick and dirty projection to Vup """
import shutil
import h5py
import pandas as pd
from __init__ import *

def main(wd, path_gnss, ref_sta='LOY2', use_mask=True):
    """ use_mask can be a path to a mask or a bool to use the waterMask.h5 """
    ## get GPS rate
    df_sta = pd.read_csv(path_gnss, index_col=None)
    df_sta = df_sta[df_sta.sta == ref_sta]
    gps_vup = df_sta.u_vel

    ## get velocity
    path_vel = op.join(wd, 'velocity.h5')
    with h5py.File(path_vel, 'r') as h5:
        ins_vlos = h5['velocity'][:]
        units    = h5.attrs['UNIT'][0]
        if units == 'm':
            ins_vlos *= 1000
        # get the reference point in radar coords; dont actually need
        ref_x = int(h5.attrs['REF_X'])
        ref_y = int(h5.attrs['REF_Y'])

    # vel_ref   = ins_vlos[ref_y, ref_x] # 0

    ## get inc angle
    path_geom = op.join(wd, 'inputs', 'geometryRadar.h5')
    with h5py.File(path_geom, 'r') as h5:
        inc   = h5['incidenceAngle'][:]

    ins_vup     = ins_vlos / np.cos(np.deg2rad(inc))
    ins_gps_vup = ins_vup + gps_vup.item()

    ## try and get the waterMask
    if use_mask:
        src = op.join(wd, 'waterMask.h5') if isinstance(use_mask, bool) else use_mask

        if op.exists(src):
            print ('Found mask:', src)
            with h5py.File(src, 'r') as h5:
                try:
                    mask = h5['mask'][:]
                    mask = np.where(np.isclose(mask, 0), np.nan, mask)
                except Exception as E:
                    print ('Could not load "mask" dataset in:', src)
                    print (E)
                    mask = np.ones_like(inc)
        else:
            mask = np.ones_like(inc)
    else:
        mask = np.ones_like(inc)

    ## write array as mintpy readable object
    dst         = op.join(wd, 'Vup.h5')
    shutil.copy(path_vel, dst)
    with h5py.File(dst, 'r+') as h5_new:
        data    = h5_new['velocity']
        data[:] = mask*ins_gps_vup/1000

        # just for masking
        data    = h5_new['velocityStd']
        data[:] = mask * data[:]
        h5_new.attrs['FILE_PATH'] = op.join(dst)

    print ('Wrote combined Vup to:', dst)
    return

if __name__ == '__main__':
    vlm_dir = op.join(os.getenv('dataroot'), 'VLM')
    wd      = op.join(vlm_dir, 'Sentinel1', 'HR', 'MintPy_11_7')#_2alks_5rlks')
    path_gnss = op.join(vlm_dir, *'GNSS UNR Midas_GPS-HR_cln.csv'.split())
    main(wd, path_gnss)
