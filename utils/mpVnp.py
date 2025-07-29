#!/usr/bin/env python

"""
Compare the trend from a timeseries computed by mintpy vs numpy polyfit
"""
import argparse
from VLM.bzFRInGE import *
from mintpy.objects import timeseries
from mintpy.utils import readfile
from BZ.bbTS import date2dec, bzSeason

log = logging.getLogger('BZ')

def createParser():
    parser = argparse.ArgumentParser(description='Compare the trend from a timeseries computed by mintpy vs numpy polyfit')
    parser.add_argument('path_ts', type=str, help='path to timeseries h5')
    parser.add_argument('path_vel', type=str, help='path to velocity h5')
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter the debugger on completion')
    parser.add_argument('--yx', default='449,1089', help='pixel coordinate as y,x')
    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)


def main(inps):
    tsObj = timeseries(inps.path_ts)
    arr_ts = tsObj.read() * 1000

    dates  = pd.to_datetime(tsObj.get_date_list())
    decyr  = date2dec(dates)
    y, x   = [int(yx) for yx in inps.yx.split(',')]

    arr_vel, meta = readfile.read(inps.path_vel)
    arr_vel *= 1000
    if meta['mintpy.timeFunc.polynomial'] == '1' and \
        meta['mintpy.timeFunc.periodic'] == '[]':
        coeffs = np.polyfit(decyr, arr_ts[:, y,x], 1)
        rate   = coeffs[0]
    elif meta['mintpy.timeFunc.polynomial'] == '1' and \
        meta['mintpy.timeFunc.periodic'] == '[1.0, 0.5]':
        log.info('Fitting periodic...')
        rate = bzSeason.fit_lsq(decyr, arr_ts[:, y, x])[1][-2]

    else:
        raise Exception('Other models not supported yet.')

    los = 'Vup' if 'Vup' in inps.path_vel else 'LoS'
    log.info(f'Rate from Velocity {los} File": {arr_vel[y, x]:.3f}')
    log.info(f'Rate from Timeseries LoS File": {rate:.3f}')
    log.info(f'Residual: {arr_vel[y,x]-rate:.3f}')

    inc = readfile.read('./geo_geometryRadar.h5', datasetName='incidenceAngle')[0]
    inc1 = inc[y, x]
    gps_vel = -1.21
    cos_inc = np.cos(np.deg2rad(inc1))

    if 'Vup' in inps.path_vel:
        log.warning(f'Undo Vup: {(arr_vel[y,x]-gps_vel)*cos_inc:.3f}')
        # convert the rate from the timeseries (LOS) to Vup
        log.warning(f'New Residual Vup: {((arr_vel[y,x]-gps_vel)*cos_inc)-rate:.3f}')

        ## now try and use the GPS trend

        gps_trend  = np.polyval([gps_vel, 0], decyr)
        gps_trend -= gps_trend[0]
        ts_vup = arr_ts[:, y, x]/cos_inc + gps_trend
        coeffs = np.polyfit(decyr, ts_vup, 1)
        rate_vup_ts = coeffs[0]

        print (f'MP Vup: {arr_vel[y,x]:.3f}')
        print (f'TS Vup: {rate_vup_ts:.3f}')
        print (f'Residual by hand Vup: {(arr_vel[y,x]-rate_vup_ts):.3f}')

    ## calculate plate motion differnece
    arr_itrf_los = readfile.read('geo_velocity_recon20_ITRF14.h5')[0]
    arr_noitrf_los = readfile.read('geo_velocity_recon20_noitrf.h5')[0]
    arr_pm_vup = 1000*(arr_itrf_los-arr_noitrf_los)[y, x] / cos_inc
    print (arr_pm_vup)

    if 'PM' in inps.path_vel.upper():
        arr_pm0, meta1 = readfile.read('ITRF14.h5')
        assert 'REF_X' in meta1.keys(), 'Need to reference the plate motion'
        arr_pm1 = (1000 * arr_pm0 / np.cos(np.deg2rad(inc))).reshape(-1)
        inters = np.zeros_like(arr_pm1)
        # arr_pm = np.polyval([arr_pm.reshape(-1), np.zeros_like(arr_pm.reshape(-1))], decyr)
        arr_pm1 = (arr_pm1[:, np.newaxis] * decyr) + inters[:, np.newaxis]
        arr_pm  = arr_pm1.transpose(1, 0).reshape(arr_ts.shape)
        print (np.polyfit(decyr, arr_pm[:, y, x], 1)[0])
    breakpoint()



    # print('Vup:', (arr_vel[y, x]/cos_inc)+gps_vel)

if __name__ == '__main__':
    inps = cmdLineParse()
    assert os.path.exists(inps.path_ts), 'Incorrect path to timeseries file'
    assert os.path.exists(inps.path_vel), 'Incorrect path to velocity file'
    main(inps)

    if inps.interact:
        try:
            embed(colors='neutral')
        except:
            breakpoint()