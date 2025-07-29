""" Combine the unwrapped residual and unwrapped multilooked """
from __init__ import *
import re

def main(path_ml, path_resid, plot=False):
    ## upsample the multilooked unwrapped ifg to the full res wrapped ifg
    ds_ml    = gdal.Open(path_ml)
    ds_resid = gdal.Open(path_resid)

    ## make readable if it's not already; not robust
    ## this doesn't work; actually it might
    if ds_resid is None:
        from makeUnwVRT import write_xml
        ifg       = re.search('20[0-9][0-9]*_20[0-9][0-9]*', path_resid).group()
        path_wrap = op.join(op.dirname(path_ml), f'{ifg}.int')
        write_xml(path_resid, bands=2, path_wrap=path_wrap)
        ds_resid = gdal.Open(path_resid)

    ## upsample multilooked
    ds_unw_0 = gdal.Translate('', ds_ml, format='VRT', resampleAlg='lanczos',
                    width=ds_resid.RasterXSize, height=ds_resid.RasterYSize)

    ## sum results
    arr_unw0, arr_unw1  = ds_unw_0.ReadAsArray(), ds_resid.ReadAsArray()
    arr_new             = arr_unw0 + arr_unw1

    ## save as gdal readable
    dst    = path_resid.replace('resid', 'sum')
    ds_new = gdal.GetDriverByName('GTiff').CreateCopy(dst, ds_resid)
    ds_new.GetRasterBand(1).WriteArray(arr_new[0])
    ds_new.GetRasterBand(2).WriteArray(arr_new[1])
    # ds_out.GetRasterBand(1).SetNoDataValue(inps.nodata)

    print ('Combined multilooked and residual unwrapped ifgs:', dst)
    if plot:
        arrs      = [arr_unw0, arr_unw1, arr_new]
        tis       = ['Upsampled Unwrapped', 'Residual Unwrapped', 'Sum Unwrapped']
        for i, kind in enumerate(['Amplitude', 'Phase']):
            cm  = 'binary' if i == 0 else 'coolwarm'
            fig, axes = plt.subplots(figsize=(10,10), nrows=3, sharex=True)
            for j, ax in enumerate(axes.ravel()):
                arr = arrs[j][i]
                im  = ax.imshow(arr, cmap=cm, origin='lower', interpolation='nearest')
                imshow_cbar(im, pad=0.2)
                ax.set_title(f'{tis[j]} {kind}')
            fig.set_label(f'Combined_{kind}')
        savefigs(op.dirname(dst), True, True)

    return ds_new


if __name__ == '__main__':
    root      = os.getenv('dataroot')
    path_root = op.join(root, *'VLM Sentinel1 HR_Test_steps_local'.split())
    path_wd   = op.join(path_root, 'PS_DS2', 'unwrap')

    # unwrapped, multilooked
    path_unw_ml = op.join(path_wd, '20190124_20190217_masked.7alks_19rlks.unw')

    # unwrapped, full resolution, resid
    path_unw_fr = op.join(path_wd, '20190124_20190217_masked_resid.unw')

    main(path_unw_ml, path_unw_fr, plot=True)
    plt.show()
