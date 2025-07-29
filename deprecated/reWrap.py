import re
from __init__ import *


def main(path_unw_ml, path_wrap_fr, path_wrap_ml=None, plot=True):
    ## rewrap the multilooked, unwrapped ifg
    ds_unw_ml     = gdal.Open(path_unw_ml)
    if ds_unw_ml is None:
        raise Exception('File is not gdal readble; use makeUnwVrt.py')
    arr_unw_ml    = ds_unw_ml.GetRasterBand(2).ReadAsArray()
    arr_rewrap_ml = np.arctan2(np.sin(arr_unw_ml), np.cos(arr_unw_ml)) / np.pi

    ## create a gdal dataset out of the numpy rewrapped multilooked ifg
    ds_rewrap_ml  = gdal.GetDriverByName('MEM').Create('', arr_rewrap_ml.shape[1],
                                arr_rewrap_ml.shape[0], 1, gdalconst.GDT_Float32)
    ds_rewrap_ml.GetRasterBand(1).WriteArray(arr_rewrap_ml)

    ## upsample the multilooked rewrapped ifg to the full res wrapped ifg w/ gdal
    ds_wrap_fr   = gdal.Open(path_wrap_fr) # for full size
    ds_rewrap_fr = gdal.Translate('', ds_rewrap_ml, format='VRT',
                    width=ds_wrap_fr.RasterXSize, height=ds_wrap_fr.RasterYSize,
                                                        resampleAlg='lanczos')

    ## subtract the full resolution rewrapped from full resolution wrapped
    # maintain amplitude from original (lost in unwrapping?)
    arr_rewrap_fr = ds_rewrap_fr.ReadAsArray()
    arr_wrap_fr   = ds_wrap_fr.ReadAsArray()
    arr_phs_resid = arr_wrap_fr.imag - arr_rewrap_fr
    arr_resid     = arr_wrap_fr.real + arr_phs_resid*1j

    ## write the new wrapped to an ISCE object for snaphu
    dst    = '{}_resid{}'.format(*op.splitext(path_wrap_fr))
    ds_out = gdal.GetDriverByName('ISCE').CreateCopy(dst, ds_wrap_fr)
    ds_out.GetRasterBand(1).WriteArray(arr_resid)

    nodata = ds_wrap_fr.GetRasterBand(1).GetNoDataValue()
    if nodata is not None:
        ds_out.GetRasterBand(1).SetNoDataValue()

    print ('Wrote residual wrapped interferogram:', dst)
    if plot:
        if path_wrap_ml is not None:
            arr_wrap_ml   = gdal.Open(path_wrap_ml).ReadAsArray().imag
            arrs = [arr_unw_ml, arr_wrap_ml, arr_rewrap_ml,
                                        np.abs(arr_wrap_ml - arr_rewrap_ml)]
            tis  = ['Unwrapped', 'Wrapped', 'Rewrapped', 'Residual Wrapped']
            fig, axes = plt.subplots(figsize=(11,10), nrows=2, ncols=2,
                                                    sharex=True, sharey=True)
            for i, ax in enumerate(axes.ravel()):
                cm = 'coolwarm' if i+1 < len(axes.ravel()) else 'afmhot'
                ax.set_title(f'{tis[i]} Phase')
                im = ax.imshow(arrs[i], interpolation='nearest', cmap=cm,
                                                     origin='lower')
                imshow_cbar(im, pad=0.04)
            sti = re.search('\dalks.*rlks', path_wrap_ml).group()
            fig.suptitle(f'Multilooked: {sti}', y=0.925)
            fig.set_label('Multilooked_Rewrapping')

        arrs = [arr_wrap_fr.imag, arr_rewrap_fr, ds_out.ReadAsArray().imag]
        tis  = ['Wrapped', 'Rewrapped', 'Residual Wrapped']
        fig, axes = plt.subplots(figsize=(11,10), nrows=len(arrs), sharex=True)
        for i in range(len(axes)):
            cm = 'coolwarm'
            axes[i].set_ylabel(f'{tis[i]} Phase')
            im  = axes[i].imshow(arrs[i], interpolation='nearest', cmap=cm, origin='lower')
            imshow_cbar(im, pad=0.04)
        fig.suptitle('Full Resolution', y=0.925)
        fig.set_label('Full_Resolution_Rewrapping')
        savefigs(op.dirname(path_wrap_fr), True, True)

    del ds_rewrap_ml, ds_wrap_fr, ds_rewrap_fr, ds_out
    return


def run_all(path_wrap_dir, nrg=5, naz=2):
    files = os.listdir(path_wrap_dir)
    for f in files:
        ifg  = re.match(f'20[0-9]*_20[0-9]*_ISCE[.]{naz}alks_{nrg}rlks.int$', f)
        if ifg is not None:
            wrap_ml   = ifg.group()
            wrap_fr   = '{}_{}.int'.format(*wrap_ml.split('_')[:2])
            path_unw_ml = op.join(path_wrap_dir, 'unwrap', wrap_ml.replace('int', 'unw'))
            main(path_unw_ml, op.join(path_wrap_dir, wrap_fr), plot=False)
    return


def main_new(path_unw_ml, path_wrap_fr, path_wrap_ml=None, plot=True):
    ## rewrap the multilooked, unwrapped ifg
    ds_unw_ml     = gdal.Open(path_unw_ml)
    if ds_unw_ml is None:
        raise Exception('File is not gdal readble; use makeUnwVrt.py')
    # band 1 is coherence?
    arr_unw_ml    = ds_unw_ml.GetRasterBand(2).ReadAsArray()
    y = np.sin(arr_unw_ml)
    x = np.cos(arr_unw_ml)

    arr_rewrap_ml = np.arctan2(y, x)

    arr_unw_ml1 = np.arctan2(np.sin(-arr_rewrap_ml), -np.cos(-arr_rewrap_ml))
    print (np.allclose(arr_unw_ml, arr_unw_ml1))
    plt.imshow(np.abs(arr_unw_ml - arr_unw_ml1))
    plt.colorbar(); plt.show()
    X()

    # how to undo this... operation
    # magnitude     = np.sqrt(dx**2 + dy**2)

    ##
    arr_wrap_ml   = gdal.Open(path_wrap_ml).ReadAsArray()
    arr_wrap_ml   = np.angle(arr_wrap_ml)

    # identical
    # arr_wrap_ml  = np.arctan2(arr_wrap_ml.imag, arr_wrap_ml.real)



    ## create a gdal dataset out of the numpy rewrapped multilooked ifg
    ds_rewrap_ml  = gdal.GetDriverByName('MEM').Create('', arr_rewrap_ml.shape[1],
                                arr_rewrap_ml.shape[0], 1, gdalconst.GDT_Float32)
    ds_rewrap_ml.GetRasterBand(1).WriteArray(arr_rewrap_ml)

    ## upsample the multilooked rewrapped ifg to the full res wrapped ifg w/ gdal
    ds_wrap_fr   = gdal.Open(path_wrap_fr) # for full size
    ds_rewrap_fr = gdal.Translate('', ds_rewrap_ml, format='VRT',
                    width=ds_wrap_fr.RasterXSize, height=ds_wrap_fr.RasterYSize,
                                                        resampleAlg='lanczos')

    ## subtract the full resolution rewrapped from full resolution wrapped
    # maintain amplitude from original (lost in unwrapping?)
    arr_rewrap_fr = ds_rewrap_fr.ReadAsArray() # floats
    arr_wrap_fr   = ds_wrap_fr.ReadAsArray()   # complex
    embed(colors='neutral'); X()
    # arr_phs_resid = arr_wrap_fr.imag - arr_rewrap_fr
    # arr_resid     = arr_wrap_fr.real + arr_phs_resid*1j

    # ## write the new wrapped to an ISCE object for snaphu
    # dst    = '{}_resid{}'.format(*op.splitext(path_wrap_fr))
    # ds_out = gdal.GetDriverByName('ISCE').CreateCopy(dst, ds_wrap_fr)
    # ds_out.GetRasterBand(1).WriteArray(arr_resid)

    nodata = ds_wrap_fr.GetRasterBand(1).GetNoDataValue()
    if nodata is not None:
        ds_out.GetRasterBand(1).SetNoDataValue()

    print ('Wrote residual wrapped interferogram:', dst)
    if plot:
        if path_wrap_ml is not None:
            arr_wrap_ml   = gdal.Open(path_wrap_ml).ReadAsArray().imag
            arrs = [arr_unw_ml, arr_wrap_ml, arr_rewrap_ml,
                                        np.abs(arr_wrap_ml - arr_rewrap_ml)]
            tis  = ['Unwrapped', 'Wrapped', 'Rewrapped', 'Residual Wrapped']
            fig, axes = plt.subplots(figsize=(11,10), nrows=2, ncols=2,
                                                    sharex=True, sharey=True)
            for i, ax in enumerate(axes.ravel()):
                cm = 'coolwarm' if i+1 < len(axes.ravel()) else 'afmhot'
                ax.set_title(f'{tis[i]} Phase')
                im = ax.imshow(arrs[i], interpolation='nearest', cmap=cm,
                                                     origin='lower')
                imshow_cbar(im, pad=0.04)
            sti = re.search('\dalks.*rlks', path_wrap_ml).group()
            fig.suptitle(f'Multilooked: {sti}', y=0.925)
            fig.set_label('Multilooked_Rewrapping')

        arrs = [arr_wrap_fr.imag, arr_rewrap_fr, ds_out.ReadAsArray().imag]
        tis  = ['Wrapped', 'Rewrapped', 'Residual Wrapped']
        fig, axes = plt.subplots(figsize=(11,10), nrows=len(arrs), sharex=True)
        for i in range(len(axes)):
            cm = 'coolwarm'
            axes[i].set_ylabel(f'{tis[i]} Phase')
            im  = axes[i].imshow(arrs[i], interpolation='nearest', cmap=cm, origin='lower')
            imshow_cbar(im, pad=0.04)
        fig.suptitle('Full Resolution', y=0.925)
        fig.set_label('Full_Resolution_Rewrapping')
        savefigs(op.dirname(path_wrap_fr), True, True)

    del ds_rewrap_ml, ds_wrap_fr, ds_rewrap_fr, ds_out
    return

if __name__ == '__main__':
    root      = os.getenv('dataroot')
    path_root = op.join(root, *'VLM Sentinel1 HR'.split())
    path_wd   = op.join(path_root, 'Unwrap_Exps')
    # run_all(path_wd)


    # unwrapped, multilooked
    path_unw_ml  = op.join(path_wd, 'unwrap', '20150310_20211028_ISCE.2alks_5rlks_filled.unw')

    # original wrapped for comparison
    path_wrap_ml = op.join(path_wd, '20150310_20211028_ISCE.2alks_5rlks_filled.int')


    # full resolution wrapped
    path_wrap_fr = op.join(path_wd, '20150310_20211028.int')


    # main_new(path_unw_ml, path_wrap_fr, path_wrap_ml, plot=False)
    # plt.show()
    import numpy as np
    theta = np.random.uniform(low=-3*np.pi, high=3*np.pi, size=400)
    theta = theta.reshape(20,20)

    y     = np.sin(theta)
    x     = np.cos(theta)

    phi = -np.arctan2(y, x)


    x = np.cos(-theta)

    theta1 = np.arctan2(y, x)
    print (np.allclose(theta, theta1))
