#!/usr/bin/env python
from __init__ import *

def createParser():
    parser = argparse.ArgumentParser(description='Apply a mask to another file (need to be same shape).',
                                     epilog='Examples of use:\n\t applyAmpMask.py ./PS_DS/20190124_20190217.int ./masks/ampMask.msk'\
                                             '\n\t applyAmpMask.py ./PS_DS/20150310_20211028_ISCE.2alks_5rlks.int ./masks/waterMask_2alks_5rlks.envi',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path_input',  type=str, help='Path of file to mask')
    parser.add_argument('path_mask',  type=str, help='Path of maskfile')
    parser.add_argument('--nodata', type=float, default=0, help='nodata value for water. Default is 0.')
    parser.add_argument('-v', '--view', action='store_true', help='View the mask')
    parser.add_argument('-s', '--save', action='store_true', help='Save to maskAmp.png')
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter an interactive session on completion')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 3:
        parser.print_help()
        os.sys.exit(1)
    inps   = parser.parse_args(args=iargs)
    return inps


def main(inps):
    """ Mask a wrapped interferogram or temporal correlation """
    assert op.exists(inps.path_input), f'{inps.path_input} does not exist'
    assert op.exists(inps.path_mask), f'{inps.path_mask} does not exist'

    ## to mask
    ds2m    = gdal.Open(inps.path_input)
    arr2m   = ds2m.ReadAsArray()
    # gdal_dt = gdal.GetDataTypeName(ds2m.GetRasterBand(1).DataType)

    ## mask
    dsm  = gdal.Open(inps.path_mask)
    arrm = dsm.ReadAsArray()

    ## nan out, then replace with nodata
    # nodata     = inps.nodata+inps.nodata*1j if 'CFloat' in gdal_dt else inps.nodata
    arr_masked = np.where(np.isclose(arrm, 0), inps.nodata, arr2m)

    dst    = '{}_masked{}'.format(*op.splitext(inps.path_input))

    ds_out = gdal.GetDriverByName('ENVI').CreateCopy(dst, ds2m)
    ds_out.GetRasterBand(1).WriteArray(arr_masked)
    # ds_out.GetRasterBand(1).SetNoDataValue(inps.nodata)
    print ('Wrote masked file:', dst)

    ## quick plot
    fig, axes = plt.subplots(figsize=(9,9), nrows=3, sharex=True)
    arrs = [arrm, arr2m, arr_masked]
    tis  = 'Mask Original AppliedMask'.split()
    for ax, ti, arr in zip(axes.ravel(), tis, arrs):
        arr  = arr.real
        im   = ax.imshow(arr, interpolation='nearest', origin='lower',
                            cmap='binary')
        if ti == 'Mask':
            ticks = np.array([0, 1])
            imshow_cbar(im, **dict(ticks=ticks, pad=0.25)).set_ticks(ticks)
        else:
            imshow_cbar(im, pad=0.25)

    fig.set_label(op.basename(dst))
    if inps.save or (inps.view and BACKEND == 'agg'):
        savefigs(inps.path_input, figures=True, overwrite=True)

    embed(colors='neutral') if inps.interact else ''
    return dst


if __name__ == '__main__':
    inps = cmdLineParse()
    main(inps)


    plt.show() if inps.view and not BACKEND == 'agg' else ''
