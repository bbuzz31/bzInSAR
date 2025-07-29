#!/usr/bin/env python
from __init__ import *

def createParser():
    parser = argparse.ArgumentParser(description='Use gdal_fillnodata to fill nans in wrapped IFGs',
                                     epilog='Examples of use:\n\t gapFillWrapped.py ./PS_DS',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path_ifgs', type=str, help='Path to PS_DS')
    parser.add_argument('path_mask', type=str, help='Path to mask')
    parser.add_argument('--naz', default=7, type=int, help='Number of azimuth looks (default=7)')
    parser.add_argument('--nrg', default=19, type=int, help='Number of range looks (default=19)')
    # parser.add_argument('--annual', action='store_true', help='Only do annual pairs')
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter an interactive session on completion')
    parser.add_argument('--si', default=0, type=int, help='Optionally perform smoothing')
    parser.add_argument('--overwrite', default=1, type=int, help='Overwrite existing gap filled files. Default is True.')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(1)

    inps = parser.parse_args(args=iargs)
    return inps


def progress_cb(complete, message, cb_data):
    """ Emit progress report in numbers for 10% intervals and dots for 3%

    https://stackoverflow.com/questions/68025043/adding-a-progress-bar-to-gdal-translate
    """
    if int(complete*100) % 10 == 0:
        print(f'{complete*100:.0f}', end='', flush=True)
    elif int(complete*100) % 3 == 0:
        print(f'{cb_data}', end='', flush=True)
    return


def main(inps):
    if op.isdir(inps.path_ifgs):
        print (f'Searching {inps.path_ifgs} for wrapped ifgs')
        if inps.naz > 1 or inps.nrg > 1:
            path_ifgs = glob.glob(op.join(inps.path_ifgs,
                                    f'20*{inps.naz}alks_{inps.nrg}rlks.int'))
        else:
            path_ifgs = glob.glob(op.join(inps.path_ifgs, '20*[0-9].int'))

        print (f'gapFilling {len(path_ifgs)} wrapped ifgs.')

    else:
        path_ifgs = [inps.path_ifgs]


    ## use the mask file for matching the sizes
    mask     = gdal.Open(inps.path_mask)
    maskband = mask.GetRasterBand(1) # for gdal
    maskshp  = maskband.ReadAsArray().shape

    for i, path_ifg in enumerate(path_ifgs):
        dst  = '{}.gf{}'.format(*op.splitext(path_ifg))
        if op.exists(dst) and not inps.overwrite:
            continue

        if i % 20 == 0:
            print (f'gapFilling ifg {i} of {len(path_ifgs)} ifgs.')
        ds_ifg  = gdal.Open(path_ifg)
        arr_ifg = ds_ifg.ReadAsArray() # complex
        assert arr_ifg.shape == maskshp, 'ifg shape != mask shape'
        # phase   = np.angle(arr_ifg)  # not sure how to convert back to complex
        arrs    = []
        for kind in 'real imag'.split():
            ds   = gdal.GetDriverByName('MEM').CreateCopy('', mask)
            band = ds.GetRasterBand(1)
            if kind == 'real':
                band.WriteArray(arr_ifg.real)
            else:
                band.WriteArray(arr_ifg.imag)
            ds.FlushCache()

            gdal.FillNodata(targetBand=band, maskBand=maskband, maxSearchDist=1500,
                    smoothingIterations=inps.si,  callback=progress_cb, callback_data='.')


            # do this here so gdal messages print to new lines
            # print (f'\nFinished {kind}')
            arrs.append(ds.ReadAsArray())
            print ('')
            del ds

        arr_fill = arrs[0] + 1j*arrs[1]
        ## now save the ifg
        ds   = gdal.GetDriverByName('ENVI').CreateCopy(dst, ds_ifg)
        ds.GetRasterBand(1).WriteArray(arr_fill)
        del ds
        # print ('Wrote:', dst)

    embed(colors='neutral') if inps.interact else ''
    return


if __name__ == '__main__':
    inps = cmdLineParse()

    main(inps)
