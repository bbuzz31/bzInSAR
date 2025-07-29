#!/usr/bin/env python
from VLM.bzFRInGE import *
from isce.applications import looks

def createParser():
    parser = argparse.ArgumentParser(description='Multilook wrapped interferogram (or geometry files) using ISCE. If given a whole directory will do wrapped ifgs',
                                     epilog='Examples of use:\n\t multilookWrapped.py ./PS_DS/20190124_20190217.int' \
                                                             '\n\t multilookWrapped.py ./PS_DS/',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_file', type=str, help='Path to wrapped IFG')
    parser.add_argument('-m', '--mask', action='store_true', help='Force the use of mask file (will try to use extension though).')
    parser.add_argument('--naz', default=7, help='Number of azimuth looks (default=7)')
    parser.add_argument('--nrg', default=19, help='Number of range looks (default=19)')
    parser.add_argument('-v', '--view', action='store_true', help='View the plot')
    parser.add_argument('-s', '--save', action='store_true', help='Save to src.png')
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter an interactive session on completion')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        print ('Not enough input arguments!')
        os.sys.exit(0)

    inps = parser.parse_args(args=iargs)

    if int(inps.naz) == 1 and int(inps.nrg) == 1:
        print ('\nNo multilooking specified!')
        os.sys.exit(0)
    return inps


def plot_result(dst):
    """ Plot the multilooked result; should be removed """
    ds   = gdal.Open(dst)
    arr  = ds.ReadAsArray()
    if arr.ndim == 3:
        amp, phs = arr[0], arr[1]

    elif 'complex' in str(arr.dtype):
        amp  = np.where(arr.real == nodata, np.nan, arr.real)
        phs  = np.where(arr.imag == nodata, np.nan, arr.imag)

    else:
        amp, phs  =  arr, np.full_like(arr, np.nan)

    fig, (axe1, axe2) = plt.subplots(figsize=(9,9), nrows=2, sharex=True)
    im   = axe1.imshow(amp, interpolation='nearest', origin='lower',
                            cmap='binary')
    axe1.set_title('Amplitude')
    imshow_cbar(im, pad=0.035)

    im   = axe2.imshow(phs, interpolation='nearest', origin='lower',
                            cmap='hsv')
    axe2.set_title('Wrapped Phase')
    imshow_cbar(im, pad=0.035)
    fig.set_label(op.basename(dst))
    fig.suptitle('Width {} height {}'.format(*amp.shape), y=0.93)
    # fig.subplots_adjust(wspace=0.01, hspace=0.01) # doesnt work

    if inps.save or (inps.view and BACKEND == 'agg'):
        savefigs(op.dirname(src), figures=True, overwrite=True)

    embed(colors='neutral') if inps.interact else ''
    return dst


## theres really no need for this
def ml_mask(inps):
    """ Multilook mask by translating using already multilooked wrapped IFG

    Also works for other files (e.g., corrv) with v. small diffs from looks.py
    """
    ext   = re.search('\dalks.*rlks', inps.input_file)
    if ext is None:
        print ('Trying to multilook mask; did you pass multilooked ifg?')
        ext == ''
    else:
        ext   = f'.{ext.group()}'

    root, ext1 = op.splitext(inps.mask)
    dst        = f'{root}{ext}{ext1}'

    ds_wrap_ml = gdal.Open(inps.input_file) # for multilook size

    ds_mask_ml = gdal.Translate(dst, inps.mask, format='ENVI',
                    width=ds_wrap_ml.RasterXSize, height=ds_wrap_ml.RasterYSize,
                                                        resampleAlg='lanczos')
    del ds_mask_ml

    cmap        = mpl.colors.ListedColormap(['white', 'black'])
    arr_mask_ml = gdal.Open(dst).ReadAsArray()
    arr_mask_fr = gdal.Open(inps.mask).ReadAsArray()

    fig, (axe1, axe2) = plt.subplots(figsize=(9,9), nrows=2)

    im   = axe1.imshow(arr_mask_fr, interpolation='nearest', origin='lower',
                            cmap=cmap)
    axe2.set_title('Original Mask')
    imshow_cbar(im, pad=0.035)

    im   = axe2.imshow(arr_mask_ml, interpolation='nearest', origin='lower',
                            cmap=cmap)
    axe2.set_title('Multilooked Mask')
    imshow_cbar(im, pad=0.035)
    fig.set_label(op.basename(dst))
    # fig.subplots_adjust(wspace=0.01, hspace=0.01) # doesnt work

    if inps.save or (inps.view and BACKEND == 'agg'):
        savefigs(op.dirname(src), figures=True, overwrite=True)
    print ('Multilooked mask to:', dst)
    return


def get_wrapped_ifgs(wrapdir):
    """ Get the wrapped files in a directory """
    print ('Multilooking directory:', inps.input_file)
    srcs = []
    for f in os.listdir(inps.input_file):
        ifg  = re.match('20[0-9]*_20[0-9]*.int$', f)
        if ifg is not None:
            src      = op.join(inps.input_file, ifg.group())
            src_test = '{}_ISCE{}'.format(*op.splitext(src))
            ## check if an ISCE version already exists
            if op.exists(src_test):
                srcs.append(src_test)
            else:
                srcs.append(src)
    return sorted(srcs)


def main(inps):
    if op.isdir(inps.input_file):
        srcs = get_wrapped_ifgs(inps.input_file)

    else:
        srcs = [inps.input_file]

    for src in srcs:
        ds  = gdal.Open(src, gdal.GA_ReadOnly)
        nodata = ds.GetRasterBand(1).GetNoDataValue()
        driver = ds.GetDriver().GetDescription()
        ## build the files for looks.py if you dont give an ISCE file
        if driver != 'ISCE':
            path, ext = op.splitext(src)
            dst       = f'{path}_ISCE{ext}'
            ds1 = gdal.Translate(dst, src, format='ISCE', noData=nodata)
            ds2 = gdal.Translate(f'{dst}.vrt', dst, format='VRT', noData=nodata)
            print ('Wrote:', dst)
            print (f'Wrote: {dst}.vrt')
            ## make the new ISCE file the src for looks.py
            src       = dst
            ds1.FlushCache(); ds2.FlushCache()
            del ds1, ds2

        ## run looks.py
        args = argparse.Namespace(infile=src, outfile=None,
                                        azlooks=inps.naz, rglooks=inps.nrg)
        dst  = looks.main(args)
        del ds

        if ('mask' in src.lower() and not 'shadow' in src.lower()) or inps.mask:
            ds   = gdal.Open(dst)
            arr  = ds.ReadAsArray()
            arr1 = np.where((arr>0) & (arr<1), 1, arr)

            driver = op.splitext(dst)[1][1:].upper()
            ds1    = gdal.GetDriverByName(driver).CreateCopy(dst, ds)
            ds1.GetRasterBand(1).WriteArray(arr1)
            del ds1
            print (f'Forced {dst} to 0 & 1')


if __name__ == '__main__':
    inps = cmdLineParse()

    main(inps)

    plt.show() if inps.view and not BACKEND == 'agg' else ''
