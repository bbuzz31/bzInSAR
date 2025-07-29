#!/usr/bin/env python
from __init__ import *

def createParser():
    parser = argparse.ArgumentParser(description='Make a mask based on the magnitude of the SLC stack. Call from parent of coreg_stack.',
                                     epilog='Examples of use:\n\t makeAmpMask.py -t 100\n\t makeAmpMask.py -w ./StudyArea',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-w', '--workdir', dest='wd', default='.', type=str, help='Path to were "coreg_stack" directory is.')
    parser.add_argument('--ifg', dest='ifg', default=None, help='Path to use a single ifg for masking (ie, a multilooked one or coherence).')
    parser.add_argument('-t', '--threshold', dest='thresh', type=float, default=None, help='Mask magnitude values below this threshold.')
    parser.add_argument('-n', '--nifgs', dest='nifgs', type=float, default=30, help='Number of SLCs to use in computation (Default=30).')
    parser.add_argument('-v', '--view', action='store_true', help='View the mask')
    parser.add_argument('-s', '--save', action='store_true', help='Save to maskAmp.png')
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter an interactive session on completion')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    inps   = parser.parse_args(args=iargs)
    return inps


def main(inps):
    stack = op.join(inps.wd, 'coreg_stack', 'slcs_base.vrt')
    ds    = gdal.Open(stack, gdal.GA_ReadOnly)
    # magnitude works MUCH better than amp
    arr   = ds.ReadAsArray()
    nifgs = arr.shape[0]
    # take 30 or less ifgs; otherwise memory problems
    if nifgs < inps.nifgs:
        arr = np.abs(arr)
    else:
        np.random.seed(10)
        idx = sorted(np.random.randint(0, nifgs, inps.nifgs))
        arr = np.abs(arr[idx])

    amp   = arr.mean(0)

    thresh = amp.mean() if inps.thresh is None else inps.thresh
    mask   = np.where(amp>=thresh, 1, 0)

    print (f'Masked {100*(mask.size - mask.sum())/mask.size:.2f}% of pixels '\
           f' below threshold: {thresh}')


    # write the mask
    path_masks = op.join(inps.wd, 'masks')
    os.makedirs(path_masks, exist_ok=True)

    basen  = 'ampMask.msk'
    dst    = op.join(path_masks, basen)
    ## can't do this because its now float and was complex
    # alternatively, may want to directly add mask band to slcs_base.vrt
    ds_out = gdal.GetDriverByName('ENVI').Create(dst, ds.RasterXSize,
                                        ds.RasterYSize, 1, gdalconst.GDT_Byte)
    ds_out.GetRasterBand(1).WriteArray(mask)
    print (f'Wrote: {dst}')

    ## quick plot
    cmap = mpl.colors.ListedColormap(['white', 'black'])
    fig, axes = plt.subplots(figsize=(9,9))
    im   = axes.imshow(mask, interpolation='nearest', origin='lower',
                            cmap=cmap)
    ticks = np.array([0, 1])
    imshow_cbar(im, **dict(ticks=ticks, pad=0.25)).set_ticks(ticks)
    fig.set_label(basen)

    if inps.save or (inps.view and BACKEND == 'agg'):
        savefigs(path_masks, figures=False, overwrite=True)

    embed(colors='neutral') if inps.interact else ''

    return mask


## TODO
def make_ml_mask(inps):
    """ Make a multilooked mask using just one ifg """
    raise Exception('Masking from multilook not yet working')
    ds    = gdal.Open(inps.ifg, gdal.GA_ReadOnly)
    # magnitude works MUCH better than amp
    arr   = np.abs(ds.ReadAsArray())
    amp   = arr.mean(0) if arr.ndim == 3 else arr

    thresh = amp.mean() if inps.thresh is None else inps.thresh
    mask   = np.where(amp>=thresh, 1, 0)

    print (f'With threshold of {thresh:.3f},')
    print (f'Masked {100*(mask.size - mask.sum())/mask.size:.2f}% of pixels')

    # write the mask
    path_masks = op.join(inps.wd, 'masks')
    os.makedirs(path_masks, exist_ok=True)


    ## add the extension to the filename
    ext   = re.search('\dalks.*rlks', inps.ifg)
    ext   = f'_{ext.group()}' if ext is not None else ''
    basen = f'ampMask{ext}.msk'
    dst   = op.join(path_masks, basen)

    ## can't do a copy because its now float and was complex
    # alternatively, may want to directly add mask band to slcs_base.vrt
    ds_out = gdal.GetDriverByName('ENVI').Create(dst, ds.RasterXSize,
                                        ds.RasterYSize, 1, gdalconst.GDT_Byte)
    ds_out.GetRasterBand(1).WriteArray(mask)
    print (f'Wrote: {dst}')

    ## quick plot
    cmap = mpl.colors.ListedColormap(['white', 'black'])
    fig, axes = plt.subplots(figsize=(9,9))
    im   = axes.imshow(mask, interpolation='nearest', origin='lower',
                            cmap=cmap)
    ticks = np.array([0, 1])
    imshow_cbar(im, **dict(ticks=ticks, pad=0.25)).set_ticks(ticks)
    axes.set_title(basen)
    fig.set_label(basen)

    if inps.save or (inps.view and BACKEND == 'agg'):
        savefigs(path_masks, figures=False, overwrite=True)

    embed(colors='neutral') if inps.interact else ''

    return mask

## maybe it will work better when ts is longer
def make_corr_mask(inps):
    """ Make a multilooked mask using correlation file """
    ds    = gdal.Open(inps.ifg, gdal.GA_ReadOnly)
    arr   = ds.ReadAsArray()
    arr   = arr.mean(0) if arr.ndim == 3 else arr

    thresh = 0.9 if inps.thresh is None else inps.thresh
    mask   = np.where(arr>=thresh, 1, 0)

    print (f'With threshold of {thresh:.3f},')
    print (f'Masked {100*(mask.size - mask.sum())/mask.size:.2f}% of pixels')

    # write the mask
    path_masks = op.join(inps.wd, 'masks')
    os.makedirs(path_masks, exist_ok=True)


    ## add the extension to the filename
    ext   = re.search('\dalks.*rlks', inps.ifg)
    ext   = f'_{ext.group()}' if ext is not None else ''
    basen = f'corrMask{ext}.msk'
    dst   = op.join(path_masks, basen)

    ## can't do a copy because its now float and was complex
    # alternatively, may want to directly add mask band to slcs_base.vrt
    ds_out = gdal.GetDriverByName('ENVI').Create(dst, ds.RasterXSize,
                                        ds.RasterYSize, 1, gdalconst.GDT_Byte)
    ds_out.GetRasterBand(1).WriteArray(mask)
    print (f'Wrote: {dst}')

    ## quick plot
    cmap = mpl.colors.ListedColormap(['white', 'black'])
    fig, axes = plt.subplots(figsize=(9,9))
    im   = axes.imshow(mask, interpolation='nearest', origin='lower',
                            cmap=cmap)
    ticks = np.array([0, 1])
    imshow_cbar(im, **dict(ticks=ticks, pad=0.25)).set_ticks(ticks)
    axes.set_title(basen)
    fig.set_label(basen)

    if inps.save or (inps.view and BACKEND == 'agg'):
        savefigs(path_masks, figures=False, overwrite=True)

    embed(colors='neutral') if inps.interact else ''

    return mask


if __name__ == '__main__':
    inps = cmdLineParse()
    # main(inps)
    # make_ml_mask(inps)
    make_corr_mask(inps)


    plt.show() if inps.view and not BACKEND == 'agg' else ''
