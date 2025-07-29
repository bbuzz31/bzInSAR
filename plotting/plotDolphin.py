#!/usr/bin/env python
import cmocean
from collections import OrderedDict
from VLM.bzFRInGE import *
from BZ import bbPlot
import rioxarray as xrr
import pyproj
import cartopy.crs as ccrs
from cartopy.io import img_tiles as cimgt


gdal.UseExceptions()

## example on how to 0 out the amplitude in an interferogram
# imageMath.py -e '0+imag(a)*J' -t cfloat --a 20190124_20190217_ISCE.7alks_19rlks0.int -o test.int

def createParser():
    parser = argparse.ArgumentParser(description='Program to view/save with GDAL some datasets produced by FRINGE.',
                                     epilog='Examples of use:\n\t plotFringe.py ./KS2/count \n\t plotFringe.py ./EVD/20190124.slc',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', metavar='input_file', type=Path, help='Path to file to plot')
    parser.add_argument('-m', '--mask', dest='path_mask', default=None, help='Path to mask of same size to apply')
    parser.add_argument('-v', '--view', default=True, type=bool, help='View the plot')
    parser.add_argument('--vmin', default=None)
    parser.add_argument('--vmax', default=None)
    parser.add_argument('-s', '--save', action='store_true', help='Save to src.png')
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter an interactive session on completion')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(0)

    inps = parser.parse_args(args=iargs)
    return inps


def configure(cmd, figsize=(6,6), labelsize=10):
    iargs = cmd.split()[1:]
    inps  = cmdLineParse(iargs)

    if op.isdir(inps.input):
        f, a = plot_isce(inps, figsize, labelsize)

    else:
        f, a = main(inps, figsize, labelsize)

    return f, a


def default_parms(vmin, vmax):
    vmin = float(vmin) if isinstance(vmin, str) else vmin
    vmax = float(vmax) if isinstance(vmax, str) else vmax
    norm = mpl.colors.TwoSlopeNorm(0, vmin, vmax)
    dct_parms = {
                'cmap': 'binary',    # colorbar
                'norm': norm
            }

    return dct_parms

## todo
def default_parms_cbar():
    dct_parms = {'':''}
    return dct_parms


def main(inps, figsize=(11, 11), ls=15):
    """ figsize / labelsize for plotting in jupyter """
    path     = inps.input.resolve()
    assert path.exists(), f'{path} doesnt exist'
    basen    = op.basename(path).replace('.vrt', '')
    ti, ext  = op.splitext(basen)

    ds   = xrr.open_rasterio(path)
    band = 2 if ext == '.unw' else 1 # untested

    da = ds.sel(band=band)
    da = da.rio.reproject('EPSG:4326') # the warning doesnt matter
    da = da.where(da != da.rio.nodata, np.nan)

    DCT_PARMS  = default_parms(inps.vmin, inps.vmax)
    clim       = None

    if inps.path_mask is not None:
        mask = xrr.open_rasterio(inps.path_mask).ReadAsArray()
        print ('CHECK THE CRS ON THE MASK')
        assert len(np.unique(mask)) == 2, 'Mask must be only 0 (water) and 1 (land)'
        da  *= xr.where(mask==0, np.nan, 1)

    # if nodata is not None:
    #     nans = np.nan+np.nan*1j if 'CFloat' in gdal_dtype else np.nan
    #     arr  = np.where(np.abs(arr) == nodata, nans, arr)

    if ext == '.slc':
        da   = np.abs(da)
        clab = 'magnitude'

    # something wrong with topophase.flat; its scaled or something
    # elif ext in '.int .flat'.split():
    elif ext == '.int':
        if da.dtype in [np.float32, np.float64]:
            clab = 'amplitude'
        else:
            clab    = 'phase'
            da.data = np.angle(da)
            # DCT_PARMS['cmap'] = 'cmo.phase'
            DCT_PARMS['cmap'] = 'turbo'

    elif ext == '.phs' or basen.endswith('phase'):
        clab = 'phase'
        DCT_PARMS['cmap'] = 'cmo.phase'

    elif ext == '.unw' or ti == 'unw':
        clab = 'unwrapped phase'
        # arr  = np.where(np.abs(arr)>999, np.nan, arr)
        DCT_PARMS['cmap'] = 'coolwarm'

    elif ext == '.conncomp':
        clab = 'Conn Comp'
        # cc    = np.unique(np.ma.masked_invalid(arr)).compressed()
        cc    = np.unique(da)
        cc    = cc[~np.isnan(cc)]
        ticks = list(range(len(cc)+1))
        # clim  = (-0.5, len(cc)-1.5)
        cmap = plt.get_cmap('Dark2', len(cc))
        DCT_PARMS['cmap'] = cmap
        # DCT_PARMS.pop('norm', None)
        DCT_PARMS['norm'] = mpl.colors.BoundaryNorm(ticks, cmap.N, clip=True)

    elif ext in '.msk .wbd'.split() or 'swbd' in ti or 'waterMask' in ti or 'OSM' in ti:
        clab = ''
        # in case there onlys 1 and nans
        if len(np.unique(da)) == 2:
            cmap = mpl.colors.ListedColormap(['white', 'black'])
            DCT_PARMS['cmap'] = cmap
            DCT_PARMS['norm'] = mpl.colors.BoundaryNorm([0, 1, 2], cmap.N, clip=True)

        else:
            DCT_PARMS['cmap'] = plt.cm.binary

    elif ti.startswith('temporal_c') or ext == '.cor':
        clab = 'correlation'
        cmap = plt.cm.magma
        # cmap   = mpl.colors.ListedColormap(colors)
        DCT_PARMS['cmap'] = cmap
        DCT_PARMS['norm'] = mpl.colors.BoundaryNorm(np.linspace(0, 1, 11),
                                                    cmap.N, clip=True)
        # otherwise multilooked images gets weird

    elif ti == 'ampdispersion':
        clab = 'amplitude'

    elif ti.startswith('ps'):
        clab = 'ps'
        cmap = mpl.colors.ListedColormap(['white', 'black'])
        DCT_PARMS['cmap'] = cmap
        DCT_PARMS['norm'] = mpl.colors.BoundaryNorm([0, 1, 2], cmap.N, clip=True)

    elif ti.startswith('hgt') or ti.startswith('hei') or ext == '.dem':
        clab              = 'height'
        DCT_PARMS['cmap'] = 'coolwarm'

    else:
        raise Exception(f'ti={ti} with ext={ext} is not yet a supported filetype')


    fig, axes = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
    im   = axes.pcolormesh(da.x, da.y, da, shading='auto',  **DCT_PARMS)
    im.set_clim(clim) if clim is not None else ''
    cbar = bbPlot.cartopy_cbar(im, xlabel=clab)
    cbar.ax.set_yscale('linear')
    cbar.ax.tick_params(labelsize=ls/2)
    axes.tick_params(axis='both', which='major', labelsize=ls/2)
    axes.set_title(basen, fontsize=ls)
    gl = axes.gridlines(draw_labels=True)
    bbPlot.fmt_gridlines(gl, bottom=True)



    fig.set_label(basen)

    if inps.save:# or (inps.view and BACKEND == 'agg'):
        savefigs(op.dirname(path), figures=True, overwrite=True)
        # ds_new = geocode(ds, arr)
        # dst = op.join(op.dirname(path), title.replace(' ', '_'))
        # ds_tif = gdal.Translate('{}.tif'.format(dst), ds_new)
        # ds_vrt = gdal.BuildVRT('{}.vrt'.format(dst), ds_tif)
        # ds_new = ds_tif = ds_vrt = None

    embed(colors='neutral') if inps.interact else ''
    return fig, axes


def plot_isce(inps, figsize=(11, 11), ls=15):
    """ Make a four paneled plot of the ISCE results """
    tis     = ['Correlation', 'Wrapped Phase', 'Unwrapped Phase', 'Connected Components']
    dct_res = OrderedDict()
    n       = 0
    for ti, ext in zip(tis, 'cor int unw unw.conncomp'.split()):
        path = op.join(inps.input, f'filt_fine.{ext}')
        norm = None

        if op.exists(path):
            ds   = gdal.Open(path)
            band = ds.GetRasterBand(2) if ext == 'unw' else ds.GetRasterBand(1)
            arr  = np.flipud(band.ReadAsArray()) if inps.flip else band.ReadAsArray()

            if ext == 'cor':
                cmap = plt.cm.magma
                norm =  mpl.colors.BoundaryNorm(np.linspace(0, 1, 11), cmap.N, clip=True)


            elif ext == 'int':
                if arr.dtype in [np.float32, np.float64]:
                    clab = 'amplitude'
                else:
                    arr  = np.angle(arr)
                    cmap = 'cmo.phase'

            elif ext == 'unw':
                cmap = 'coolwarm'

            elif ext == 'unw.conncomp':
                cmap = plt.cm.tab20c
                norm =  mpl.colors.BoundaryNorm(np.linspace(0, 20, 21), cmap.N, clip=True)


            dct_res[ti] = [arr, cmap, norm]
            n += 1

        else:
            print (path, 'does not exist, skipping...')

    nrows, ncols = (2, 2) if n == 4 else (n, 1)
    fig, axes    = plt.subplots(figsize=figsize, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True)
    axe = axes.ravel()
    for i, (ti, lst) in enumerate(dct_res.items()):
        im = axe[i].imshow(lst[0], cmap=lst[1], norm=lst[2], interpolation='nearest')
        divider = make_axes_locatable(axe[i])
        cbar_ax = divider.append_axes('right', size='5%', pad=0.1)
        cbar    = plt.colorbar(im, cax=cbar_ax)
        axe[i].set_title(ti, fontsize=ls)


    fig.subplots_adjust(hspace=0.01, wspace=0.2)
    fig.suptitle(op.basename(inps.input), y=0.825, fontsize=ls)
    fig.tight_layout()
    return fig, axes


if __name__ == '__main__':
    inps = cmdLineParse()

    if op.isdir(inps.input):
        plot_isce(inps)

    else:
        main(inps)

    plt.show() if inps.view and not BACKEND == 'agg' else ''
