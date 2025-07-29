#!/usr/bin/env python
import shutil
from VLM.bzFRInGE import *
import multilookWrapped
"""
This script moves around the geometry files and multilooks to make everything happy for mintpy
"""

def createParser():
    parser = argparse.ArgumentParser(description='Setup the multilook geometry for MintPy',
                                     epilog='Examples of use:\n\t setup_multilook.py ./HR --naz 2 --nrg 5',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path_root', type=str, help='Path to home of geometry files')
    parser.add_argument('--path_ml', default=None, type=str,help= 'Path to store multilooked'\
                                 'geometry files. Default: path_root/geometry_multilook')
    parser.add_argument('--naz', default=7, help='Number of azimuth looks (default=7)')
    parser.add_argument('--nrg', default=19, help='Number of range looks (default=19)')
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter an interactive session on completion')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(1)

    inps = parser.parse_args(args=iargs)

    if int(inps.naz) == 1 and int(inps.nrg) == 1:
        print ('\nNo multilooking specified!')
        os.sys.exit(0)

    return inps


def make_fake_vrts(src):
    for k in 'shadowMask incLocal los'.split():
        bands = 1 if k == 'shadowMask' else 2
        dst   = op.join(op.dirname(src), f'{k}.vrt')
        with open(dst, 'w') as fhw:
            for band in range(1, bands+1):
                with open(src, 'r') as fh:
                    for line in fh:
                        if 'Filename' in line:
                            ## replace the filename prefix
                            f0   = op.splitext(op.basename(src))[0]
                            line = line.replace(f0, k)
                        if 'band' in line.lower():
                            line = line.replace('1', str(band))
                        if '/VRTDataset' in line and band < bands:
                            line = ''
                        if 'VRTDataset ' in line and band > 1:
                            line = ''
                        fhw.write(line)


def rename_multilook_geom(vrts):
    """ rename the full resolution vrts to full and multilook """
    for src in vrts:
        dst = f'{op.splitext(src)[0]}_full'
        os.rename(src, dst)
        inps_ml = argparse.Namespace(input_file=dst, naz=inps.naz,
                                            nrg=inps.nrg, mask=False)
        multilookWrapped.main(inps_ml)

        # correct its name (had to be left off previously)
        os.rename(dst, f'{dst}.vrt')
    return


def symlink_full(vrts):
    """ symlink the multilooked files to 'full', which actually get read """
    for vrt in vrts:
        root = op.splitext(vrt)[0]
        # src  = f'{root}.{inps.naz}alks_{inps.nrg}rlks'
        src  = f'{op.basename(root)}.{inps.naz}alks_{inps.nrg}rlks'
        dst  = op.join(op.dirname(root), f"{op.basename(root).split('_')[0]}.rdr.full")
        os.symlink(src, dst)
        os.symlink(f'{src}.vrt', f'{dst}.vrt')
    return


def multilook_ps_pixels(path_ps):
    inps_ml  = argparse.Namespace(input_file=path_ps, naz=inps.naz,
                                        nrg=inps.nrg, mask=False)
    multilookWrapped.main(inps_ml)
    return


def add_src_rect(path_geom_ml):
    """ Add the src rect field to path_lat_vrt for mintpy to read """
    path_lat = op.join(path_geom_ml, 'lat.vrt')
    # first get the sizes to write
    ds       = gdal.Open(path_lat)
    xs, ys   = ds.RasterXSize, ds.RasterYSize
    del ds
    line0    = f'        <SrcRect xOff="0" yOff="0" xSize="{xs}" ySize="{ys}"/>\n'
    line1    = line0.replace('Src', 'Dst')


    # now copy the original vrt to a temporary one
    path_tmp = op.join(path_geom_ml, 'lat_orig.vrt')
    os.rename(path_lat, path_tmp)

    # now write out lat.vrt with the rect
    with open(path_lat, 'w') as fhw:
        with open(path_tmp, 'r') as fh:
            for line in fh:
                fhw.write(line)
                if 'LineOffset' in line:
                    fhw.write(line0)
                    fhw.write(line1)
    # remove the original
    os.remove(path_tmp)
    return


def main(inps):
    path_geom = op.join(inps.path_root, 'geometry')
    ml        = f'{inps.naz}alks_{inps.nrg}rlks'
    assert op.exists(path_geom), f'Cannot find the geometry path: {path_geom}'

    # make the destination
    path_geom_ml = op.join(path_geom, f'geometry_{ml}')
    os.makedirs(path_geom_ml, exist_ok=True)

    geom_files   = os.listdir(path_geom_ml)
    [os.remove(op.join(path_geom_ml, f)) for f in geom_files]

    # ugly hack to create the shadowmask/los/inc for mintpy
    make_fake_vrts(op.join(path_geom, 'lat.vrt'))

    # copy lat/lon/hgt vrts to the multilooked dir
    vrts = glob.glob(op.join(path_geom, '*.vrt'))

    for vrt in vrts:
        dst = op.join(path_geom_ml, op.basename(vrt))
        shutil.copy(vrt, dst)


    # handle the geometry files
    vrts = glob.glob(op.join(path_geom_ml, '*.vrt'))
    rename_multilook_geom(vrts)
    vrts = glob.glob(op.join(path_geom_ml, '*ISCE.vrt'))
    symlink_full(vrts)

    # not doing here because I dont want to pass nx ny
    # multilook the ps_pixels as well
    # path_ps  = op.join(inps.path_root, 'ampDispersion', 'ps_pixels')
    # multilook_ps_pixels(path_ps)


    ## make the original lat.vrt which gets read for properties by prep_fringe
    src = op.join(path_geom_ml, 'lat.rdr.full.vrt')
    dst = op.join(path_geom_ml, 'lat.vrt')
    shutil.copy(src, dst)

    ## manually edit the src rect in lat.vrt
    add_src_rect(path_geom_ml)

    print ('\nWrote multilooked geometry files to:', path_geom_ml)


if __name__ == '__main__':
    inps = cmdLineParse()
    main(inps)
