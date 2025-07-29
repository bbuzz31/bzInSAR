#!/usr/bin/env python
import h5py
from VLM.bzFRInGE import *

def createParser():
    parser = argparse.ArgumentParser(description='Pull a datalayer out of an hdf5 and write it to gdal',
                                     epilog='Examples of use:\n\t hdf52gdal.py mask.h5 -of ENVI'\
                                            '\n\t makeh5.py velocity.h5 Vup.h5 -k velocity',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', type=str, help='Path to file to convert')
    parser.add_argument('-l', '--layer', type=str, default=None, help='Layer in the h5 file to convert')
    parser.add_argument('--of', dest='fmt', default='GTiff', help='GDAL writeable format')
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter an interactive session on completion')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        print ('Not enough input arguments.')
        os.sys.exit(0)

    inps = parser.parse_args(args=iargs)
    return inps


def configure(cmd):
    iargs = cmd.split()[1:]
    inps  = cmdLineParse(iargs)
    main(inps)
    return


def main(inps):

    with h5py.File(inps.input, 'r') as h5:
        if inps.layer is None:
            keys = list(h5.keys())
            assert len(keys) == 1, 'Must pass a dataset key since there are multiple options'
        arr = h5[keys[0]][:]

    if 'swbd' in op.basename(inps.input):
        print ('I think this is an swbd mask file. Will write out 0 (water) and 1 (land)')
        arr = np.where(arr < 0, 0, 1)

    wd  = op.dirname(inps.input) # presumably, this is mintpy working dir
    dst = op.join(wd, f'{op.splitext(op.basename(inps.input))[0]}.{inps.fmt}')

    ds_out = gdal.GetDriverByName(inps.fmt).Create(dst, arr.shape[1],
                                        arr.shape[0], 1, gdalconst.GDT_Float32)
    ds_out.GetRasterBand(1).WriteArray(arr)
    del ds_out

    print ('Wrote:', dst)
    embed(colors='neutral') if inps.interact else ''
    return


if __name__ == '__main__':
    inps = cmdLineParse()
    main(inps)
