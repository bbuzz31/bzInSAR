#!/usr/bin/env python
import shutil
import h5py
from __init__ import *

def createParser():
    parser = argparse.ArgumentParser(description='Make an HDF5 with the necessary metadata for MintPy to read it.',
                                     epilog='Examples of use:\n\t makeH5.py maskTempCoh.h5 ../masks/ampMask.msk '\
                                            '\n\t makeh5.py velocity.h5 Vup.h5 -k velocity',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('data_file', type=str, help='Path with data to use')
    parser.add_argument('attrib_file', type=str, help='Path to MintPy readable h5 file')
    parser.add_argument('-k', '--key', dest='key', default=None, help='H5 dataset to replace. Must be given if more than one.')
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter an interactive session on completion')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 3:
        parser.print_help()
        os.sys.exit(1)

    inps = parser.parse_args(args=iargs)
    return inps


def configure(cmd):
    iargs = cmd.split()[1:]
    inps  = cmdLineParse(iargs)
    main(inps)
    return


def main(inps):
    wd  = op.dirname(inps.attrib_file) # presumably, this is mintpy working dir
    dst = op.join(wd, f'{op.splitext(op.basename(inps.data_file))[0]}.h5')

    ## copy the mintpy object
    shutil.copy(inps.attrib_file, dst)

    if op.splitext(inps.data_file) == '.h5':
        with h5py.File(dst, 'r') as h5_orig:
            if inps.key is None:
                keys = list(h5_new.keys())
                assert len(keys) == 1, 'Must pass a dataset key since there are multiple options'
            arr = h5_orig[keys][:]

    else:
        arr = gdal.Open(inps.data_file).ReadAsArray()


    with h5py.File(dst, 'r+') as h5_new:
        if inps.key is None:
            keys = list(h5_new.keys())
            assert len(keys) == 1, 'Must pass a dataset key since there are multiple options'

        data    = h5_new[keys[0]]
        data[:] = arr
        h5_new.attrs['FILE_PATH'] = op.join(dst)

    # you can compare them with view.py attrib_file, dst
    print ('Wrote:', dst)
    embed(colors='neutral') if inps.interact else ''
    return


if __name__ == '__main__':
    inps = cmdLineParse()
    main(inps)
