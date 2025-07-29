#!/usr/bin/env python

import argparse
import h5py
from BZ import *

def createParser():
    parser = argparse.ArgumentParser(description='Compare two mintpy h5 timeseries',
                                     epilog='Examples of use:\n\t read_hdf5.py timeseries.h5 timeseries_SR.h5',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inp1', type=str, help='Path to h5 1')
    parser.add_argument('inp2', type=str, help='Path to h5 2')
    parser.add_argument('-l', '--layer', type=str, default='timeseries', help='Layer (timeseries)')

    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 3:
        parser.print_help()
        os.sys.exit(1)

    inps = parser.parse_args(args=iargs)
    return inps


def main(inps):
    with h5py.File(inps.inp1, 'r') as h5:
        arr1 = h5[inps.layer][:]

    with h5py.File(inps.inp2, 'r') as h5:
        arr2 = h5[inps.layer][:]

    print ('Allclose?', np.allclose(arr1, arr2, equal_nan=True))
    print ('Max Resid:', np.nanmax(np.abs(arr1-arr2)))

    try:
        embed(colors='neutral')
    except:
        breakpoint()

    return

if __name__ == '__main__':
    inps = cmdLineParse()
    main(inps)
