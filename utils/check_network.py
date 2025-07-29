#!/usr/bin/env python3

import argparse
import os, sys
import time
import h5py
import numpy as np
from scipy import linalg
from mintpy.objects import ifgramStack
from mintpy.utils import readfile


################################################################################################
EXAMPLE = """example:
  check_inversion.py inputs/ifgramStack.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Check whether network of interferograms is connected.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    # input dataset
    parser.add_argument('ifgramStackFile', help='interferograms stack file to be inverted')

    parser.add_argument('-i','-d', '--dset', dest='obsDatasetName', type=str,
                        help='dataset name of unwrap phase / offset to be used for inversion'
                             '\ne.g.: unwrapPhase, unwrapPhase_bridging, ...')

    parser.add_argument('--ref-date_idx', '-ri', dest='refi', default=0, help='Reference date index, first date by default.')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input file type
    atr = readfile.read_attribute(inps.ifgramStackFile)
    if atr['FILE_TYPE'] not in ['ifgramStack']:
        raise ValueError('input is {} file, support ifgramStack file only.'.format(atr['FILE_TYPE']))

    return inps


def main(path_ifgramStack, refi=0):
    """ Check if network is connected """
    refi      = int(refi)
    stack_obj = ifgramStack(path_ifgramStack)
    stack_obj.open(print_msg=False)

    ## mostly for removing the first few
    date_list = stack_obj.get_date_list(dropIfgram=True)
    ref_date  = date_list[refi]
    date_list = date_list[refi:]

    date12_list0 = stack_obj.get_date12_list(dropIfgram=True)
    date12_list  = []

    for dt in date12_list0:
        st, en = dt.split('_')
        if (st in date_list) and (en in date_list):
            date12_list.append(dt)

    length, width = stack_obj.length, stack_obj.width

    # 1.2 design matrix
    A = stack_obj.get_design_matrix4timeseries(date12_list)[0]
    num_pair, num_date = A.shape[0], A.shape[1]+1

    print(f'number of acquisitions  : {num_date}')
    print(f'number of interferograms: {num_pair}')
    print(f'number of lines   : {length}')
    print(f'number of columns : {width}')
    print(f'Reference date : {ref_date}')


    # 1.3 print key setup info
    msg = '-------------------------------------------------------------------------------\n'
    if np.linalg.matrix_rank(A) < A.shape[1]:
        msg += '***WARNING: the network is NOT fully connected.\n'
        msg += '\tInversion result can be biased!\n'
        msg += '\tContinue to use SVD to resolve the offset between different subsets.\n'
    else:
        msg += '***GOOD: The network is fully connected.\n'

    msg += '-------------------------------------------------------------------------------'

    print(msg)

    return


################################################################################################
if __name__ == '__main__':
    parser = create_parser()
    inps   = parser.parse_args()
    main(inps.ifgramStackFile, inps.refi)
