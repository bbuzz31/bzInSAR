#!/usr/bin/env python
"""
Custom script to make the output from snaphu readable with GDAL
Modified from: unwrap_fringe.py
    Works correctly. Band 1 is (i think) amplitude.
"""

from VLM.bzFRInGE import *
import isce
import isceobj

def createParser():
    parser = argparse.ArgumentParser(description = 'Make metafiles for external unwrapped data',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('unwrapFile', help='input unwrapped file')
    parser.add_argument('-w', '--wrapFile', default=None,
                    help='A file readable by gdal with same size as unwrapFile')

    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 1:
        parser.print_help()
        print ('Not enough arguments!')
        os.sys.exit(0)
    inps = parser.parse_args(args=iargs)


    ## unw will be two bands; first band is coherence i think
    if inps.wrapFile is None or inps.unwrapFile.endswith('conncomp'):
        inps.band     = 1
    else:
        inps.band     = 2

    return inps


def write_xml(src, bands, path_wrap=None):
    path_wrap     = src if path_wrap is None else path_wrap
    length, width = getSize(path_wrap)

    if inps.unwrapFile.endswith('conncomp'):
        dataType  = 'BYTE'
        scheme    = 'BIP'
    else:
        dataType  = 'FLOAT'
        scheme    = 'BIL'

    img = isceobj.createImage()
    img.setFilename(src)
    img.setWidth(width)
    img.setLength(length)
    img.setAccessMode('READ')
    img.bands    = bands
    img.dataType = dataType
    img.scheme   = scheme
    img.renderHdr()
    img.renderVRT()

    bname = op.basename(src)
    print (f'Wrote .xml and .vrt for: {bname}')
    return None


def getSize(path):
    ds = gdal.Open(path, gdal.GA_ReadOnly)
    if ds is None:
        print ('\nGDAL could not open:', path)
        print ('Probably you need to pass a wrapped file for the size')
        os.sys.exit()
    length = ds.RasterYSize
    width  = ds.RasterXSize
    ds     = None
    print ('\nLength:', length, 'Width:', width)
    return length, width


if __name__ == '__main__':
    inps = cmdLineParse()
    write_xml(inps.unwrapFile, inps.band, inps.wrapFile)
