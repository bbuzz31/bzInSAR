#!/usr/bin/env python
# Author: Marin Govorcin
# Caltech NASA/JPL, Dec 05, 2021

import os
import argparse
from isce3.io import Raster
from isce3.unwrap import Phass
from osgeo import gdal

def cmdLineParser():
    '''
    Command Line Parser
    Script to use Phass in the isce3 environment
    Note, if you get Image not found error, make sure to run 'python' unwrap_phass_MG.py ...
    TODO: include the option to use input Power file
    '''
    parser = argparse.ArgumentParser(description='Phass Unwrapper', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--interferogram_file',type=str, dest='phase',
                        required=True,help='Input interferogram file in radians, example: file.int')
    parser.add_argument('-c','--correlation_file',type=str,dest='correlation',
                        required=True,help='Input correlation file, example: tcorr_ds_ps.bin')
    parser.add_argument('-o','--output_dir',type=str,dest='outputDir',
                        default='./',help='Output directory path')
    parser.add_argument('-ct','--correlation_threshold',type=float,dest='correlation_threshold',
                        default=0.2,help='Phass parameter - Correlation threshold')
    parser.add_argument('-gc','--good_correlation',type=float,dest='good_correlation',
                        default=0.7,help='Phass parameter - Good Correlation,')
    parser.add_argument('-mp','--minPixels',type=int,dest='min_pixels_region',
                        default=200,help='Phass parameter - Minimum pixels per region')

    return parser.parse_args()

######################################################################

def cpx2rad(infile):
    '''
    Function to extract interferogram phase from complex file
    needed when unwrapping interferogram in radians
    '''

    ## Find the input interferogram dtype, if complex32 or complex64
    ## create vrt file that extract phase values onthefly

    ds     = gdal.Open(infile+'.vrt',gdal.GA_ReadOnly)
    ## BB - dont use vrt
    # ds     = gdal.Open(infile, gdal.GA_ReadOnly)
    width  = ds.RasterXSize
    height = ds.RasterYSize
    stype  = ds.GetRasterBand(1).DataType

    # Get path to directory of inputfile and inputfile name without extension
    inDir  = os.path.dirname(os.path.abspath(infile))
    inName = os.path.splitext(os.path.basename(infile))[0]

    GDAL_DATATYPE = {
         6 : "Float32",
         7 : "Float64",
         10 : "CFloat32",
         11 : "CFloat64"
    }

    #Input file is COMPLEX, create VRT with default pixelfunction : phase to extract phase in radian
    if stype == 10 or stype == 11:

        vrttmpl='''<VRTDataset rasterXSize="{width}" rasterYSize="{height}">
    <VRTRasterBand dataType="{DTYPE}" band="1" subClass="VRTDerivedRasterBand">
        <Description>Wrapped interferogram in radians on-the-fly</Description>
        <SimpleSource>
           <SourceFilename relativeToVRT="1">{PATH}</SourceFilename>
           <ImageOffset>0</ImageOffset>
           <PixelOffset>4</PixelOffset>
           <LineOffset>{linewidth}</LineOffset>
           <ByteOrder>LSB</ByteOrder>
        </SimpleSource>
        <PixelFunctionType>phase</PixelFunctionType>
        <SourceTransferType>{STYPE}</SourceTransferType>
    </VRTRasterBand>
</VRTDataset>'''

        if stype == 10:
            dtype = GDAL_DATATYPE[6]
        elif stype == 11:
            dtype = GDAL_DATATYPE[7]

        outfile = os.path.join(inDir, '{0}.phase.vrt'.format(inName))

        with open(outfile, 'w') as fid:
            fid.write( vrttmpl.format(width=width,
                                     height=height,
                                     DTYPE=dtype,
                                     STYPE=GDAL_DATATYPE[stype],
                                     PATH=infile,
                                     linewidth=4*width))
        new_infile  = outfile

    else:
        new_infile  = os.path.abspath(infile)

    # Return the path to VRT file as input to Phass
    return new_infile

#def write_xml

#########################################################
# obsolete readFile and plotUnw function, not in use
def readFile(infile,band):
    print('Read file: {0}'.format(infile))
    ds     = gdal.Open(infile,gdal.GA_ReadOnly)
    data   = ds.GetRasterBand(band).ReadAsArray()
    #Map extent
    trans  = ds.GetGeoTransform()
    xsize  = ds.RasterXSize
    ysize  = ds.RasterYSize
    extent = [trans[0], trans[0] + xsize * trans[1],
          trans[3] + ysize*trans[5], trans[3]]
    ds     = None

    return data,extent

def plotUnw(phase,correlation,unw_igram,label):
    # Function to plot wrapped and unwrapped interferogram,
    #          correlation and connected components
    from matplotlib import pyplot as plt
    import numpy as np

    plt.figure('Phass Unwrapping')
    plt.subplot(1,4,1,aspect='equal')
    data,extent1 = readFile(phase,1)
    plt.imshow(np.angle(data), clim=[-np.pi, np.pi], extent=extent1, cmap='jet')
    plt.title('Wrapped phase')

    plt.subplot(1,4,2,aspect='equal')
    data,extent1 = readFile(correlation,1)
    plt.imshow(data, clim=[0,1],extent=extent1,cmap='gray')
    plt.title('Correlation')

    plt.subplot(1,4,3,aspect='equal')
    data, extent1 = readFile(unw_igram,1)
    plt.imshow(data, clim=[-5,5],extent=extent1,cmap='Spectral')
    plt.title('Unwrapped phase')

    plt.subplot(1,4,4,aspect='equal')
    data,extent1 = readFile(label,1)
    plt.imshow(data,clim=[-5,5],extent=extent1,cmap='gray')
    plt.title('Connected Components')


    plt.show()

####################################################################################
####################################################################################
if __name__ == '__main__':
    '''
    Main driver
    '''

    inps = cmdLineParser()


    # Get the input interferogram name

    inpName = os.path.splitext(os.path.basename(inps.phase))[0]

    ######################### Phass input and output files ########################
    '''
    Input:
        phase       :: isce3.io.Raster,  Wrapped phase in radians [interferogram *.int in PS_DS folder for fringe]
        correlation :: isce3.io.Raster,  Correlation file [tcorr_ds_ps.bin for fringe]
    Output:
        unw_igram   :: isce3.io.Raster, Unwrapped interferogram, type: isce3.io.Raster
        label       :: isce3.io.Raster, Connected component labels, type: isce3.io.Raster
    '''

    # Check if phase is COMPLEX format
    Phase       = cpx2rad(inps.phase)

    phase       = Raster(Phase)
    correlation = Raster(inps.correlation+'.vrt')
    unw_igram   = Raster(os.path.join(inps.outputDir,inpName+'.unw'),phase.width,
                        phase.length, 1, gdal.GDT_Float32, "GeoTiff") #works also with ISCE format
    label       = Raster(os.path.join(inps.outputDir,inpName+'_int.conn'), phase.width,
                        phase.length, 1, gdal.GDT_Byte, "GeoTiff")

    print('\nPhass Unwrapping:')
    print('Input:')
    print('  Phase: ' + Phase)
    print('  Correlation: ' + inps.correlation+'.vrt')

    ######################## Set-up Phass object and Unwrap ############################

    phass = Phass()

    #Set parameters for Phass unwrapping
    phass.correlation_threshold = inps.correlation_threshold
    phass.good_correlation      = inps.good_correlation
    phass.min_pixels_region     = inps.min_pixels_region

    #Unwrap
    phass.unwrap(phase,correlation,unw_igram,label)

    ################################################################################
    print('#############################################')
    print('Finished unwrapping, output files:')
    print('  Unwrapped interferogram:' + os.path.join(inps.outputDir,inpName+'.unw'))
    print('  Connected components:' + os.path.join(inps.outputDir,inpName+'_int.conn'))


    #### Mising the generation of xml file, can be done with isce.object,
    # did not write that function as isce3 is installed in a different environment
    # for now use gdal2isce_xml.py
