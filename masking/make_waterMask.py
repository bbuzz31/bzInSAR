import os.path as op
from osgeo import gdal
from mintpy.utils import readfile, writefile, utils as ut


def createParser():
    parser = argparse.ArgumentParser(description=
        'Use MintPy tools to convert a gdal-readable file into a waterMask.h5.'
        '\nRequres inputs/geometry to exist',
         epilog='Examples of use:\n\t make_waterMask.py ./OSM_wmask.msk ../MintPy_LOY2',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path_mask', type=str, help='Path to gdal readable mask file')
    parser.add_argument('-mp, --path_mp', type=str, default=None, help='Path to directory containing inputs/*geo*')
    parser.add_argument('-dd, --dst_dir', type=str, default=None, help='Directory to write waterMask.h5.')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(1)

    inps = parser.parse_args(args=iargs)
    return inps


def main(path_mask, path_mp, dst_dir):
    """ path_mp should house the inputs you want to copy """
    path_wd = op.dirname(path_mask) if path_mp is None else path_mp
    dst_dir = path_wd if dst_dir is None else dst_dir

    stack_file, geom_file = ut.check_loaded_dataset(path_wd, print_msg=True)[:2]

    water_mask_file = op.join(dst_dir, 'waterMask.h5')
    _, atr     = readfile.read(geom_file, datasetName=None)
    water_mask = gdal.Open(path_mask).ReadAsArray().astype(bool)

    # ignore no-data pixels in geometry files
    ds_name_list = readfile.get_dataset_list(geom_file)
    for ds_name in ['latitude','longitude']:
        if ds_name in ds_name_list:
            print('set pixels with 0 in {} to 0 in waterMask'.format(ds_name))
            ds = readfile.read(geom_file, datasetName=ds_name)[0]
            water_mask[ds == 0] = 0

    atr['FILE_TYPE'] = 'mask'
    writefile.write(water_mask, out_file=water_mask_file, metadata=atr)

    return

    
if __name__ == '__main__':
    # path_mask = '/Users/buzzanga/data/VLM/Sentinel1/track_004/Results/HR/mask/OSM_wmask.msk'
    # path_mp = '/Users/buzzanga/data/VLM/Sentinel1/track_004/Results/HR/LOY2_Test'
    # main(path_mask, path_mp)

    inps = cmdLineParse()
    main(inps.path_mask, inps.path_mp, inps.dst_dir)
