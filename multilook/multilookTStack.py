"""
Multilook the output of topsStack by hand
"""
from contrib import *
import shutil
from isce.applications import looks

def main(path_ifgs, naz=7, nrg=19):
    # will multilook into this, maintaining same file struct as original
    path_dst = op.join(op.dirname(path_ifgs), 'interferograms')
    shutil.rmtree(path_dst)
    os.makedirs(path_dst)
    for ifg_dir in os.listdir(path_ifgs):
        # this is the ifg directory
        os.makedirs(op.join(path_dst, ifg_dir))
        for f in 'filt_fine.int filt_fine.cor'.split():
            src  = op.join(path_ifgs, ifg_dir, f)
            dst  = op.join(path_dst, ifg_dir, f)
            args = argparse.Namespace(infile=src, outfile=dst,
                                            azlooks=naz, rglooks=nrg)
            dst  = looks.main(args)
    print ('Finished multilooking ISCE')
    return


if __name__ == '__main__':
    Exp = ExpBase(HR_Base)
    path_ifgs = Exp.path_stack_ifgs
    # should have moved the full res to a new dir
    path_ifgs += '_FR'
    main(path_ifgs)
