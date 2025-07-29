from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase
from VLM.bzFRInGE import *
import h5py


def move_sbas(Exp):
    """ Move wrapped/unwrapped SBAS ifgs from 'all' to correct Exp.path_wraps/unwraps

    Idea is to just get the two nearest neighbor pairs to compare to ISCE
    """
    Exp         = ExpBase(Exp)
    path_aria   = f'{Exp.path_vlm}/Sentinel1/HR_ARIA/unwrappedPhase'
    ifgs2match0 = [op.splitext(op.basename(ifg))[0] for ifg in glob.glob(op.join(path_aria, '20*.vrt'))]
    print (f'Found {len(ifgs2match0)} ARIA ifgs')
    # reverse order
    ifgs2match  = []
    for ifg in ifgs2match0:
        sec, ref = ifg.split('_')
        ifgs2match.append(f'{ref}_{sec}')

    ## these have all 800
    path_wraps_all   = op.join(Exp.path_wraps, 'wrap_SBAS-2_all')
    path_unwraps_all = op.join(Exp.path_wraps, 'unwrap_SBAS-2_1')

    ## symlink the files that exist and match ARIA
    ifgs_poss     = np.unique([ifg[:17] for ifg in os.listdir(path_wraps_all)])
    i = 0
    for ifg in ifgs_poss:
        if ifg in ifgs2match:
            i+=1
            paths_wraps   = glob.glob(op.join(path_wraps_all, f'{ifg}*gf*'))
            paths_unwraps = glob.glob(op.join(path_unwraps_all, f'{ifg}*gf*'))

            # move the wraps
            for path in paths_wraps:
                path_new = op.join(Exp.path_wraps, op.basename(path))
                try:
                    os.symlink(path, path_new)
                except:
                    pass

            for path in paths_unwraps:
                path_new = op.join(Exp.path_unwraps, op.basename(path))
                try:
                    os.symlink(path, path_new)
                except:
                    pass
    print (f'Found {i} matching Fringe IFGS')

def dates_for_gacos(Exp):
    """ Print out the dates for copy/paste into GACOS 20 at a time"""
    path = Exp.path_ts_geo
    with h5py.File(path, 'r') as h5:
        dates = h5['date'][:]
        dates = [dt.decode('utf-8') for dt in dates]

    f = op.join(op.expanduser('~'), f'{Exp.reg}_gacos_dates.txt')
    with open(f, 'w') as fh:
        [fh.write(i + '\n') for i in dates]
    print ('Wrote dates to:', f)


if __name__ == '__main__':
    # move_sbas(HR_SBAS_90)
    dates_for_gacos(ExpBase(NYC_SR, 'SET_PM_Fast', 'NYBK', 30))

