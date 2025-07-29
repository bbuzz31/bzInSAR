from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase

from BZ import bbPlot, bbGIS
from mintpy.utils import readfile, writefile, utils as ut
from mintpy.cli import view, tsview

def plot_ts_interactive(Expi):
    # note that the GPS value is not added
    cmd  = f'{Expi.path_ts_geo} --ms 4 -u mm -v -25 25 --ylim -25 25 -c roma_r '
    cmd += f'--mask {Expi.path_mask_vup} '#--no-multilook '
    cmd += f'--periodic 1.0 0.5 '
    # cmd += f'--lalo {lalo1} '

    tsview.main(cmd.split())


def plot_ts_postage(Expi, n=None):
    # note that the GPS value is not added
    if n is not None:
        if not isinstance(n, str):
            n = f"-n {' '.join([str(nn) for nn in n])}"
    else:
        n = ''

    cmd = f'{Expi.path_ts_geo} -u mm -v -25 25 -c roma_r --noverbose -m {Expi.path_mask_vup}'
    cmd = f'{cmd} --ncols 10 {n}'

    view.main(cmd.split())


if __name__ == '__main__':
    Exp0    = ExpBase(DC_SR, 'ERA5_SET_PM_ex_Fast', 'USN7', neofs=20)
    print (f'\nPlotting: {Path(Exp0.path_ts_geo).stem}\n')

    # plot_ts_interactive(Exp0)

    n = np.arange(100, 130)
    plot_ts_postage(Exp0, n=n)
