""" Write a file that has the index of dates you want to exclude  """
# from experiments import *
from contrib import *
from contrib.utils.mp_readers import stack2df, get_ifg_dates


class EditDates(object):
    def __init__(self, path_stack, nconn=2, max_tbase=36):
        self.path_stack = path_stack
        self.df = stack2df(self.path_stack, drop=False)[0]
        self.ifgs = get_ifg_dates(self.path_stack)
        self.ann_pairs = get_annual_pairs(self.ifgs)
        self.sbas_pairs = get_sbas_pairs(self.ifgs, nconn, max_tbase)

    def __call__(self):
        """ Write the files with the dates to exclude in them """
        path_wd = op.dirname(op.dirname(self.path_stack))
        ## write the annual only network
        dst_ann = op.join(path_wd, 'ANN.txt')
        with open(dst_ann, 'w') as fh:
            for i, ifg in enumerate(self.ifgs):
                if ifg in self.ann_pairs:
                    continue
                print (f'{i} #{ifg}', file=fh)

        dst_sbas = op.join(path_wd, 'SBAS.txt')
        with open(dst_sbas, 'w') as fh:
            for i, ifg in enumerate(self.ifgs):
                if ifg in self.sbas_pairs:
                    continue
                print (f'{i} #{ifg}', file=fh)

        dst_both = op.join(path_wd, 'ANN+SBAS.txt')
        with open(dst_both, 'w') as fh:
            for i, ifg in enumerate(self.ifgs):
                if ifg in self.sbas_pairs or ifg in self.ann_pairs:
                    continue
                print (f'{i} #{ifg}', file=fh)

        [print ('Wrote:', dst) for dst in [dst_ann, dst_sbas, dst_both]]
        return


## these next two should be combined
class SymlinkCustom(ExpBase):
    """ Symlink the ifgs in date_list from path_src to path_dst

    For running a custom prep fringe with SR, SBAS and ANN
    """

    def __init__(self, exp, path_srcs, path_dst, date_list, ext='.unw'):
        super().__init__(exp, v=False)
        self.path_srcs = path_srcs
        self.path_dst  = path_dst
        self.ext       = ext
        self.date_list = self._read_date_list(date_list)
        os.makedirs(self.path_dst, exist_ok=True)


    def __call__(self):
        path_srcs = glob.glob(op.join(self.path_srcs, f'*{self.ml}*{self.gap_fill}*{self.ext}'))

        for ifg in self.date_list:
            for f in path_srcs:
                if ifg in f:
                    for ext in ['', '.vrt', '.xml']:
                        src = f+ext
                        dst = op.join(self.path_dst, op.basename(src))
                        # skip if this extension not there or its linking to itself
                        if not op.exists(src) or src == dst:
                            continue

                        os.unlink(dst) if op.islink(dst) else ''
                        os.symlink(src, dst)


    def _read_date_list(self, date_list):
        with open(date_list, 'r') as fh:
            dates = [line.strip() for line in fh]
        return dates


class WriteRefDate(ExpBase):
    """ Write a file that contains the ifgs you want to give to prep_fringe """
    def __init__(self, exp):
        super().__init__(exp)
        # self.make_all()
        # self.make_aria()


    def make_all(self):
        """" Get all ifgs everywhere """
        unwrap_dirs = get_unw_dirs(self.path_ps_ds)

        ifgs = []
        ## iterate over all ifgs, adding if it doesnt already exist
        for unw_dir in unwrap_dirs:
            pot_ifgs = [f[:17] for f in os.listdir(unw_dir)]
            for ifg in pot_ifgs:
                if not ifg in ifgs:
                    ifgs.append(ifg)

        ifgs     = sorted(ifgs)
        path_ref = op.join(self.path_ps_ds, 'date12_all.txt')
        with open(path_ref, 'w') as fh:
            [print (ifg, file=fh) for ifg in ifgs]
        print (f'Wrote all {len(ifgs)} ifgs to:', path_ref)
        return


    def make_aria(self):
        """ make the date12 using an ifgramStack from ARIA """
        from contrib.utils import mp_readers
        reg = self.reg if self.reg == 'HR' else 'Charleston'
        mp  = f'MintPy_{DCT_REG[reg][3][0]}'

        path_aria  = op.join(op.dirname(self.path_wd), f'{reg}_ARIA')
        path_stack = op.join(path_aria, mp, 'inputs', 'ifgramStack.h5')

        if op.exists(path_stack):
            print ('Found stack in ARIA directory:', path_aria)
        else:
            print ('Could not find stack in ARIA directory:', path_aria)
            return

        df_stack = mp_readers.stack2df(path_stack, drop=True)[0]
        ifgs     = df_stack.apply(lambda x: f"{x['ref']}_{x['sec']}", 1)

        path_ref = op.join(self.path_ps_ds, 'date12_ARIA.txt')
        ifgs.to_csv(path_ref, index=None, header=None)
        # with open(path_ref, 'w') as fh:
        #     [print (ifg, file=fh) for ifg in ifgs]
        print (f'Wrote all {len(ifgs)} ifgs to:', path_ref)

        return path_ref


if __name__ == '__main__':
    # path_ifgramStack = f'{os.getenv("dataroot")}/VLM/Sentinel1/'\
    #                     f'HR_ARIA/MintPy_LOY2/inputs/ifgramStack.h5'
    # EditDates(path_ifgramStack)()

    root         = '/u/leffe-data2/buzzanga/data/VLM/Sentinel1/HR2/PS_DS_33_15'
    path_unwraps = op.join(root, 'SBAS', 'unwrap_SBAS-2')
    path_dst     = op.join(root, 'PS_DS_33_15', 'Custom')
    path_ref     = op.join(root, 'date12_163.txt')
    # SymlinkCustom(HR_90c, path_unwraps, path_dst, path_ref)()

    WriteRefDate(SC_Base)
