import shutil
import h5py
from VLM.bzFRInGE import *
""" Mask using gradient of unwrapped interferograms; for unwrapping errors """


class Masker(object):
    def __init__(self, path_mp_exp, corr='SR'):
        super().__init__()
        self.path_mp_exp = path_mp_exp
        self.corr        = f'_{corr}'
        self.unw_phs     = self.load_stack()


    # 95% might work for all
    def __call__(self, thresh=90, thresh_n=85, path_mask=None, show=False):
        """
        Pixelwise, get the number of times the gradient is < its 'thresh'
        Get a number of dates it must be valid using percentage in thresh_n
        Keep (in t, y, x) the pixels where grad is valid more than tresh_n
        """
        if op.exists(path_mask):
            with h5py.File(path_mask, 'r') as h5:
                try:
                    mask0 = h5['waterMask'][:]
                except:
                    mask0 = h5['mask'][:]
            mask0 = np.where(np.isclose(mask0, 0), np.nan, mask0)
        else:
            mask0 = np.ones(self.unw_phs.shape[1:])

        # thresh = np.pi*thresh
        # phase2range = (0.055) / ( 4*np.pi)
        # unw_phs = self.unw_phs / phase2range
        unw_phs = self.unw_phs * mask0

        # date threshold; min number of dates pixel must be valid to keep
        th_n   = np.percentile(range(unw_phs.shape[0]), thresh_n)

        unw_phs = unw_phs.astype(np.float32)
        grad    = np.abs(np.gradient(unw_phs))  # lst of dt, dy, dx
        # grad1 = np.abs(np.gradient(unw_phs[0])) # arr[(dy,dx), nrows, ncols]

        th_t   = np.nanpercentile(grad[0], thresh)

        cnts_t = np.where(grad[0]<=th_t, 1, 0).sum(0) # <= important due to nans
        mask_t = np.where(cnts_t<th_n, 0, 1)

        th_y   = np.nanpercentile(grad[1], thresh)
        cnts_y = np.where(grad[1]<=th_y, 1, 0).sum(0)
        mask_y = np.where(cnts_y<th_n, 0, 1)

        th_x   = np.nanpercentile(grad[2], thresh)
        cnts_x = np.where(grad[2]<=th_x, 1, 0).sum(0)
        mask_x = np.where(cnts_x<th_n, 0, 1)

        ## apply all the masks
        mask     = mask_t * mask_y * mask_x

        src_mask = op.join(self.path_mp_exp, 'waterMask.h5')
        dst_mask = op.join(self.path_mp_exp, f'gradientMask{self.corr}.h5')
        shutil.copy(src_mask, dst_mask)

        with h5py.File(dst_mask, 'r+') as h5_new:
            try:
                data    = h5_new['mask']
                data[:] = mask
            except:
                data = h5_new['waterMask'][:]
                dset = h5_new.create_dataset('mask', data=data, dtype=bool)
                del h5_new['waterMask']

            h5_new.attrs['FILE_PATH'] = dst_mask

        print ('Wrote:', dst_mask)

        if show:
            fig, axes = plt.subplots(figsize=(10, 10), ncols=2, nrows=2,
                                                sharex=True, sharey=True)
            parms     = {'cmap': 'binary', 'interpolation': 'nearest', 'origin': 'lower'}
            masks     = [mask, mask_t, mask_y, mask_x]
            tis       = 'Combined t y x'.split()
            for ax, msk, ti in zip(axes.ravel(), masks, tis):
                ax.imshow(mask*mask0, **parms)
                ax.set_title(ti)

        return dst_mask


    def load_stack(self, use_geo=False):
        corrs      = self.corr.split('_')
        ts         = corrs[1]
        corrs      = f'_'.join(*corrs[2:]) if len(corrs)>2 else ''
        geo        = 'geo_' if use_geo else ''
        path_stack = op.join(self.path_mp_exp, 'inputs', f'{geo}ifgramStack_{ts}.h5')
        path_stack = path_stack.replace('_SBAS', '') if not op.exists(path_stack) else path_stack
        print ('Getting unwrapPhase from:', op.basename(path_stack))
        with h5py.File(path_stack, 'r') as h5:
            keep = h5['dropIfgram'][:]
            idx  = np.arange(len(keep))[keep]
            arr  = h5[f'unwrapPhase{corrs}'][idx]

        return arr



if __name__ == '__main__':
    path_mp   = op.join(op.expanduser('~'),
                    *'data VLM Sentinel1 HR MintPy_2alks_5rlks_33_15'.split())
    path_mask = op.join(path_mp, 'waterMask.h5')
    Masker(path_mp)(path_mask=path_mask, show=True)
    plt.show()
