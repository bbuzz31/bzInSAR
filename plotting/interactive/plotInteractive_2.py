"""
A tool for identifying and adjusting unwrapping errors
"""
import os, os.path as op
import time
import shutil
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from roipoly import RoiPoly

# for within stack
class MultiRoi2:
    def __init__(self, fig=None, ax=None, arr=None, parms=None, init=False, idx=0):
        """
        Parameters
        ----------
        fig: matplotlib figure
            Figure on which to draw the ROIs
        ax: matplotlib axes
           Axes on which to draw the ROIs
        """

        if fig is None:
            fig = plt.gcf()
        if ax is None:
            ax = fig.gca()

        self.fig    = fig
        self.ax     = ax
        self.arr    = arr
        self.arr_up = arr.copy()
        self.parms  = parms
        self.rois   = [{}]*self.arr.shape[0] # each dictionary corresponds to an ifg
        self.npi    = [{}]*self.arr.shape[0]
        self.txt1   = self.ax.text(-0.05, -0.1, 'Click "New ROI" to begin',
                            transform=self.ax.transAxes, fontsize=14)
        # self.ifg_idx = idx
        self.make_buttons(idx)

    ## make two more buttons that add/subtract 2pi
    def make_buttons(self, ifg_idx):
        i = ifg_idx+1
        self.ax_add_btn = plt.axes([0.7, 0.02, 0.1, 0.04], label=i)
        self.ax_finish_btn = plt.axes([0.81, 0.02, 0.1, 0.04], label=i+1)
        self.ax_down_btn  = plt.axes([0.34, 0.02, 0.1, 0.04], label=i+2)
        self.ax_up_btn    = plt.axes([0.45, 0.02, 0.1, 0.04], label=i+3)
        self.ax_clear_btn = plt.axes([0.56, 0.02, 0.1, 0.04], label=i+4)
        self.set_buttons()
        return

    def set_buttons(self, ifg_idx=0):
        btn_finish = Button(self.ax_finish_btn, 'Save')
        btn_finish.on_clicked(lambda event: self.finish(event, 'tst'))

        btn_add = Button(self.ax_add_btn, 'New ROI')
        btn_add.on_clicked(lambda event: self.add(event, ifg_idx))

        btn_up   = Button(self.ax_up_btn, 'Add 2$\pi$')
        btn_up.on_clicked(lambda event: self.increase(event, ifg_idx))

        btn_down = Button(self.ax_down_btn, 'Subtract 2$\pi$')
        btn_down.on_clicked(lambda event: self.decrease(event, ifg_idx))

        btn_clear = Button(self.ax_clear_btn, 'Clear')
        btn_clear.on_clicked(lambda event: self.clear(event, ifg_idx))
        # embed(colors='neutral')

        plt.draw()
        plt.show(block=True)

    # store the text in self.rois dict?
    def add(self, event, ifg_idx):
        """" Add a new ROI """
        # Only draw a new ROI if the previous one is completed
        if self.rois[ifg_idx]:
            if not all(r.completed for r in self.rois.values()):
                return

        roi_name = len(self.rois[ifg_idx]) # starts at 0, increments each new ROI

        logger.critical(f"Creating new ROI {roi_name}")

        ## clear the text
        self.txt1.remove()
        self.txt1  = self.ax.text(-0.05, -0.1, 'Adjust # of 2$\pi$ increments',
                                transform=self.ax.transAxes, fontsize=14)
        plt.draw()


        roi = RoiPoly(self.fig, self.ax, color='k', close_fig=False, show_fig=False)

        self.rois[ifg_idx] = {**self.rois[ifg_idx], roi_name:roi}
        self.npi[ifg_idx]  = {**self.npi[ifg_idx], roi_name:0}

        # overwrite original array with new data; for multiple aois
        self.arr[ifg_idx]            = self.arr_up[ifg_idx]


    def finish(self, event, txt):
        print (txt)
        logger.critical("Updated data??")
        # self.fig.canvas.mpl_disconnect(event)
        # plt.show(block=False)
        cid = self.fig.canvas.mpl_connect('button_press_event', self)
        event.canvas.mpl_disconnect(cid)
        return


    def increase(self, event, ifg_idx):
        """ Increase integer (of 2pi cycles) to add """
        keys  = self.rois[ifg_idx].keys()
        if len(keys) < 1:
            return
        else:
            roi_idx  = max(keys)

        self.npi[ifg_idx][roi_idx] += 1

        ## npi to add
        npi = self.npi[ifg_idx][roi_idx]
        mask_aoi    = self.rois[ifg_idx][roi_idx].get_mask(self.arr_up[ifg_idx])
        self.arr_up[ifg_idx] = np.where(mask_aoi, self.arr[ifg_idx]+(2*npi*np.pi), self.arr[ifg_idx])

        # update plot
        self.ax.imshow(self.arr_up[ifg_idx], **self.parms)
        try:
            mask_aoi    = self.rois[ifg_idx][roi_idx].get_mask(self.arr_up[ifg_idx])
            self.arr_up[ifg_idx] = np.where(mask_aoi, self.arr[ifg_idx]+(2*npi*np.pi), self.arr[ifg_idx])

            # update plot
            self.ax.imshow(self.arr_up[ifg_idx], **self.parms)
            # plt.draw()
        except:
            logger.info('No data to actually update plot')

        logger.debug(f'AOI {[ifg_idx][roi_idx]}: + {self.npi[ifg_idx][roi_idx]}*2pi')
        # print ('increase idx', self.ifg_idx)

        try: self.txt1.remove()
        except: pass

        ti       = f'Adding: {self.npi[ifg_idx][roi_idx]} $\\times 2\pi$'
        self.txt1 = self.ax.text(-0.05, -0.1, ti, transform=self.ax.transAxes, fontsize=14)
        plt.draw()


    def decrease(self, event, ifg_idx):
        """ Decrease integer (of 2pi cycles) to add """
        # get the index based on the roi
        keys  = self.rois[ifg_idx].keys()
        if len(keys) < 1:
            return
        else:
            keys    = [int(key) for key in keys]
            roi_idx = max(keys)
        self.npi[ifg_idx][roi_idx] -= 1

        ## npi to add
        npi = self.npi[ifg_idx][roi_idx]
        try:
            mask_aoi             = self.rois[ifg_idx][roi_idx].get_mask(self.arr_up[ifg_idx])
            self.arr_up[ifg_idx] = np.where(mask_aoi, self.arr[ifg_idx]+(2*npi*np.pi), self.arr[ifg_idx])

            # update plot
            self.ax.imshow(self.arr_up[ifg_idx], **self.parms)
            # plt.draw()
        except:
            logger.info('No data to actually update plot')

        try: self.txt1.remove()
        except: pass
        logger.debug(f'AOI {roi_idx}: -{self.npi[ifg_idx][roi_idx]}*2pi')
        ti       = f'Adding: {self.npi[ifg_idx][roi_idx]} $\\times 2\pi$'
        self.txt1 = self.ax.text(-0.05, -0.1, ti, transform=self.ax.transAxes, fontsize=14)
        plt.draw()


    def clear(self, event, ifg_idx):
        keys  = self.rois[ifg_idx].keys()
        # if there isnt an aoi yet
        if len(keys) < 1:
            return
        # get the newest one
        else:
            roi_idx  = max(keys)

        roi  = self.rois[ifg_idx][roi_idx]
        # no easy way to remove just the current ones (and while drawing)
        [line.remove() for line in roi.ax.get_lines()]

        npi             = self.npi[ifg_idx][roi_idx]
        ## also undo the addition/subtraction
        try:
            mask_aoi    = self.rois[ifg_idx][roi_idx].get_mask(self.arr[ifg_idx])
            self.arr_up[ifg_idx] = np.where(mask_aoi, self.arr_up[ifg_idx]-(2*npi*np.pi), self.arr[ifg_idx])
            self.ax.imshow(self.arr_up, **self.parms)
        except:
            logger.info('No data to actually update plot')
        plt.draw()


        # reset the npi?
        self.rois[ifg_idx].pop(roi_idx)
        self.npi[ifg_idx][roi_idx] = 0

assert int(mpl.__version__.replace('.', '')) < 351, 'Use a lower matplotlib version for interactivity'

def quit(event):
    print ('Early Exit.')
    print (f'Progress saved in {PATH_STACK}:unwrapPhase_temp')
    os.sys.exit()

    
def cont(event):
    plt.close()


def plot_updated(ax):
    """ Plot the updated after adjusting """
    ax.cla() # clear
    return


## ----------------------------------------------------------------- set globals
LAYER       = 'unwrapPhase' # unwrapPhase connectComponent
FS          = (12, 8) # figsize
PARMS       = dict(norm=mpl.colors.TwoSlopeNorm(0, None, None),
                     cmap='coolwarm', interpolation='nearest')
PATH_MP     = '/Users/buzzanga/data/VLM/Sentinel1/HR/MintPy_7alks_19rlks_33_15'
PATH_MASK   = op.join(PATH_MP, 'waterMask.h5')
PATH_STACK  = op.join(PATH_MP, 'inputs', 'ifgramStack.h5')
PATH_STACK0 = op.join(PATH_MP, 'inputs', 'ifgramStack_orig.h5')
# shutil.copy(PATH_STACK0, PATH_STACK)
shutil.copy(PATH_STACK, PATH_STACK0) if not op.exists(PATH_STACK0) else ''

## ----------------------------------------------------------------- get data
with h5py.File(PATH_MASK, 'r') as h5:
    mask = h5['mask'][:]
    mask = np.where(mask, 1, np.nan)


with h5py.File(PATH_STACK, 'r') as h5:
    stack  = h5[LAYER][:]
    stack *= mask
    dates  = h5['date'][:]
    orb    = h5.attrs['ORBIT_DIRECTION']

PARMS['origin'] = 'lower' if orb == 'ASCENDING' else 'upper'
IFGS = [f'{dt[0].decode("utf-8")}_{dt[1].decode("utf-8")}' for dt in dates]

## ----------------------------------------------------------------- plot data

new_stack = np.zeros_like(stack)

# have to get this out of loop for slider

for i, ifg in enumerate(IFGS):
    data      = stack[i]

    ## temporary hack for conncomp
    if LAYER == 'connectComponent':
        cc    = np.unique(data)
        ticks = range(len(cc))
        clim  = (-0.5, len(cc)-1.5)
        PARMS['cmap']   = plt.get_cmap('Dark2', len(cc)-1)
        PARMS.pop('norm', None)
    else:
        ticks, clim = None, None


    ## make the initial plot
    fig, axes = plt.subplots(figsize=FS)
    im        = axes.imshow(data, **PARMS)
    divider = make_axes_locatable(axes)
    cbar_ax = divider.append_axes('right', size='5%', pad=0.7)
    cbar_ax.set_xlabel('rad')
    cbar = plt.colorbar(im, cax=cbar_ax, ticks=ticks)
    im.set_clim(clim)

    axes.text(0.01, 1.02, ifg, transform=axes.transAxes, fontsize=14)
    axes.text(0.94, 1.02, f'{i}/{len(IFGS)}', transform=axes.transAxes, fontsize=12)

    # unbind default key bindings to prevent warning
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)

    ## make the box to take # of 2pi on top
    # x0, y0, w, h = axes.get_position().bounds

    ## begin drawing of polygon
    aois = MultiRoi(fig, axes, data, PARMS)


    if not aois.npi or all(v == 0 for v in aois.npi.values()):
        print ('No adjustments made, continuing...')
        new_stack[i] = data
        continue
    os.sys.exit()
    # else:
    #     # new plot side by side
    #     fig, axe1 = plt.subplots(figsize=FS)
    #     im        = axe1.imshow(data_new, **PARMS)
    #     axe1.set_title('Adjusted')
    #     fig.subplots_adjust(wspace=0.05)
    #     plt.show()
    #     break

    new_stack[i] = aois.arr_up

    ## save temporary stack
    with h5py.File(PATH_STACK, 'r+') as h5:
        try: del h5['unwrapPhase_temp']
        except: pass
        dset = h5.create_dataset('unwrapPhase_temp', data=new_stack)


    ## add the quit and continue buttons here
    # ax_quit_btn = plt.axes([0.40, 0.02, 0.1, 0.04])
    # btn_quit = Button(ax_quit_btn, 'Quit')
    # btn_quit.on_clicked(quit)
    #
    # ax_cont_btn = plt.axes([0.50, 0.02, 0.1, 0.04])
    # btn_cont = Button(ax_cont_btn, 'Continue')
    # btn_cont.on_clicked(cont)



## save final stack
with h5py.File(PATH_STACK, 'r+') as h5:
    try: del h5['unwrapPhase_adj']
    except: pass
    dset = h5.create_dataset('unwrapPhase_adj', data=new_stack)
    del h5['unwrapPhase_temp']

print ('Wrote unwrapPhase_adj into:', PATH_STACK)
## real time/preview
## go back and forth in stack (jump to certain one)
## undo, unclick
## colorblock in 2pi increments? grayscale
