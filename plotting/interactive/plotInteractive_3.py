"""
A tool for identifying and adjusting unwrapping errors
"""
import os, os.path as op
import copy
import time
import shutil
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
from mpl_toolkits.axes_grid1 import make_axes_locatable
from contrib.plotting.roipoly import RoiPoly, MultiRoi


assert int(mpl.__version__.replace('.', '')) < 351, 'Use a lower matplotlib version for interactivity'


class Plotter(object):
    def __init__(self, path_stack, layer='unwrapPhase', figsize=(12, 8), noedit=False):
        ## ----------------------------------------------------------------- set globals
        self.layer  = layer
        self.fs     = figsize
        self.parms  = dict(norm=mpl.colors.TwoSlopeNorm(0, None, None),
                             cmap='coolwarm', interpolation='nearest')
        self.path_stack = path_stack
        self.path_mp    = op.dirname(op.dirname(self.path_stack))
        self.path_mask  = op.join(self.path_mp, 'waterMask.h5')
        self.noedit     = noedit # just scroll through the stack
        self._set_data()


    def main(self):
        self.fig, self.axes = self.make_initial_plot(self.stack[self.idx])
        self.fig.suptitle('Scroll left/right to begin')
        self.text = self.axes.text(0.01, 1.02, self.ifgs[self.idx], fontsize=14,
                                    transform=self.axes.transAxes)
        self.fig.canvas.draw()

        self.ax_tslider = self.fig.add_axes([0.325, 0.89, 0.475, 0.025])
        self.plot_init_time_slider(init_idx=self.idx)
        self.tslider.on_changed(self.update_slider_mouse)
        self.fig.canvas.mpl_connect('key_press_event', self.update_slider_key)
        self.fig.canvas.mpl_connect('key_press_event', self.save_quit_key)

        # self.Obj = MultiRoi2(self.fig, self.axes, self.stack, idx=self.idx)
        # self.Obj.set_buttons()

        ## initialize the object with the data
        # self.fig2, self.axes2 = self.make_initial_plot(self.stack[self.idx])
        self.make_save_buttons()


    def make_initial_plot(self, data):
        """ make the initial plot """
        fig, axes = plt.subplots(figsize=self.fs)
        im        = axes.imshow(data, **self.parms)
        divider = make_axes_locatable(axes)
        cbar_ax = divider.append_axes('right', size='5%', pad=0.7)
        cbar_ax.set_xlabel('rad')
        cbar = plt.colorbar(im, cax=cbar_ax, ticks=self.ticks)
        im.set_clim(self.clim)

        # set the initial plot to be drawn on
        self.data_img = data

        return fig, axes


    def plot_init_time_slider(self, init_idx=0):
        """ Initialize the slider """
        val_min, val_max = 0, self.nifgs-1
        self.tslider     = Slider(self.ax_tslider, label='',
                        valinit=0, valmin=val_min, valmax=val_max, valstep=1)

        ## makes a bar plot to fill in the slider
        bar_width = 1 / 8.

        self.tslider.ax.bar(np.arange(self.nifgs), np.ones(self.nifgs), bar_width, align='center', facecolor='black', ecolor=None)
        self.tslider.ax.set_xticks(np.round(np.linspace(val_min, val_max, num=5)))
        self.tslider.ax.set_xlim([val_min, val_max])
        self.tslider.ax.xaxis.tick_top()
        self.tslider.ax.set_yticks([])
        self.tslider.valtext.set_visible(False)   #hide slider values
        return self.tslider


    def update_slider_mouse(self, val):
        """Update Displacement Map using Slider"""
        # initialize another MultiRoi causes problems

        # idx = np.argmin(np.abs(np.array(self.yearList) - self.tslider.val))
        idx = self.tslider.val

        # replace original data
        # data_new = self.axes2.get_images()[0].get_array()
        # self.stack[self.idx] = data_new
        # plt.draw()

        self.update_title(self.ifgs[idx])

        # read  and update data
        data_img = np.array(self.stack[idx])
        im = self.axes.get_images()[0]
        im.set_data(data_img)
        im.autoscale()
        self.idx = idx
        self.fig.canvas.draw()

        ## new plot for drawing
        # plt.close(self.fig2)
        # self.fig2, self.axes2 = self.make_initial_plot(self.stack[idx])

        ## begin drawing of polygon
        # aois = MultiRoi2(self.fig2, self.axes2, self.data_img, self.parms, self.ifgs[idx])
        # MultiRoi2(self.fig, self.axes, self.stack, idx=self.idx)

        # plt.close('all')
        return


    def update_slider_key(self, event):
        """Slide images with left/right key on keyboard"""
        ## from previous event; use this to write older array
        # print (self.idx) use this to overwrite
        self.fig.suptitle('')
        if event.inaxes and event.inaxes.figure == self.fig:
            idx = None
            if event.key == 'left':
                idx = max(self.idx - 1, 0)
            elif event.key == 'right':
                idx = min(self.idx + 1, self.stack.shape[0] - 1)

            if idx is not None and idx != self.idx:
                # replace original data (except on first)
                if hasattr(self, 'fig2'):
                    data_new = self.axes2.get_images()[0].get_array()
                    self.stack[self.idx] = data_new
                    plt.close(self.fig2)
                    plt.draw()

                # update title
                self.update_title(self.ifgs[idx])

                # read data
                self.data_img = np.array(self.stack[idx])

                # update slider and image
                self.tslider.set_val(idx)
                im = self.axes.get_images()[0]
                im.set_data(self.data_img)
                im.autoscale()
                self.idx = idx

                self.fig.canvas.draw()

                if not self.noedit:
                    ## new plot for drawing
                    self.fig2, self.axes2 = self.make_initial_plot(self.stack[idx])
                    ## begin drawing of polygon
                    Obj = MultiRoi(self.fig2, self.axes2, self.data_img,
                                                                self.ifgs[idx])


        return


    def update_hdf5(self, event=None, finalize=False):
        """ Write temporaray dataset to the hdf5 file """
        temp_layer = f'{self.layer}_temp'
        with h5py.File(self.path_stack, 'r+') as h5:
            # ensure original dataset exists
            if not f'{self.layer}_og' in h5.keys():
                h5[f'{self.layer}_og'] = h5[self.layer][:]

            # delete existing temporary data
            if temp_layer in h5.keys():
                del h5[temp_layer]

            # write the temporary data
            h5[temp_layer] = self.stack
            print ('Wrote updated data to disk.')

            # write the new stack into the regular arr
            if finalize:
                del h5[temp_layer]
                del h5[self.layer]
                h5[self.layer] = self.stack
                plt.close('all')
                print ('Wrote final data to disk.')

        return


    def save_quit_key(self, event):
        """Slide images with left/right key on keyboard"""
        if event.inaxes and event.inaxes.figure == self.fig:
            if event.key == 'u':
                self.update_hdf5()
                print ('Saved temp data')
            elif event.key == 'f':
                self.update_hdf5(finalize=True)
                print ('Wrote final data.')
        return


    def update_title(self, ti):
        """ Update the interferoram title """
        self.text.remove()
        self.text = self.axes.text(0.01, 1.02, ti, fontsize=14,
                                    transform=self.axes.transAxes)

        return


    def make_save_buttons(self):
        """ Buttons dont actually work, just guidelines """
        ax_finalize_btn = self.fig.add_axes([0.31, 0.02, 0.1, 0.04])
        ax_update_btn   = self.fig.add_axes([0.21, 0.02, 0.1, 0.04])

        btn_update = Button(ax_update_btn, 'Update ("u")')
        btn_update.on_clicked(lambda event: self.update_hdf5(event, False))

        btn_finalize = Button(ax_finalize_btn, 'Finalize ("f")')
        btn_finalize.on_clicked(lambda event: self.update_hdf5(event, True))
        self.fig.canvas.draw()
        plt.show(block=True)
        return


    def _set_data(self):
        ## ----------------------------------------------------------------- get data
        with h5py.File(self.path_mask, 'r') as h5:
            try:
                mask = h5['mask'][:]
            except:
                mask = h5['waterMask'][:]
            mask = np.where(mask, 1, np.nan)


        with h5py.File(self.path_stack, 'r') as h5:
            stack  = h5[self.layer][:]
            stack *= mask
            dates  = h5['date'][:]
            orb    = h5.attrs['ORBIT_DIRECTION']

        self.parms['origin'] = 'lower' if orb == 'ASCENDING' else 'upper'
        self.ifgs  = [f'{dt[0].decode("utf-8")}_{dt[1].decode("utf-8")}' for dt in dates]
        self.stack = stack
        self.stack_OG = self.stack.copy()
        self.nifgs = self.stack.shape[0]
        self.mask  = mask
        self.idx   = 0
        # self.yearList = bbTS.date2dec(dates)

        if self.layer == 'connectComponent':
            cc         = np.unique(self.stack)
            self.ticks = range(len(cc))
            self.clim  = (-0.5, len(cc)-1.5)
            self.parms['cmap']   = plt.get_cmap('Dark2', len(cc)-1)
            self.parms.pop('norm', None)
        else:
            self.ticks, self.clim = None, None

        return


if __name__ == '__main__':
    path_mp    = '/Users/buzzanga/data/VLM/Sentinel1/HR/MintPy_7alks_19rlks_33_15'
    path_stack = op.join(path_mp, 'inputs', 'ifgramStack_SR.h5')

    # path_mp = '/Users/buzzanga/data/VLM/Sentinel1/track_004/Results/HR/LOY2_MintPy'
    # path_stack = op.join(path_mp, 'inputs', 'ifgramStack.h5')
    Plotter(path_stack, noedit=False).main()
    # plt.show()
