"""
Interactively plot wrapped/unwrapped stack

This should automatically be able to plot:
    MintPy:
        ifgramStack.h5
            bperp/coherence/connectComponent/dropIfgram/unwrapPhase
                and rewrap unwrapPhase

    ARIA:
        cohStack.vrt
        connCompStack.vrt
        unwrapStack.vrt

    FRINGE
        wrapped raw

    ISCE
        wrapped/unwrapped/coherence/conncomp

"""
import h5py
import time, copy, shutil, socket
from matplotlib.widgets import Button, Slider
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
import cmocean, cmcrameri as cmc

from VLM.bzFRInGE import *
from VLM.bzFRInGE.plotting.roipoly import RoiPoly, MultiRoi
from VLM.bzFRInGE.FRINGEBase import ExpBase

FIG_PATH = f'{op.expanduser("~")}/Desktop/VLM' if \
            socket.gethostname().startswith('MT') else op.expanduser('~')

# assert int(mpl.__version__.replace('.', '')) < 351, 'Use a lower matplotlib version for interactivity'

if int(mpl.__version__.replace('.', '')) < 351:
    print('Use a lower matplotlib version for interactivity')

def prep_mask(path_mask):
    try:
        with h5py.File(path_mask, 'r') as h5:
            keys = [key for key in list(h5.keys()) if 'mask' in key.lower()]

            if len(keys) > 1:
                keys = [key for key in keys if 'water' in key.lower()]

            if keys:
                mask = h5[keys[0]][:]
                mask = np.where(mask, 1, np.nan)
            else:
                mask = None
                print (f'Couldnt find appropriate key in {path_mask}')

    except Exception as E:
        pass

    if mask is None:
        try:
            # in case it tries to open a h5 file
            mask = gdal.Open(path_mask).ReadAsArray()
            mask = np.where(mask, 1, np.nan) if mask is not None else 1

        except Exception as E:
            print ('Couldnt access mask:')
            print (E)
            mask = 1
    return mask


def imshow_cbar(im, **parms):
    divider = make_axes_locatable(im.axes)
    cbar_ax = divider.append_axes('right', size='5%', pad=0.2)
    cbar_ax.set_xlabel('rad')
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar_ax.set_xlabel('mm')
    return


## have classes that just plot one layer, and then others that combine two layers or ifgs

class unwPlotter_ARIA(ExpBase):
    """ Plot unwrapped FRInGE v ARIA """
    def __init__(self, figsize=(12, 12)):
        super().__init__(HR_SBAS_90)
        self.layer  = 'unwrapPhase'
        # use equal aspect ratio
        self.fs     = figsize
        self._set_data()


    def main(self):
        self.make_initial_plot()
        self.text = self.axes[0].text(0.01, 1.02, self.ifgs[self.idx], fontsize=14,
                                    transform=self.axes[0].transAxes)
        self.fig.canvas.draw()
        self.ax_tslider = self.fig.add_axes([0.42, 0.67, 0.475, 0.025])
        self.plot_init_time_slider(init_idx=self.idx)
        self.tslider.on_changed(self.update_slider_mouse)
        self.fig.canvas.mpl_connect('key_press_event', self.update_slider_key)
        return


    def make_initial_plot(self):
        """ make the initial plot """
        data      = self.stacka[self.idx], self.stack[self.idx]
        lbls      = 'ARIA FRInGE'.split()
        fig, axes = plt.subplots(figsize=self.fs, ncols=2, sharey=True)
        for i, (dat, ax) in enumerate(zip(data, axes)):
            parms = {'cmap': 'coolwarm', 'interpolation':'nearest',
                    'norm': mpl.colors.TwoSlopeNorm(0), 'origin': 'lower'}

            im      = ax.imshow(dat, **parms)
            ax.set_xlabel(lbls[i])

            divider      = make_axes_locatable(ax)
            cbar_ax = divider.append_axes('right', size='5%', pad=0.2)
            cbar_ax.set_xlabel('rad')
            cbar = plt.colorbar(im, cax=cbar_ax)
            # im.set_clim(clim)
            plt.tight_layout()

        # set the initial plot to be drawn on
        self.data_img1 = data[0]
        self.data_img2 = data[1]

        self.fig      = fig
        self.axes     = axes
        self.cbar_ax  = cbar_ax

        return


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

        idx = self.tslider.val

        self.update_title(self.ifgs[idx])

        # read  and update data
        im0 = self.axes[0].get_images()[0]
        im0.set_data(self.stacka[idx])
        im0.autoscale()

        im1 = self.axes[1].get_images()[0]
        im1.set_data(self.stack[idx])
        im1.autoscale()

        self.idx = idx
        self.fig.canvas.draw()

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
                # update title
                self.update_title(self.ifgs[idx])

                # update slider and image
                self.tslider.set_val(idx)

                im0 = self.axes[0].get_images()[0]
                im0.set_data(self.stacka[idx])
                im0.autoscale()

                im1 = self.axes[1].get_images()[0]
                im1.set_data(self.stack[idx])
                im1.autoscale()
                self.idx = idx

                self.fig.canvas.draw()

        return


    def update_title(self, ti):
        """ Update the interferogram title """
        self.text.remove()
        self.text = self.axes[0].text(0.01, 1.02, ti, fontsize=14,
                                    transform=self.axes[0].transAxes)

        return


    def _set_data(self):
        with h5py.File(self.path_mask_mp, 'r') as h5:
            mask = h5['mask'][:]
            mask = np.where(mask, 1, np.nan)

        with h5py.File(self.path_ifgramStack, 'r') as h5:
            stack  = h5[self.layer][:]
            stack *= mask
            dates  = h5['date'][:]
            orb    = h5.attrs['ORBIT_DIRECTION']

        path_aria  = op.join(self.path_vlm, *'Sentinel1 track_004 Results '\
                                             'HR LOY2_MintPy'.split())

        with h5py.File(op.join(path_aria, 'inputs', 'ifgramStack.h5'), 'r') as h5:
            stack_aria  = h5[self.layer][:]
            stack_aria  = np.where(np.isclose(stack_aria, 0), np.nan, stack_aria)
            stack_aria  = np.fliplr(stack_aria)
            # stack_aria *= mask
            dates_aria  = h5['date'][:]


        self.ifgs  = [f'{dt[0].decode("utf-8")}_{dt[1].decode("utf-8")}' for dt in dates]
        ifgs_aria  = [f'{dt[0].decode("utf-8")}_{dt[1].decode("utf-8")}' for dt in dates_aria]
        idx_fr = []
        idx_ar = []
        for i, ifg in enumerate(self.ifgs):
            if ifg in ifgs_aria:
                idx_fr.append(i)
                idx_ar.append(ifgs_aria.index(ifg))

        self.ifgs   = np.array(self.ifgs)[idx_fr].tolist()
        self.stack  = stack[idx_fr]
        self.stacka = stack_aria[idx_ar]

        self.nifgs = len(self.ifgs)
        self.mask  = mask
        self.idx   = 0
        return


class AdjustUnwrap(object):
    def __init__(self, path_stack, path_mask=None, layer='unwrapPhase', figsize=(12, 8), noedit=False):
        ## ----------------------------------------------------------------- set globals
        self.layer  = layer
        self.fs     = figsize
        self.parms  = dict(norm=mpl.colors.TwoSlopeNorm(0, None, None),
                             cmap='coolwarm', interpolation='nearest')
        self.path_stack = path_stack
        self.path_mp    = op.dirname(op.dirname(self.path_stack))
        self.mask       = prep_mask(path_mask)
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
        mask = self.mask


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


class PlotStack(object):
    """ This is the same as adjustUnwrap in "noedit" mode """
    def __init__(self, path_stack, layer='unwrapPhase', figsize=(12, 8)):
        ## ----------------------------------------------------------------- set globals
        self.layer  = layer
        self.fs     = figsize
        self.parms  = dict(norm=mpl.colors.TwoSlopeNorm(0, None, None),
                             cmap='coolwarm', interpolation='nearest')
        self.path_stack = path_stack
        self.path_mp    = op.dirname(op.dirname(self.path_stack))
        self.path_mask  = op.join(self.path_mp, 'waterMask.h5')
        self._set_data()


    def main(self):
        self.fig, self.axes = self.make_initial_plot(self.stack[self.idx])
        self.text = self.axes.text(0.01, 1.02, self.ifgs[self.idx], fontsize=14,
                                    transform=self.axes.transAxes)
        self.fig.canvas.draw()

        self.ax_tslider = self.fig.add_axes([0.34, 0.89, 0.470, 0.025])
        self.plot_init_time_slider(init_idx=self.idx)
        self.tslider.on_changed(self.update_slider_mouse)
        self.fig.canvas.mpl_connect('key_press_event', self.update_slider_key)


    def make_initial_plot(self, data):
        """ make the initial plot """
        fig, axes = plt.subplots(figsize=self.fs)

        im        = axes.imshow(data, **self.parms)
        divider = make_axes_locatable(axes)
        cbar_ax = divider.append_axes('right', size='5%', pad=0.3)
        cbar_ax.set_xlabel(f'rad')

        cbar = plt.colorbar(im, cax=cbar_ax, ticks=self.ticks)
        cbar.set_label(self.layer, rotation=270, labelpad=15, fontsize=15)
        im.set_clim(self.clim)

        # set the initial plot to be drawn on
        self.data_img = data

        return fig, axes


    def plot_init_time_slider(self, init_idx=0):
        """ Initialize the slider """
        val_min, val_max = 0, self.nifgs-1
        self.tslider     = Slider(self.ax_tslider, label='',
                        valinit=0, valmin=val_min, valmax=val_max, valstep=1)

        return self.tslider


    def update_slider_mouse(self, val):
        """Update Displacement Map using Slider"""
        # initialize another MultiRoi causes problems
        idx = self.tslider.val

        self.update_title(self.ifgs[idx])

        # read  and update data
        data_img = np.array(self.stack[idx])
        im = self.axes.get_images()[0]
        im.set_data(data_img)
        im.autoscale()
        self.idx = idx
        self.fig.canvas.draw()

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



        return


    def update_title(self, ti):
        """ Update the interferogram title """
        self.text.remove()
        self.text = self.axes.text(0.01, 1.02, ti, fontsize=14,
                                    transform=self.axes.transAxes)

        return


    def _set_data(self):
        ## ----------------------------------------------------------------- get data
        if op.exists(self.path_mask):
            with h5py.File(self.path_mask, 'r') as h5:
                try:
                    mask = h5['mask'][:]
                except:
                    mask = h5['waterMask'][:]
                mask = np.where(mask, 1, np.nan)
        else:
            print (f'Mask at {self.path_mask} does not exist')
            mask = 1


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
            cc    = np.unique(self.stack)
            self.ticks = range(len(cc)) # breaks scroll bar
            # self.clim  = (-0.5, len(cc)-1.5)
            # self.parms['cmap']   = plt.get_cmap('Dark2', len(cc)-1)
            # self.parms.pop('norm', None)

            cmap = plt.cm.tab20c
            # norm =  mpl.colors.BoundaryNorm(np.linspace(0, 20, 21), cmap.N, clip=True)

            norm  = mpl.colors.BoundaryNorm(np.linspace(np.nanmin(self.stack),
                        np.nanmax(self.stack), len(cc)+1), cmap.N, clip=True)

            self.ticks, self.clim = None, None


        else:
            self.ticks, self.clim = None, None

        return


# pretty good
class PlotWrapUnwrap(object):
    """ Fringe/MintPy/ARIA. Interactively Plot the wrapped/unwrapped in an ifgramStack file """
    def __init__(self, path_stack, path_mask=None, figsize=(8,8), show_dropped=True):
        # use equal aspect ratio
        self.fs     = figsize
        self.parms  = dict(norm=mpl.colors.TwoSlopeNorm(0, None, None),
                             cmap='coolwarm', interpolation='nearest')
        self.show_dropped = show_dropped
        self.mask         = prep_mask(path_mask)
        self._set_data(path_stack)


    def main(self):
        self.make_initial_plot()
        # ifg date
        self.text = self.axes[0].text(0.01, 1.02, self.ifgs[self.idx], fontsize=14,
                                    transform=self.axes[0].transAxes)
        self.fig.canvas.draw()
        # self.ax_tslider = self.fig.add_axes([0.42, 0.82, 0.475, 0.025]) # for 2 cols
        self.ax_tslider = self.fig.add_axes([0.22, 0.48, 0.575, 0.02])   # for 2 rows
        self.plot_init_time_slider(init_idx=self.idx)
        self.tslider.on_changed(self.update_slider_mouse)
        self.fig.canvas.mpl_connect('key_press_event', self.update_slider_key)
        self.fig.canvas.mpl_connect('key_press_event', self.keep_drop_key)

        self.make_buttons()
        return


    def make_initial_plot(self):
        """ make the initial plot """
        data      = self.stackWrap[self.idx], self.stackUnw[self.idx]
        fig, axes = plt.subplots(figsize=self.fs, nrows=2, sharey=True, sharex=True)
        # wrapped
        for i, (dat, ax) in enumerate(zip(data, axes)):
            if i == 0:
                clim  = [-np.pi, np.pi]
                parms = {**self.parms, 'cmap': 'cmc.cork',
                        'norm': mpl.colors.TwoSlopeNorm(0, *clim)}
                lbl   = 'Wrapped Phase'
            else:
                clim  = [np.nanquantile(dat, 0.05), np.nanquantile(dat, 0.95)]
                # in case all pos/neg
                try:
                    norm = mpl.colors.TwoSlopeNorm(0, *clim)
                except:
                    norm = mpl.colors.TwoSlopeNorm(np.mean(clim), *clim)
                parms = {**self.parms, 'cmap': 'coolwarm', 'norm': norm}
                lbl   = 'Unwrapped Phase'

            im      = ax.imshow(dat, **parms)
            im.autoscale()
            ax.set_ylabel(lbl)

            divider      = make_axes_locatable(ax)
            cbar_ax = divider.append_axes('right', size='5%', pad=0.2)
            cbar_ax.set_xlabel('rad')
            cbar = plt.colorbar(im, cax=cbar_ax)
            # im.set_clim(clim)
            # plt.tight_layout()

        # set the initial plot to be drawn on
        self.data_img1 = data[0]
        self.data_img2 = data[1]

        self.fig      = fig
        self.axes     = axes
        self.cbar_ax  = cbar_ax
        # unbind default keybindings
        self.update_title_keep()
        self.fig.subplots_adjust(hspace=0.5)
        self.fig.canvas.mpl_disconnect(self.fig.canvas.manager.key_press_handler_id)

        return


    def plot_init_time_slider(self, init_idx=0):
        """ Initialize the slider """
        val_min, val_max = 0, self.nifgs-1
        self.tslider     = Slider(self.ax_tslider, label='',
                        valinit=0, valmin=val_min, valmax=val_max, valstep=1)

        ## makes a bar plot to fill in the slider
        # bar_width = 1 / 8.

        # self.tslider.ax.bar(np.arange(self.nifgs), np.ones(self.nifgs), bar_width, align='center', facecolor='black', ecolor=None)
        self.tslider.ax.set_xticks(np.round(np.linspace(val_min, val_max, num=5)))
        self.tslider.ax.set_xlim([val_min, val_max])
        self.tslider.ax.xaxis.tick_top()
        self.tslider.ax.set_yticks([])
        self.tslider.valtext.set_visible(False)   #hide slider values
        return self.tslider


    def update_slider_mouse(self, val):
        """Update Displacement Map using Slider"""
        # initialize another MultiRoi causes problems

        idx = self.tslider.val

        self.update_title(self.ifgs[idx])

        # read  and update data
        self.axes[0].get_images()[0].set_data(self.stackWrap[idx])

        dat = self.stackUnw[idx]
        im  = self.axes[1].get_images()[0]
        im.set_data(dat)

        clim  = [np.nanquantile(dat, 0.05), np.nanquantile(dat, 0.95)]
        # in case all pos/neg
        try:
            norm = mpl.colors.TwoSlopeNorm(0, *clim)
        except:
            norm = mpl.colors.TwoSlopeNorm(np.mean(clim), *clim)

        # update the colobar also
        # clim  = np.nanquantile(self.stackUnw[idx], 0.01), np.nanquantile(self.stackUnw[idx], 0.99)
        # im.set_clim([unw.min(), unw.max()])
        im.set_norm(norm)
        im.autoscale()

        self.idx = idx
        self.update_title_keep()
        self.fig.canvas.draw()

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
                idx = min(self.idx + 1, self.stackUnw.shape[0] - 1)

            if idx is not None and idx != self.idx:
                # update title
                self.update_title(self.ifgs[idx])

                # update slider and image
                self.tslider.set_val(idx)
                self.axes[0].get_images()[0].set_data(self.stackWrap[idx])
                dat = self.stackUnw[idx]
                im = self.axes[1].get_images()[0]
                im.set_data(dat)

                clim  = [np.nanquantile(dat, 0.05), np.nanquantile(dat, 0.95)]
                # in case all pos/neg
                try:
                    norm = mpl.colors.TwoSlopeNorm(0, *clim)
                except:
                    norm = mpl.colors.TwoSlopeNorm(np.mean(clim), *clim)

                # update the colobar also
                # clim  = np.nanquantile(self.stackUnw[idx], 0.01), np.nanquantile(self.stackUnw[idx], 0.99)
                # im.set_clim([unw.min(), unw.max()])
                im.set_norm(norm)
                im.autoscale()

                self.idx = idx

                self.update_title_keep()
                self.fig.canvas.draw()

        return


    def update_title(self, ti):
        """ Update the interferogram title """
        self.text.remove()
        self.text = self.axes[0].text(0.01, 1.02, ti, fontsize=14,
                                    transform=self.axes[0].transAxes)

        return


    def update_title_keep(self):
        try:
            self.text1.remove()
        except:
            pass

        if not self.drop_idx[self.idx]:
            self.text1 = self.axes[0].text(0.62, 1.02, '(Dropped)', fontsize=18,
                                color='darkred', transform=self.axes[0].transAxes)

        else:
            self.text1 = self.axes[0].text(0.74, 1.02, '(Kept)', fontsize=18,
                            color='darkblue', transform=self.axes[0].transAxes)

        self.fig.canvas.draw()
        return


    ## buttons
    def make_buttons(self):
        """ Buttons for dropping / keeping ifg """
        ax_drop_btn = self.fig.add_axes([0.35, 0.025, 0.15, 0.04])
        ax_keep_btn = self.fig.add_axes([0.5, 0.025, 0.15, 0.04])

        ax_save_btn = self.fig.add_axes([0.145, 0.025, 0.1, 0.04])
        ax_kall_btn = self.fig.add_axes([0.85, 0.025, 0.1, 0.04])

        btn_drop = Button(ax_drop_btn, 'Drop Ifg ("1")')
        btn_drop.on_clicked(self.mark_dropped)

        btn_keep = Button(ax_keep_btn, 'Keep Ifg ("2")')
        # btn_finalize.on_clicked(lambda event: self.update_hdf5(event, True))
        btn_keep.on_clicked(self.mark_kept)

        btn_kall = Button(ax_kall_btn, 'Keep All')
        btn_kall.on_clicked(self.mark_all_kept)

        btn_save = Button(ax_save_btn, 'Save ("s")')
        btn_save.on_clicked(self.update_hdf5)

        self.fig.canvas.draw()

        plt.show(block=True)
        return


    def keep_drop_key(self, event):
        if event.inaxes and event.inaxes.figure == self.fig:
            if event.key == '1':
                self.mark_dropped()
            elif event.key == '2':
                self.mark_kept()
            elif event.key == 's':
                self.update_hdf5()
        return


    def mark_dropped(self, event=None):
        """ Write temporaray dataset to the hdf5 file """
        self.drop_idx[self.idx] = False
        self.update_title_keep()

        # with h5py.File(self.path_stack, 'r+') as h5:
        #     # delete existing temporary data
        #     del h5['dropIfgram']
        #
        #     # write the temporary data
        #     h5['dropIfgram'] = self.drop_idx
        #     print ('Wrote updated data to disk.')

        return


    def mark_kept(self, event=None):
        """ Write temporaray dataset to the hdf5 file """
        self.drop_idx[self.idx] = True
        self.update_title_keep()
        return


    def mark_all_kept(self, event=None):
        """ Write temporaray dataset to the hdf5 file """
        self.drop_idx = np.ones_like(self.drop_idx, dtype=bool)
        self.update_title_keep()
        return


    def update_hdf5(self, event=None):
        with h5py.File(self.path_ifgramStack, 'r+') as h5:
            # delete existing temporary data
            del h5['dropIfgram']

            # write the temporary data
            h5['dropIfgram'] = self.drop_idx
            print ('Wrote updated data to disk.')
        return


    def _set_data(self, path_stack):
        with h5py.File(path_stack, 'r') as h5:
            stackUnw  = h5['unwrapPhase'][:]
            if 'wrapPhase' in list(h5.keys()):
                stackWrap = h5['wrapPhase'][:]
            else:
                print (f"Couldn't get wrapPhase in {path_stack}")
                stackWrap = np.full_like(stackUnw, np.nan)

            stackUnw  *= self.mask
            stackWrap *= self.mask
            # print ('Flipping wrapped ISCE to match FRInGE')
            # stackWrap*=-1

            dates  = h5['date'][:]
            orb    = h5.attrs['ORBIT_DIRECTION']
            drop_idx = h5['dropIfgram'][:]
            proc     = h5.attrs.get('mintpy.load.processor', 'isce')



        self.parms['origin'] = 'lower' if (orb == 'ASCENDING' and proc.lower() != 'aria')  else 'upper'


        self.ifgs  = [f'{dt[0].decode("utf-8")}_{dt[1].decode("utf-8")}' for dt in dates]
        self.stackUnw  = stackUnw
        self.stackWrap = stackWrap
        self.idx   = 0

        self.drop_idx = drop_idx
        if not self.show_dropped:
            self.stackUnw  = self.stackUnw[drop_idx]
            self.stackWrap = self.stackWrap[drop_idx]
            self.drop_idx  = drop_idx[drop_idx]

        self.nifgs = self.stackUnw.shape[0]

        return



class PDFWrapUnwrap(object):
    """ Fringe/MintPy/ARIA. Write the stack of wrapped/unwrapped to a pdf """
    def __init__(self, path_stack, path_mask=None, figsize=(8,8), show_dropped=True):
        # use equal aspect ratio
        self.fs     = figsize
        self.parms  = dict(norm=mpl.colors.TwoSlopeNorm(0, None, None),
                             cmap='coolwarm', interpolation='nearest')
        self.show_dropped = show_dropped
        self.mask         = prep_mask(path_mask)
        self._set_data(path_stack)


    def __call__(self):
        """ Make all and write """
        dst    = op.join(op.join(FIG_PATH, f'{self.reg}_wrap_unwrap.pdf'))
        pdf    = PdfPages(dst)
        for i, ifg in enumerate(self.ifgs):
            data      = self.stackWrap[i], self.stackUnw[i]
            fig, axes = plt.subplots(figsize=self.fs, nrows=2, sharey=True, sharex=True)
            # wrapped
            for j, (dat, ax) in enumerate(zip(data, axes)):
                if j == 0:
                    clim  = [-np.pi, np.pi]
                    parms = {**self.parms, 'cmap': 'cmc.cork',
                            'norm': mpl.colors.TwoSlopeNorm(0, *clim)}
                    lbl   = 'Wrapped Phase'
                    ax.set_title(f'{ifg} (#{i})')
                else:
                    clim  = [np.nanquantile(dat, 0.05), np.nanquantile(dat, 0.95)]
                    # in case all pos/neg
                    try:
                        norm = mpl.colors.TwoSlopeNorm(0, *clim)
                    except:
                        norm = mpl.colors.TwoSlopeNorm(np.mean(clim), *clim)
                    parms = {**self.parms, 'cmap': 'coolwarm', 'norm': norm}
                    lbl   = 'Unwrapped Phase'

                im      = ax.imshow(dat, **parms)
                im.autoscale()
                ax.set_ylabel(lbl)

                divider      = make_axes_locatable(ax)
                cbar_ax = divider.append_axes('right', size='5%', pad=0.2)
                cbar_ax.set_xlabel('rad')
                cbar = plt.colorbar(im, cax=cbar_ax)
                fig.subplots_adjust(hspace=0.15)
            pdf.savefig(fig)
            plt.close()
        pdf.close()
        print ('Wrote:', dst)
        return

    def _set_data(self, path_stack):
        if 'Charleston' in  path_stack:
            self.reg = 'Charleston'
        elif 'NYC' in path_stack:
            self.reg = 'NYC'
        elif 'HR' in path_stack:
            self.reg = 'HR'

        with h5py.File(path_stack, 'r') as h5:
            stackUnw  = h5['unwrapPhase'][:]
            if 'wrapPhase' in list(h5.keys()):
                stackWrap = h5['wrapPhase'][:]
            else:
                print (f"Couldn't get wrapPhase in {path_stack}")
                stackWrap = np.full_like(stackUnw, np.nan)

            stackUnw  *= self.mask
            stackWrap *= self.mask
            # print ('Flipping wrapped ISCE to match FRInGE')
            # stackWrap*=-1

            dates  = h5['date'][:]
            orb    = h5.attrs['ORBIT_DIRECTION']
            drop_idx = h5['dropIfgram'][:]
            proc     = h5.attrs.get('mintpy.load.processor', 'isce')



        self.parms['origin'] = 'lower' if (orb == 'ASCENDING' and proc.lower() != 'aria')  else 'upper'


        self.ifgs  = [f'{dt[0].decode("utf-8")}_{dt[1].decode("utf-8")}' for dt in dates]
        self.stackUnw  = stackUnw
        self.stackWrap = stackWrap
        self.idx   = 0

        self.drop_idx = drop_idx
        if not self.show_dropped:
            self.stackUnw  = self.stackUnw[drop_idx]
            self.stackWrap = self.stackWrap[drop_idx]
            self.drop_idx  = drop_idx[drop_idx]

        self.nifgs = self.stackUnw.shape[0]

        return


class PlotATM(object):
    """ Plot the timeseries of: raw, atmospheric corrected, and atmospheric corrs """
    def __init__(self, path_mp_era5, path_mp_gacos, path_mask=None, figsize=(8,8)):
        # use equal aspect ratio
        self.fs     = figsize
        self.parms  = dict(norm=mpl.colors.TwoSlopeNorm(0, -20, 20),
                             interpolation='nearest')
        self.mask   = prep_mask(path_mask)
        self._set_data(path_mp_era5, path_mp_gacos)


    def main(self):
        self.make_initial_plot()
        # ifg date
        self.text = self.axes[0, 1].text(0.01, 1.2, self.dates[self.idx], fontsize=14,
                                    transform=self.axes[0, 1].transAxes)
        self.fig.canvas.draw()
        # self.ax_tslider = self.fig.add_axes([0.42, 0.82, 0.475, 0.025]) # for 2 cols
        self.ax_tslider = self.fig.add_axes([0.22, 0.1, 0.575, 0.02])   # for 2 rows
        self.plot_init_time_slider(init_idx=self.idx)
        self.tslider.on_changed(self.update_slider_mouse)
        self.fig.canvas.mpl_connect('key_press_event', self.update_slider_key)
        # self.fig.canvas.mpl_connect('key_press_event', self.keep_drop_key)

        return


    def make_initial_plot(self):
        """ make the initial plot """
        data = self.ERA5, self.GACOS
        cmap0 = 'cmc.roma_r'
        cmap1 = 'cmc.bam'
        tis   = 'Raw Corrected Topo'.split()

        fig, axes = plt.subplots(figsize=self.fs, nrows=2, ncols=3, sharey=True, sharex=True)

        for i in range(3):
            if i < 2:
                parms = {**self.parms, 'cmap': cmap0}
            else:
                parms = {**self.parms, 'cmap': cmap1, 'vmin': -3500, 'vmax': -2500}
                parms.pop('norm')

            im   = axes[0, i].imshow(self.ERA5[0][i], **parms)
            axes[0, i].set_title(tis[i])
            imshow_cbar(im)

            im = axes[1, i].imshow(self.GACOS[0][i], **parms)
            imshow_cbar(im)

        axes[0, 0].set_ylabel('ERA5')
        axes[1, 0].set_ylabel('GACOS')

        self.fig      = fig
        self.axes     = axes
        self.fig.subplots_adjust(hspace=-0.65, wspace=0.5)
        self.fig.tight_layout()
        self.fig.canvas.mpl_disconnect(self.fig.canvas.manager.key_press_handler_id)

        return


    def plot_init_time_slider(self, init_idx=0):
        """ Initialize the slider """
        val_min, val_max = 0, self.ndt-1
        self.tslider     = Slider(self.ax_tslider, label='',
                        valinit=0, valmin=val_min, valmax=val_max, valstep=1)

        ## makes a bar plot to fill in the slider
        # bar_width = 1 / 8.

        # self.tslider.ax.bar(np.arange(self.nifgs), np.ones(self.nifgs), bar_width, align='center', facecolor='black', ecolor=None)
        self.tslider.ax.set_xticks(np.round(np.linspace(val_min, val_max, num=5)))
        self.tslider.ax.set_xlim([val_min, val_max])
        self.tslider.ax.xaxis.tick_top()
        self.tslider.ax.set_yticks([])
        self.tslider.valtext.set_visible(False)   #hide slider values
        return self.tslider


    def update_slider_mouse(self, val):
        """Update Displacement Map using Slider"""
        # initialize another MultiRoi causes problems

        idx = self.tslider.val

        self.update_title(self.dates[idx])

        # read and update data
        for col in range(3):
            im0 = self.axes[0, col].get_images()[0].set_data(self.ERA5[col][idx])
            im1 = self.axes[1, col].get_images()[0].set_data(self.GACOS[col][idx])

        # only autoscale last one
        # fucking thing breaks it
        # for i, im in enumerate([im0, im1]):
            # lims = np.nanmin(im.get_array()), np.nanmax(im.get_array())
            # im.set_norm(mpl.colors.BoundaryNorm(lims, 256))
            # im.autoscale()


        self.idx = idx
        self.fig.canvas.draw()

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
                idx = min(self.idx + 1, self.ERA5.shape[1] - 1)

            if idx is not None and idx != self.idx:
                # update title
                self.tslider.set_val(idx)
                self.update_title(self.dates[idx])

                for col in range(3):
                    im0 = self.axes[0, col].get_images()[0].set_data(self.ERA5[col][idx])
                    im1 = self.axes[1, col].get_images()[0].set_data(self.GACOS[col][idx])

                self.idx = idx
                self.fig.canvas.draw()

        return


    def update_title(self, ti):
        """ Update the interferogram title """
        self.text.remove()
        self.text = self.axes[0, 1].text(0.01, 1.2, ti, fontsize=14,
                                    transform=self.axes[0, 1].transAxes)

        return


    def _set_data(self, path_mp_era5, path_mp_gacos):
        dct_ts = {}
        for ti, path_mp in zip(['ERA5', 'GACOS'], [path_mp_era5, path_mp_gacos]):
            lst_ts = []
            for f in f'timeseries timeseries_{ti} inputs/{ti}'.split():
                path = op.join(path_mp, f'{f}.h5')
                with h5py.File(path, 'r') as h5:
                    ts  = h5['timeseries'][:] * self.mask

                    orb = h5.attrs['ORBIT_DIRECTION']
                    proc     = h5.attrs.get('mintpy.load.processor', 'isce')
                    dates    = h5['date'][:]
                lst_ts.append(ts.astype(np.float32))
            if len(lst_ts) < 3:
                if ti == 'ERA5':
                    raise Exception ('Couldnt find all the ERA5 data')
                # if missing GACOS put in zeros
                else:
                    arr = np.zeros_like(dct_ts['ERA5'])
            ## everything goes according to plan
            else:
                arr = np.stack(lst_ts)

            dct_ts[ti] = arr*1000 # 3 x time x lat x lon


        self.parms['origin'] = 'lower' if (orb == 'ASCENDING' and proc.lower() != 'aria')  else 'upper'
        self.dates = [f'{dt.decode("utf-8")}' for dt in dates]
        self.idx   = 0

        self.ndt    = len(self.dates)
        self.ERA5   = dct_ts['ERA5']
        self.GACOS  = dct_ts['GACOS']
        return


## not done
class PlotSET(object):
    """ Plot the timeseries of: raw, atmospheric corrected, and atmospheric corrs """
    def __init__(self, path_mp_set, path_mask=None, figsize=(12,8)):
        # use equal aspect ratio
        self.fs     = figsize
        self.parms  = dict(norm=mpl.colors.TwoSlopeNorm(0, None, None),
                             interpolation='nearest')
        self.mask         = prep_mask(path_mask)
        self._set_data(path_mp_era5, path_mp_gacos)


    def main(self):
        self.make_initial_plot()
        # ifg date
        self.text = self.axes[0].text(0.01, 1.02, self.ifgs[self.idx], fontsize=14,
                                    transform=self.axes[0].transAxes)
        self.fig.canvas.draw()
        # self.ax_tslider = self.fig.add_axes([0.42, 0.82, 0.475, 0.025]) # for 2 cols
        self.ax_tslider = self.fig.add_axes([0.22, 0.48, 0.575, 0.02])   # for 2 rows
        self.plot_init_time_slider(init_idx=self.idx)
        # self.tslider.on_changed(self.update_slider_mouse)
        # self.fig.canvas.mpl_connect('key_press_event', self.update_slider_key)
        # self.fig.canvas.mpl_connect('key_press_event', self.keep_drop_key)

        return


    def make_initial_plot(self):
        """ make the initial plot """
        cmap  = 'cmc.roma_r'
        cmap1 = 'cmc.bam'

        fig, axes = plt.subplots(figsize=self.fs, ncols=3, sharey=True, sharex=True)

        im = axes[0, 0].imshow(self.ts[0][0], cmap=cmap, **self.parms)
        axes[0,0].set_title('RAW')
        imshow_cbar(im)

        im = axes[0, 1].imshow(self.ts[1][0], cmap=cmap, **self.parms)
        axes[0, 1].set_title('Corrected')
        imshow_cbar(im)

        im = axes[0, 2].imshow(self.ts[2][0], cmap=cmap1, **self.parms)
        axes[0, 1].set_title('SET')
        imshow_cbar(im)


        # set the initial plot to be drawn on
        self.data_img1 = data[0][0]
        self.data_img2 = data[0][1]
        self.data_img3 = data[0][2]

        self.fig      = fig
        self.axes     = axes
        # self.cbar_ax  = cbar_ax
        # unbind default keybindings
        # self.update_title_keep()
        self.fig.subplots_adjust(hspace=0.25, wspace=0.5)
        self.fig.canvas.mpl_disconnect(self.fig.canvas.manager.key_press_handler_id)

        return


    def plot_init_time_slider(self, init_idx=0):
        """ Initialize the slider """
        val_min, val_max = 0, self.nifgs-1
        self.tslider     = Slider(self.ax_tslider, label='',
                        valinit=0, valmin=val_min, valmax=val_max, valstep=1)

        ## makes a bar plot to fill in the slider
        # bar_width = 1 / 8.

        # self.tslider.ax.bar(np.arange(self.nifgs), np.ones(self.nifgs), bar_width, align='center', facecolor='black', ecolor=None)
        self.tslider.ax.set_xticks(np.round(np.linspace(val_min, val_max, num=5)))
        self.tslider.ax.set_xlim([val_min, val_max])
        self.tslider.ax.xaxis.tick_top()
        self.tslider.ax.set_yticks([])
        self.tslider.valtext.set_visible(False)   #hide slider values
        return self.tslider


    def update_slider_mouse(self, val):
        """Update Displacement Map using Slider"""
        # initialize another MultiRoi causes problems

        idx = self.tslider.val

        self.update_title(self.dates[idx])

        # read and update data
        for col in range(3):
            im = self.axes[0, col].get_images()[0].set_data(self.ERA5[col][idx])
            im.autoscale()
            im = self.axes[1, col].get_images()[0].set_data(self.GACOS[col][idx])
            im.autoscale()

        self.idx = idx
        self.update_title_keep()
        self.fig.canvas.draw()

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
                idx = min(self.idx + 1, self.stackUnw.shape[0] - 1)

            if idx is not None and idx != self.idx:
                # update title
                self.update_title(self.ifgs[idx])

                # update slider and image
                self.tslider.set_val(idx)
                self.axes[0].get_images()[0].set_data(self.stackWrap[idx])
                dat = self.stackUnw[idx]
                im = self.axes[1].get_images()[0]
                im.set_data(dat)

                clim  = [np.nanquantile(dat, 0.05), np.nanquantile(dat, 0.95)]
                # in case all pos/neg
                try:
                    norm = mpl.colors.TwoSlopeNorm(0, *clim)
                except:
                    norm = mpl.colors.TwoSlopeNorm(np.mean(clim), *clim)

                # update the colobar also
                # clim  = np.nanquantile(self.stackUnw[idx], 0.01), np.nanquantile(self.stackUnw[idx], 0.99)
                # im.set_clim([unw.min(), unw.max()])
                im.set_norm(norm)
                im.autoscale()

                self.idx = idx

                self.update_title_keep()
                self.fig.canvas.draw()

        return


    def update_title(self, ti):
        """ Update the interferogram title """
        self.text.remove()
        self.text = self.axes[0].text(0.01, 1.02, ti, fontsize=14,
                                    transform=self.axes[0].transAxes)

        return


    def _set_data(self, path_mp_era5, path_mp_gacos):
        lst_ts = []
        for f in f'timeseries timeseries_{ti} inputs/{ti}'.split():
            path = op.join(path_mp, f'{f}.h5')
            with h5py.File(path, 'r') as h5:
                ts  = h5['timeseries'][:] * self.mask

                orb = h5.attrs['ORBIT_DIRECTION']
                proc     = h5.attrs.get('mintpy.load.processor', 'isce')
                dates    = h5['date'][:]
            lst_ts.append(ts)

        arr = np.stack(lst_ts) * 1000 # 3 x time x lat x lon

        self.parms['origin'] = 'lower' if (orb == 'ASCENDING' and proc.lower() != 'aria')  else 'upper'
        self.dates = [f'{dt.decode("utf-8")}' for dt in dates]
        self.idx   = 50

        self.ndt    = len(self.dates)
        self.ts     = arr
        return


if __name__ == '__main__':
    # path_mp    = '/Users/buzzanga/data/VLM/Sentinel1/HR/MintPy_7alks_19rlks_33_15'
    # path_stack = op.join(path_mp, 'inputs', 'ifgramStack_SR.h5')

    # PlotObj = unwPlotter_ARIA()
    # PlotObj.main()

    # path_mp    = '/Users/buzzanga/data/VLM/Sentinel1/HR/MintPy_7alks_19rlks_33_15'
    # path_stack = op.join(path_mp, 'inputs', 'ifgramStack_SR.h5')

    path_stack = f'/Volumes/BB_4TB/data/VLM/Sentinel1/Charleston/ifgramStack_SR_OG.h5'
    path_mask  = f'/Volumes/BB_4TB/data/VLM/Sentinel1/Charleston/waterMask.h5'
    # AdjustUnwrap(path_stack, path_mask, noedit=False).main()

    # PlotStack(path_stack, layer='connectComponent').main()
    Plotter = PlotWrapUnwrap(path_stack, path_mask, figsize=(10, 10))
    Plotter.main()

    # PDFWrapUnwrap(path_stack, path_mask, figsize=(10, 10))()


    # path_base  = '/Users/buzzanga/data/VLM/Sentinel1/SC/MintPy_2alks_5rlks_33_15'
    # PlotATM(op.join(path_base, 'SCHA_ERA5'),
            # op.join(path_base, 'SCHA_GACOS'),
            # op.join(path_base, 'waterMask.h5'))



    plt.show()
