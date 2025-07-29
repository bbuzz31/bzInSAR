"""Draw polygon regions of interest (ROIs) in matplotlib images,

Similar to Matlab's roipoly function.
credit:https://github.com/jdoepfert/roipoly.py/blob/master/roipoly/roipoly.py
"""

import sys
import logging
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.widgets import Button

logger = logging.getLogger(__name__)
logging.basicConfig()
logger.setLevel(30)

warnings.simplefilter('always', DeprecationWarning)

from BZ import *

## ----------------------------------------------------------------- helper fns
def deprecation(message):
    warnings.warn(message, DeprecationWarning)


def submit(text):
    """ Call text_box.text to get the value """
    return


class RoiPoly:

    def __init__(self, fig=None, ax=None, color='b',
                show_fig=True, close_fig=True):
        """
        Parameters
        ----------
        fig: matplotlib figure
            Figure on which to create the ROI
        ax: matplotlib axes
            Axes on which to draw the ROI
        color: str
           Color of the ROI
        roicolor: str
            deprecated, use `color` instead
        show_fig: bool
            Display the figure upon initializing a RoiPoly object
        close_fig: bool
            Close the figure after finishing ROI drawing
        """
        if fig is None:
            fig = plt.gcf()
        if ax is None:
            ax = plt.gca()

        self.start_point = []
        self.end_point = []
        self.previous_point = []
        self.x = []
        self.y = []
        self.line  = None
        self.completed = False  # Has ROI drawing completed?
        self.color = color
        self.fig = fig
        self.ax = ax
        self.close_figure = close_fig

        # Mouse event callbacks
        self.__cid1 = self.fig.canvas.mpl_connect(
            'motion_notify_event', self.__motion_notify_callback)
        self.__cid2 = self.fig.canvas.mpl_connect(
            'button_press_event', self.__button_press_callback)

        if show_fig:
            self.show_figure()


    @staticmethod
    def show_figure():
        if sys.flags.interactive:
            plt.show(block=False)
        else:
            plt.show(block=True)


    def get_mask(self, image):
        """Get binary mask of the ROI polygon.
        Parameters
        ----------
        image: numpy array (2D)
            Image that the mask should be based on. Only used for determining
            the shape of the binary mask (which is made equal to the shape of
            the image)
        Returns
        -------
        numpy array (2D)
        """
        ny, nx = np.shape(image)
        poly_verts = ([(self.x[0], self.y[0])]
                      + list(zip(reversed(self.x), reversed(self.y))))
        # Create vertex coordinates for each grid cell...
        # (<0,0> is at the top left of the grid in this system)
        x, y = np.meshgrid(np.arange(nx), np.arange(ny))
        x, y = x.flatten(), y.flatten()
        points = np.vstack((x, y)).T

        roi_path = MplPath(poly_verts)
        mask = roi_path.contains_points(points).reshape((ny, nx))
        return mask


    def display_roi(self, ax=None, **linekwargs):
        line = plt.Line2D(self.x + [self.x[0]], self.y + [self.y[0]],
                          color=self.color, **linekwargs)
        ax = plt.gca() if ax is None else ax
        ax.add_line(line)
        plt.draw()


    def get_mean_and_std(self, image):
        """Get statistics about pixel values of an image inside the ROI.
        Parameters
        ----------
        image: numpy array (2D)
            Image on which the statistics should be calculated
        Returns
        -------
        list of float:
            mean and standard deviation of the pixel values inside the ROI
        """
        mask = self.get_mask(image)
        mean = np.mean(np.extract(mask, image))
        std = np.std(np.extract(mask, image))
        return mean, std


    def display_mean(self, image, **textkwargs):
        """Display statistics about pixel values of an image inside the ROI.
        Parameters
        ----------
        image: numpy array (2D)
            Image on which the statistics should be calculated
        Returns
        -------
        None
        """
        mean, std = self.get_mean_and_std(image)
        string = "%.3f +- %.3f" % (mean, std)
        plt.text(self.x[0], self.y[0],
                 string, color=self.color,
                 bbox=dict(facecolor='w', alpha=0.6), **textkwargs)


    def get_roi_coordinates(self):
        """Get co-ordinates of the ROI polygon.
        Returns
        -------
        numpy array (2D)
        """
        roi_coordinates = list(zip(self.x, self.y))
        return roi_coordinates


    def __motion_notify_callback(self, event):
        if event.inaxes == self.ax:
            x, y = event.xdata, event.ydata
            if ((event.button is None or event.button == 1) and
                    self.line is not None):
                # Move line around
                x_data = [self.previous_point[0], x]
                y_data = [self.previous_point[1], y]
                logger.debug(f"draw line x: {x_data} y: {y_data}")
                self.line.set_data(x_data, y_data)
                self.fig.canvas.draw()


    def __button_press_callback(self, event):
        if event.inaxes == self.ax:
            x, y = event.xdata, event.ydata
            ax = event.inaxes
            if event.button == 1 and not event.dblclick:
                logger.debug("Received single left mouse button click")
                if self.line is None:  # If there is no line, create a line
                    self.line = plt.Line2D([x, x], [y, y], linewidth=1,
                                           marker='o', color=self.color)
                    self.start_point = [x, y]
                    self.previous_point = self.start_point
                    self.x = [x]
                    self.y = [y]

                    ax.add_line(self.line)
                    self.fig.canvas.draw()
                    # Add a segment
                else:
                    # If there is a line, create a segment
                    x_data = [self.previous_point[0], x]
                    y_data = [self.previous_point[1], y]
                    logger.debug(f"draw line x: {x_data} y: {y_data}")
                    self.line = plt.Line2D(x_data, y_data, linewidth=1,
                                           marker='o', color=self.color)
                    self.previous_point = [x, y]
                    self.x.append(x)
                    self.y.append(y)

                    event.inaxes.add_line(self.line)
                    self.fig.canvas.draw()

            elif (((event.button == 1 and event.dblclick) or
                   (event.button == 3 and not event.dblclick)) and
                  self.line is not None):
                # Close the loop and disconnect
                logger.debug("Received single right mouse button click or "
                             "double left click")
                self.fig.canvas.mpl_disconnect(self.__cid1)
                self.fig.canvas.mpl_disconnect(self.__cid2)

                self.line.set_data([self.previous_point[0],
                                    self.start_point[0]],
                                   [self.previous_point[1],
                                    self.start_point[1]])
                ax.add_line(self.line)
                self.fig.canvas.draw()
                self.line = None
                self.lines = ax.get_lines()
                self.completed = True

                if not sys.flags.interactive and self.close_figure:
                    #  Figure has to be closed so that code can continue
                    plt.close(self.fig)

# BZ edit to make another box
class MultiRoi:
    def __init__(self, fig=None, ax=None, arr=None, ifg=None):
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
        self.rois   = {}
        self.npi    = {}
        self.txt    = self.ax.text(0.007, 1.005, ifg,
                                transform=self.ax.transAxes, fontsize=12)
        self.txt    = self.ax.text(0.45, 1.025, 'Click "New ROI" to begin',
                                transform=self.ax.transAxes, fontsize=14)

        self.make_buttons()
        return

    ## make two more buttons that add/subtract 2\pi
    def make_buttons(self):
        ax_add_btn = plt.axes([0.7, 0.02, 0.1, 0.04])
        ax_finish_btn = plt.axes([0.81, 0.02, 0.1, 0.04])

        ax_down_btn  = plt.axes([0.34, 0.02, 0.1, 0.04])
        ax_up_btn    = plt.axes([0.45, 0.02, 0.1, 0.04])
        ax_clear_btn = plt.axes([0.56, 0.02, 0.1, 0.04])

        btn_finish = Button(ax_finish_btn, 'Save & Close')
        btn_finish.on_clicked(self.finish)

        btn_add = Button(ax_add_btn, 'New ROI')
        btn_add.on_clicked(self.add)

        btn_up   = Button(ax_up_btn, 'Add $\\pi$')
        btn_up.on_clicked(self.increase)

        btn_down = Button(ax_down_btn, 'Subtract $\\pi$')
        btn_down.on_clicked(self.decrease)

        btn_clear = Button(ax_clear_btn, 'Clear')
        btn_clear.on_clicked(self.clear)

        plt.draw()
        plt.show(block=True)

    # store the text in self.rois dict?
    def add(self, event):
        """" Add a new ROI """

        # Only draw a new ROI if the previous one is completed
        if self.rois:
            if not all(r.completed for r in self.rois.values()):
                return

        roi_name = len(self.rois) # starts at 0, increments each new ROI

        logger.debug(f"Creating new ROI {roi_name}")

        ## clear the text
        self.txt.remove()
        self.txt  = self.ax.text(0.45, 1.025, 'Adjust # of $\\pi$ increments',
                                transform=self.ax.transAxes, fontsize=14)
        plt.draw()

        roi = RoiPoly(self.fig, self.ax, color='k', close_fig=False, show_fig=False)

        self.rois[roi_name] = roi
        self.npi[roi_name]  = 0
        # overwrite original array with new data; for multiple aois
        self.arr            = self.arr_up


    def finish(self, event):
        logger.debug("Stop ROI drawing")
        plt.close(self.fig)


    def increase(self, event):
        """ Increase integer (of 2pi cycles) to add """
        keys  = self.rois.keys()
        if len(keys) < 1:
            return
        else:
            idx  = max(keys)

        self.npi[idx] += 1

        ## npi to add
        npi             = self.npi[idx]
        try:
            mask_aoi    = self.rois[idx].get_mask(self.arr_up)
            self.arr_up = np.where(mask_aoi, self.arr+(npi*np.pi), self.arr)

            # update plot
            self.ax.get_images()[0].set_data(self.arr_up)
        except:
            logger.info('No data to actually update plot')

        logger.debug(f'AOI {idx}: +{self.npi[idx]}*2pi')

        try: self.txt.remove()
        except: pass

        ti       = f'Adding: {self.npi[idx]}$\\pi$'
        self.txt = self.ax.text(0.45, 1.025, ti, transform=self.ax.transAxes, fontsize=14)
        plt.draw()


    def decrease(self, event):
        """ Decrease integer (of 2pi cycles) to add """
        # get the index based on the roi
        keys  = self.rois.keys()
        if len(keys) < 1:
            return
        else:
            keys = [int(key) for key in keys]
            idx  = max(keys)
        self.npi[idx] -= 1

        ## npi to add
        npi         = self.npi[idx]
        try:
            mask_aoi    = self.rois[idx].get_mask(self.arr)
            self.arr_up = np.where(mask_aoi, self.arr+(npi*np.pi), self.arr)

            # update plot
            self.ax.get_images()[0].set_data(self.arr_up)
            # plt.draw()
        except:
            logger.info('No data to actually update plot')

        try: self.txt.remove()
        except: pass
        logger.debug(f'AOI {idx}: -{self.npi[idx]}*2pi')
        ti       = f'Adding: {self.npi[idx]}$\\pi$'
        self.txt = self.ax.text(0.45, 1.025, ti, transform=self.ax.transAxes, fontsize=14)
        plt.draw()


    def clear(self, event):
        self.txt.remove()
        self.txt  = self.ax.text(0.45, 1.025, 'Click "New ROI" to begin',
                                transform=self.ax.transAxes, fontsize=14)

        keys  = self.rois.keys()
        # if there isnt an aoi yet
        if len(keys) < 1:
            return
        # get the newest one
        else:
            idx  = max(keys)

        roi   = self.rois[idx]
        # no easy way to remove just the current ones (and while drawing)
        [line.remove() for line in roi.ax.get_lines()]

        npi             = self.npi[idx]
        ## also undo the addition/subtraction
        try:
            mask_aoi    = self.rois[idx].get_mask(self.arr)
            self.arr_up = np.where(mask_aoi, self.arr_up-(npi*np.pi), self.arr)
            self.ax.get_images()[0].set_data(self.arr_up)
        except:
            logger.info('No data to actually update plot')
        plt.draw()


        # reset the npi?
        self.rois.pop(idx)
        self.npi[idx] = 0
