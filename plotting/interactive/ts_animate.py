#!/usr/bin/env python3
import os
import re
import sys
import argparse
import numpy as np
from scipy import linalg, stats
from matplotlib import pyplot as plt, widgets, patches

from mintpy.objects import timeseries, giantTimeseries, HDFEOS
from mintpy.utils import arg_utils, ptime, time_func, readfile, utils as ut, plot as pp
from mintpy.multilook import multilook_data
from mintpy import subset, view, timeseries2velocity as ts2vel


def create_parser(subparsers=None):
    synopsis = 'Animate time series'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+',
                        help='time-series file to display\n'
                             'i.e.: timeseries_ERA5_ramp_demErr.h5 (MintPy)\n'
                             '      LS-PARAMS.h5 (GIAnT)\n'
                             '      S1_IW12_128_0593_0597_20141213_20180619.he5 (HDF-EOS5)')
    parser.add_argument('--label', dest='file_label', nargs='*',
                        help='labels to display for multiple input files')
    parser.add_argument('--ylim', dest='ylim', nargs=2, metavar=('YMIN', 'YMAX'), type=float,
                        help='Y limits for point plotting.')
    parser.add_argument('--tick-right', dest='tick_right', action='store_true',
                        help='set tick and tick label to the right')
    parser.add_argument('-l','--lookup', dest='lookup_file', type=str,
                        help='lookup table file')
    parser.add_argument('--no-show-img','--not-show-image', dest='disp_fig_img', action='store_false',
                        help='do NOT show the map figure.\n'
                             'Useful for plotting a point time series only.\n'
                             'This option requires --yx/lalo input.')

    parser.add_argument('-n', dest='idx', metavar='NUM', type=int,
                        help='Epoch/slice number for initial display.')
    parser.add_argument('--error', dest='error_file',
                        help='txt file with error for each date.')

    # time info
    parser.add_argument('--start-date', dest='start_date', type=str,
                        help='start date of displacement to display')
    parser.add_argument('--end-date', dest='end_date', type=str,
                        help='end date of displacement to display')
    parser.add_argument('--exclude', '--ex', dest='ex_date_list', nargs='*', default=['exclude_date.txt'],
                        help='Exclude date shown as gray.')
    parser.add_argument('--zf', '--zero-first', dest='zero_first', action='store_true',
                        help='Set displacement at first acquisition to zero.')
    parser.add_argument('--off','--offset', dest='offset', type=float,
                        help='Offset for each timeseries file.')

    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing.')

    # temporal model fitting
    parser.add_argument('--nomodel', '--nofit', dest='plot_model', action='store_false',
                        help='Do not plot the prediction of the time function (deformation model) fitting.')
    parser.add_argument('--plot-model-conf-int', '--plot-fit-conf-int', dest='plot_model_conf_int', action='store_true',
                        help='Plot the time function prediction confidence intervals.\n'
                             '[!-- Preliminary feature alert! --!]\n'
                             '[!-- This feature is NOT throughly checked. '
                             'Read the code before use. Interpret at your own risk! --!]')

    parser = arg_utils.add_timefunc_argument(parser)

    # pixel of interest
    pixel = parser.add_argument_group('Pixel Input')
    pixel.add_argument('--yx', type=int, metavar=('Y', 'X'), nargs=2,
                       help='initial pixel to plot in Y/X coord')
    pixel.add_argument('--lalo', type=float, metavar=('LAT', 'LON'), nargs=2,
                       help='initial pixel to plot in lat/lon coord')

    pixel.add_argument('--marker', type=str, default='o',
                       help='marker style (default: %(default)s).')
    pixel.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=6.0,
                       help='marker size (default: %(default)s).')
    pixel.add_argument('--lw', '--linewidth', dest='linewidth', type=float, default=0,
                       help='line width (default: %(default)s).')
    pixel.add_argument('--ew', '--edgewidth', dest='edge_width', type=float, default=1.0,
                       help='Edge width for the error bar (default: %(default)s)')

    # other groups
    parser = arg_utils.add_data_disp_argument(parser)
    parser = arg_utils.add_dem_argument(parser)
    parser = arg_utils.add_figure_argument(parser)
    parser = arg_utils.add_gps_argument(parser)
    parser = arg_utils.add_mask_argument(parser)
    parser = arg_utils.add_map_argument(parser)
    parser = arg_utils.add_memory_argument(parser)
    parser = arg_utils.add_reference_argument(parser)
    parser = arg_utils.add_save_argument(parser)
    parser = arg_utils.add_subset_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if inps.flip_lr or inps.flip_ud:
        inps.auto_flip = False

    if inps.gps_component:
        msg = '--gps-comp is not supported for {}'.format(os.path.basename(__file__))
        raise NotImplementedError(msg)

    if inps.file_label:
        if len(inps.file_label) != len(inps.file):
            raise Exception('input number of labels != number of files.')

    if (not inps.disp_fig or inps.outfile) and not inps.save_fig:
        inps.save_fig = True
    if inps.ylim:
        inps.ylim = sorted(inps.ylim)
    if inps.zero_mask:
        inps.mask_file = 'no'

    # default value
    inps.disp_unit = inps.disp_unit if inps.disp_unit else 'cm'
    inps.colormap = inps.colormap if inps.colormap else 'jet'
    inps.fig_size = inps.fig_size if inps.fig_size else [8.0, 4.5]

    # verbose print using --noverbose option
    global vprint
    vprint = print if inps.print_msg else lambda *args, **kwargs: None

    if not inps.disp_fig_img:
        if not inps.yx and not inps.lalo:
            inps.disp_fig_img = True
            print('WARNING: NO --yx/lalo input found for --no-show-img, turn it OFF and continue')

    if not inps.disp_fig:
        plt.switch_backend('Agg')

    return inps


###########################################################################################
def read_init_info(inps):
    # Time Series Info
    atr = readfile.read_attribute(inps.file[0])
    atr['DATA_TYPE'] = atr.get('DATA_TYPE', 'float32')

    inps.key = atr['FILE_TYPE']
    if inps.key == 'timeseries':
        obj = timeseries(inps.file[0])
    elif inps.key == 'giantTimeseries':
        obj = giantTimeseries(inps.file[0])
    elif inps.key == 'HDFEOS':
        obj = HDFEOS(inps.file[0])
    else:
        raise ValueError('input file is {}, not timeseries.'.format(inps.key))
    obj.open(print_msg=inps.print_msg)
    inps.seconds = atr.get('CENTER_LINE_UTC', 0)

    if not inps.file_label:
        inps.file_label = []
        for fname in inps.file:
            fbase = os.path.splitext(os.path.basename(fname))[0]
            fbase = fbase.replace('timeseries', '')
            inps.file_label.append(fbase)

    # default mask file
    if not inps.mask_file and 'msk' not in inps.file[0]:
        dir_name = os.path.dirname(inps.file[0])
        if 'Y_FIRST' in atr.keys():
            inps.mask_file = os.path.join(dir_name, 'geo_maskTempCoh.h5')
        else:
            inps.mask_file = os.path.join(dir_name, 'maskTempCoh.h5')
        if not os.path.isfile(inps.mask_file):
            inps.mask_file = None

    ## date info
    inps.date_list = obj.dateList
    inps.num_date = len(inps.date_list)
    if inps.start_date:
        inps.date_list = [i for i in inps.date_list if int(i) >= int(inps.start_date)]
    if inps.end_date:
        inps.date_list = [i for i in inps.date_list if int(i) <= int(inps.end_date)]
    inps.num_date = len(inps.date_list)
    inps.dates, inps.yearList = ptime.date_list2vector(inps.date_list)

    (inps.ex_date_list,
     inps.ex_dates,
     inps.ex_flag) = read_exclude_date(inps.ex_date_list, inps.date_list)

    # reference date/index
    if not inps.ref_date:
        inps.ref_date = atr.get('REF_DATE', None)
    if inps.ref_date:
        inps.ref_idx = inps.date_list.index(inps.ref_date)
    else:
        inps.ref_idx = None

    # date/index of interest for initial display
    if not inps.idx:
        if (not inps.ref_idx) or (inps.ref_idx < inps.num_date / 2.):
            inps.idx = inps.num_date - 2
        else:
            inps.idx = 2

    # Display Unit
    (inps.disp_unit,
     inps.unit_fac) = pp.scale_data2disp_unit(metadata=atr, disp_unit=inps.disp_unit)[1:3]

    # Read Error List
    inps.ts_plot_func = plot_ts_scatter
    inps.error_ts = None
    inps.ex_error_ts = None
    if inps.error_file:
        # assign plot function
        inps.ts_plot_func = plot_ts_errorbar

        # read error file
        error_fc = np.loadtxt(inps.error_file, dtype=bytes).astype(str)
        inps.error_ts = error_fc[:, 1].astype(np.float)*inps.unit_fac

        # update error file with exlcude date
        if inps.ex_date_list:
            e_ts = inps.error_ts[:]
            inps.ex_error_ts = e_ts[inps.ex_flag == 0]
            inps.error_ts = e_ts[inps.ex_flag == 1]

    # Zero displacement for 1st acquisition
    if inps.zero_first:
        inps.zero_idx = min(0, np.min(np.where(inps.ex_flag)[0]))

    # default lookup table file and coordinate object
    if not inps.lookup_file:
        inps.lookup_file = ut.get_lookup_file('./inputs/geometryRadar.h5')
    inps.coord = ut.coordinate(atr, inps.lookup_file)

    ## size and lalo info
    inps.pix_box, inps.geo_box = subset.subset_input_dict2box(vars(inps), atr)
    inps.pix_box = inps.coord.check_box_within_data_coverage(inps.pix_box)
    inps.geo_box = inps.coord.box_pixel2geo(inps.pix_box)
    data_box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))
    vprint('data   coverage in y/x: '+str(data_box))
    vprint('subset coverage in y/x: '+str(inps.pix_box))
    vprint('data   coverage in lat/lon: '+str(inps.coord.box_pixel2geo(data_box)))
    vprint('subset coverage in lat/lon: '+str(inps.geo_box))
    vprint('------------------------------------------------------------------------')

    # Map info - coordinate unit
    inps.coord_unit = atr.get('Y_UNIT', 'degrees').lower()
    inps = view.check_map_projection(inps, metadata=atr, print_msg=inps.print_msg)

    # calculate multilook_num
    # ONLY IF:
    #   inps.multilook is True (no --nomultilook input) AND
    #   inps.multilook_num ==1 (no --multilook-num input)
    # Note: inps.multilook is used for this check ONLY
    # Note: multilooking is only applied to the 3D data cubes and their related operations:
    # e.g. spatial indexing, referencing, etc. All the other variables are in the original grid
    # so that users get the same result as the non-multilooked version.
    if inps.multilook and inps.multilook_num == 1:
        inps.multilook_num = pp.auto_multilook_num(inps.pix_box, inps.num_date,
                                                   max_memory=inps.maxMemory,
                                                   print_msg=inps.print_msg)

    ## reference pixel
    if not inps.ref_lalo and 'REF_LAT' in atr.keys():
        inps.ref_lalo = (float(atr['REF_LAT']), float(atr['REF_LON']))
    if inps.ref_lalo:
        # set longitude to [-180, 180)
        if inps.coord_unit.lower().startswith('deg') and inps.ref_lalo[1] >= 180.:
            inps.ref_lalo[1] -= 360.
        # ref_lalo --> ref_yx if not set in cmd
        if not inps.ref_yx:
            inps.ref_yx = inps.coord.geo2radar(inps.ref_lalo[0], inps.ref_lalo[1], print_msg=False)[0:2]

    # use REF_Y/X if ref_yx not set in cmd
    if not inps.ref_yx and 'REF_Y' in atr.keys():
        inps.ref_yx = (int(atr['REF_Y']), int(atr['REF_X']))

    # ref_yx --> ref_lalo if in geo-coord
    # for plotting purpose only
    if inps.ref_yx and 'Y_FIRST' in atr.keys():
        inps.ref_lalo = inps.coord.radar2geo(inps.ref_yx[0], inps.ref_yx[1], print_msg=False)[0:2]

    # do not plot native reference point if it's out of the coverage due to subset
    if (inps.ref_yx and 'Y_FIRST' in atr.keys()
        and inps.ref_yx == (int(atr.get('REF_Y',-999)), int(atr.get('REF_X',-999)))
        and not (    inps.pix_box[0] <= inps.ref_yx[1] < inps.pix_box[2]
                 and inps.pix_box[1] <= inps.ref_yx[0] < inps.pix_box[3])):
        inps.disp_ref_pixel = False
        print('the native REF_Y/X is out of subset box, thus do not display')

    ## initial pixel coord
    if inps.lalo:
        inps.yx = inps.coord.geo2radar(inps.lalo[0], inps.lalo[1], print_msg=False)[0:2]
    try:
        inps.lalo = inps.coord.radar2geo(inps.yx[0], inps.yx[1], print_msg=False)[0:2]
    except:
        inps.lalo = None

    ## figure settings
    # Flip up-down / left-right
    if inps.auto_flip:
        inps.flip_lr, inps.flip_ud = pp.auto_flip_direction(atr, print_msg=inps.print_msg)

    # Transparency - Alpha
    if not inps.transparency:
        # Auto adjust transparency value when showing shaded relief DEM
        if inps.dem_file and inps.disp_dem_shade:
            inps.transparency = 0.7
        else:
            inps.transparency = 1.0

    ## display unit and wrap
    # if wrap_step == 2*np.pi (default value), set disp_unit_img = radian;
    # otherwise set disp_unit_img = disp_unit
    inps.disp_unit_img = inps.disp_unit
    if inps.wrap:
        inps.vlim = inps.wrap_range

        if (inps.wrap_range[1] - inps.wrap_range[0]) == 2*np.pi:
            inps.disp_unit_img = 'radian'

        if inps.disp_unit_img == 'radian':
            inps.range2phase = -4. * np.pi / float(atr['WAVELENGTH'])
            if   'cm' == inps.disp_unit.split('/')[0]:   inps.range2phase /= 100.
            elif 'mm' == inps.disp_unit.split('/')[0]:   inps.range2phase /= 1000.
            elif 'm'  == inps.disp_unit.split('/')[0]:   inps.range2phase /= 1.
            else:
                raise ValueError('un-recognized display unit: {}'.format(inps.disp_unit))

    inps.cbar_label = 'Amplitude' if atr['DATA_TYPE'].startswith('complex') else 'Displacement'
    inps.cbar_label += '[{}]'.format(inps.disp_unit_img)

    ## fit a suite of time func to the time series
    inps.model = time_func.inps2model(inps, date_list=inps.date_list, print_msg=inps.print_msg)

    # dense TS for plotting
    inps.date_list_fit = ptime.get_date_range(inps.date_list[0], inps.date_list[-1])
    inps.dates_fit = ptime.date_list2vector(inps.date_list_fit)[0]
    inps.G_fit = time_func.get_design_matrix4time_func(
        date_list=inps.date_list_fit,
        model=inps.model,
        seconds=inps.seconds)

    return inps, atr


def read_timeseries_data(inps):
    """Read data of time-series files
    Parameters: inps : Namespace of input arguments
    Returns:    ts_data : list of 3D np.array in size of (num_date, length, width)
                mask : 2D np.array in size of (length, width)
                inps : Namespace of input arguments
    """
    ## read list of 3D time-series
    ts_data = []
    for fname in inps.file:
        msg = f'reading timeseries from file {fname}'
        msg += f' with step of {inps.multilook_num} by {inps.multilook_num}' if inps.multilook_num > 1 else ''
        vprint(msg)
        data, atr = readfile.read(fname,
                                  datasetName=inps.date_list,
                                  box=inps.pix_box,
                                  xstep=inps.multilook_num,
                                  ystep=inps.multilook_num)
        if atr['DATA_TYPE'].startswith('complex'):
            vprint('input data is complex, calculate its amplitude and continue')
            data = np.abs(data)


        if inps.ref_idx is not None:
            vprint('reference to date: {}'.format(inps.date_list[inps.ref_idx]))
            data -= np.tile(data[inps.ref_idx, :, :], (inps.num_date, 1, 1))

        # Display Unit
        (data,
         inps.disp_unit,
         inps.unit_fac) = pp.scale_data2disp_unit(data,
                                                  metadata=atr,
                                                  disp_unit=inps.disp_unit)
        ts_data.append(data)

    ## mask file: input mask file + non-zero ts pixels - ref_point
    mask = pp.read_mask(inps.file[0],
                        mask_file=inps.mask_file,
                        datasetName='displacement',
                        box=inps.pix_box,
                        xstep=inps.multilook_num,
                        ystep=inps.multilook_num,
                        print_msg=inps.print_msg)[0]
    if mask is None:
        mask = np.ones(ts_data[0].shape[-2:], np.bool_)

    ts_stack = np.nansum(ts_data[0], axis=0)
    mask[np.isnan(ts_stack)] = False
    # keep all-zero value for unwrapError time-series
    if atr['UNIT'] not in ['cycle']:
        mask[ts_stack == 0.] = False
    del ts_stack

    # do not mask the reference point
    if inps.ref_yx and inps.ref_yx != (int(atr.get('REF_Y', -1)), int(atr.get('REF_X', -1))):
        (ry, rx) = subset_and_multilook_yx(inps.ref_yx, inps.pix_box, inps.multilook_num)
        mask[ry, rx] = True

    ## default vlim
    inps.dlim = [np.nanmin(ts_data[0]), np.nanmax(ts_data[0])]
    if not inps.vlim:
        inps.cmap_lut, inps.vlim = pp.auto_adjust_colormap_lut_and_disp_limit(ts_data[0],
                                                                              num_multilook=10,
                                                                              print_msg=inps.print_msg)[:2]
    vprint('data    range: {} {}'.format(inps.dlim, inps.disp_unit))
    vprint('display range: {} {}'.format(inps.vlim, inps.disp_unit))

    ## default ylim
    num_file = len(inps.file)
    if not inps.ylim:
        ts_data_mli = multilook_data(np.squeeze(ts_data[-1]), 4, 4)
        if inps.zero_first:
            ts_data_mli -= np.tile(ts_data_mli[inps.zero_idx, :, :], (inps.num_date, 1, 1))
        ymin, ymax = (np.nanmin(ts_data_mli[inps.ex_flag != 0]),
                      np.nanmax(ts_data_mli[inps.ex_flag != 0]))
        ybuffer = (ymax - ymin) * 0.05
        inps.ylim = [ymin - ybuffer, ymax + ybuffer]
        if inps.offset:
            inps.ylim[1] += inps.offset * (num_file - 1)
        del ts_data_mli

    return ts_data, mask, inps
