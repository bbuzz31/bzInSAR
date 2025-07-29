#!/usr/bin/env python

from __init__ import *

def createParser():
    parser = argparse.ArgumentParser(description='Make the snaphu commands for unwrapping multilooked ifgs.',
                                     epilog='Examples of use:\n\t prepUnwrap.py ./PS_DS --naz 2 --nrg 5',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_dir', type=str, help='Path to PS_DS (or a single file)')
    parser.add_argument('--net', type=str, default='SR', help='network label for writing to correct dir')
    parser.add_argument('-e', '--coarse', type=str,  default='', help='pass naz,nrg to use a multilooked ifg')
    parser.add_argument('--naz', default=7, type=int, help='Number of azimuth looks (default=7)')
    parser.add_argument('--nrg', default=19, type=int, help='Number of range looks (default=19)')
    parser.add_argument('--ntile', default=1, type=int, help='Number of tiles to use for unwrapping (default=1)')
    parser.add_argument('--ovl', default=500, type=int, help='Amount of pixel overlap for tiled unwrapping (default=500)')
    parser.add_argument('--extra',  default=' ',  choices={'', '.gf', '.resid', '.ps'},
                                    help='Use for additional filename descriptor')
    parser.add_argument('--update',  action='store_true',
                                help='Check if unwrapped files exist and skip if so')

    # parser.add_argument('--annual', action='store_true', help='Only do annual pairs')
    parser.add_argument('-o', '--dst', default=None, help='Optional path to write file'),
    parser.add_argument('-i', '--interact',  action='store_true', help='Enter an interactive session on completion')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(1)

    inps = parser.parse_args(args=iargs)
    return inps


def main(inps):
    """ Write the snaphu commands for multilooked files in input_dir (PS_DS)

    Also writes the commands for VRT post processing (so gdal can read them)
    """
    ## get the unwrapped ifgs
    if op.isdir(inps.input_dir):
        path_ifg = False
    else:
        path_ifg       = inps.input_dir
        inps.input_dir = op.dirname(op.abspath(inps.input_dir))

    os.chdir(inps.input_dir)

    ## filter multilooked/full res
    if int(inps.naz) > 1 or int(inps.nrg) > 1:
        ## filter filled/resid or not; do it this way or it can get too many files
        if inps.extra:
            pat  = f'20*{inps.naz}alks*{inps.extra}.int'
            ifgs = glob.glob(pat)
        else:
            pat  = f'20*{inps.nrg}rlks.int'
            ifgs = glob.glob(pat)

    else:
        if inps.extra:
            pat  = f'20*[0-9]{inps.extra}.int'
            ifgs = glob.glob(pat)

        else:
            pat = f'20*[0-9].int'
            ifgs = glob.glob(pat)

    if not ifgs:
        print (f'Could not find ifgs matching {pat} in {os.getcwd()}')
    # if inps.annual:
    #     ifgs = get_annual_pairs(ifgs, keep_sr=False)

    # for one ifg
    if path_ifg:
        ifgs = [ifg for ifg in ifgs if op.basename(path_ifg) in ifg]

    ifgs = sorted(ifgs)

    if 'dolphin' in inps.input_dir.lower():
        path_unw_dir = op.join(op.dirname(inps.input_dir), f'unwrap_ALL')
        path_corr0 = op.join(inps.input_dir, 'temporal_coherence.tif')
        path_corr  = op.join(inps.input_dir, 'temporal_coherence_ISCE.bin')
        gdal.Translate(path_corr, path_corr0, format='ISCE')
    else:
        path_unw_dir = op.join(op.dirname(inps.input_dir), 'ALL', f'unwrap_ALL')
        path_corr = op.join(inps.input_dir, 'tcorr_ds_ps.bin')



    if inps.update:
        ifgs = _filter_exist(ifgs, path_unw_dir)


    ## set up the label and get correlation (for later)
    # lbl = f'{op.basename(inps.input_dir)}'
    lbl = inps.net
    if int(inps.naz) > 1  or int(inps.nrg) > 1:
        ml        = f'.{inps.naz}alks_{inps.nrg}rlks'
        path_corr = op.join(inps.input_dir, f'tcorr_ds_ps_ISCE{ml}.bin')
        lbl       = f'{lbl}{ml}'

    lbl = f'{lbl}{inps.extra}'

    unw_cmds, isce_cmds = [],[]
    for wrap in ifgs:
        path_wrap = op.join(inps.input_dir, wrap)
        path_unw  = op.join(path_unw_dir, f'{op.splitext(wrap)[0]}.unw')
        path_conn = f'{path_unw}.conncomp'

        ## option to use the -e option with multilooked, upsampled, unwrapped phs
        coarse_opt = ''
        if inps.coarse:
            ifg      = f'{wrap[:17]}_ISCE'
            naz, nrg = inps.coarse.split(',')
            path_coarse = op.join(inps.input_dir, f'unwrap_{inps.net}'
                            f'{ifg}.{naz}alks_{nrg}rlks{inps.extra}.est.in')
            coarse_opt  = f'-e {path_coarse}'
            path_unw    = '{}c{}'.format(*op.splitext(path_unw))
            path_conn   = '{}c{}'.format(*op.splitext(path_conn))


        cmd = f'snaphu -f snaphu_{lbl}.conf {path_wrap} -o {path_unw} -g {path_conn} {coarse_opt} &&'
        unw_cmds.append(cmd)

        cmd = f'makeUnwVRT.py {path_unw} -w {path_wrap} &&'
        isce_cmds.append(cmd)

        cmd = f'makeUnwVRT.py {path_conn} -w {path_wrap} &&'
        isce_cmds.append(cmd)

    assert unw_cmds, 'Could not find any interferograms...'


    unw_cmds[-1]  = unw_cmds[-1].replace(' &&', '') # get rid of last && or bash yells
    isce_cmds[-1] = isce_cmds[-1].replace(' &&', '') # get rid of last && or bash yells

    wd       = op.dirname(op.realpath(inps.input_dir)) # where snaphu.conf is
    if inps.dst is None:
        dst  = op.join(wd, f'run_unwrap_{lbl}.sh')
    else:
        dst = inps.dst

    with open(dst, 'w') as fh:
        [print (cmd, file=fh) for cmd in unw_cmds]
    print ('Wrote snaphu commands to:', dst)

    dst      = op.join(wd, f'makeUnwVRT_{lbl}.sh')
    with open(dst, 'w') as fh:
        [print (cmd, file=fh) for cmd in isce_cmds]
    print ('Wrote vrt commands to:', dst)

    write_conf(path_corr, lbl, inps.ntile, inps.ovl)

    embed(colors='neutral') if inps.interact else ''
    print (f'Can now unwrap: {len(ifgs)} ifgs.')
    return


def write_conf(path_corr, lbl, ntile, ovl):
    """ Update the snaphu.conf file with the correct width and corr """
    assert op.exists(path_corr), f'Could not find correlation file: {path_corr}'
    width = gdal.Open(path_corr).ReadAsArray().shape[1]
    # parent of PS_DS
    dst     = op.join(op.dirname(op.dirname(path_corr)), f'snaphu_{lbl}.conf')
    content = f"""
########################
# Unwrapping parameters
#########################

STATCOSTMODE SMOOTH # DEFO
INITMETHOD MCF
VERBOSE FALSE

###############
# Input files
###############
LINELENGTH {width}

CORRFILE {path_corr}

################
# Output files
################

LOGFILE snaphu.log

################
# File formats
################

CORRFILEFORMAT FLOAT_DATA

################
# Tile control
################

NTILEROW {ntile} #10
NTILECOL {ntile} # 10

ROWOVRLP {ovl} # 500
COLOVRLP {ovl} # 500

NPROC {ntile}
MINREGIONSIZE 300

########################
# Other default parameters set by unwrap_fringe.py / isce
########################

ALTITUDE 700128.0009852069 # 800000.0
EARTHRADIUS 6358564.612282333 # 6371000.0
LAMBDA 0.056
NLOOKSRANGE 1
NLOOKSAZ 1
DEFOMAX_CYCLE 2.0
MAXNCOMPS 20 #100
"""
    with open(dst, 'w') as fh:
        fh.write(content)

    print ('Wrote:', dst)
    return


def _filter_exist(lst_wraps, path_unwraps):
    unws = os.listdir(path_unwraps)
    if not unws:
        print ('No unwrapped ifgs found, will write cmds to unwrap all')
        return

    ifgs_todo = []
    for ifg in lst_wraps:
        if not ifg.replace('.int', '.unw') in unws:
            ifgs_todo.append(ifg)
    print (f'Update mode, writing commands for {len(ifgs_todo)} of {len(lst_wraps)} ifgs.')
    return ifgs_todo


if __name__ == '__main__':
    inps = cmdLineParse()

    main(inps)
