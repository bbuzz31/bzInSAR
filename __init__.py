from VLM import *
from VLM.bzFRInGE.experiments import *
from osgeo import gdal, gdalconst
import re

warnings.filterwarnings("ignore", message="Mean of empty slice")
gdal.UseExceptions()

## hacks that dont actually work
os.sys.path.append(op.join(op.dirname(op.abspath(__file__)), 'multilook'))

BACKEND      = mpl.get_backend()

PATH_RES     = Path.home() / 'Desktop' if 'MT-' in \
                socket.gethostname() else Path.home() / 'VLM'

## probably should put tracks in here
DCT_ISCE     = {'HR': [36.75, 37.225, -76.655, -75.928],
                'NYC': [40.3341, 41.2099, -74.1259, -72.0931],
                # 'NYC1': [[40.35, 40.9, -74.55, -73.65], # this is what I actualyl want; need to re run isce
                'NNJ' : [40.3471, 40.9906, -74.4847, -73.9442],
                'NJ': [38.3667007, 41.8881242, -75.4336235, -71.626235],
                'PA+': [38.9182, 40.738, -76.6752, -73.8103],
                # 'DC': ,
                'Charleston':  [32.141118, 33.5, -80.5, -79.732874],
                'Savannah': [31.82, 32.265, -81.380640, -80.65],
                # 'Savannah': [31.879823, 32.262992, -81.380640, -80.836137],
                # 'Savannah': [31.379823, 32.762992, -82.880640, -80.336137],
                'Miami'   :  [25.0, 27.87, -80.56, -79.93],
                'Miami1'   : [25.375, 27.0, -80.565, -79.99],
                'Philly' : [39.25, 40.8, -76.0, -74.35],
                'Kiribati' : [1.6, 2.06, -157.51, -157.0],
                'SF' : [37.2244, 38.2736, -123.109, -121.6385],
                'Houston': [29.285256,  30.06908, -95.8574, -94.725740],
                'Kennedy': [28.222410, 28.988438, -81.998648, -80.294788],
                'Orlando': [28.222410, 28.988438, -81.998648, -80.294788], # symlink kennedy
                'DC': [38.78, 39.075, -77.22, -76.801],
                }

# reg to SNWE, track, frame
DCT_FRINGE  = {'HR':  [[36.75, 37.225, -76.655, -75.928], '004', ('115'), ['SPVA', 'VAHP']],
            'NYC': [[40.530, 40.822960, -74.200, -73.756498], '033', ('130'), ['NYBK']],
            # 'NNJ' : [[40.387010,  40.976408, -74.472474, -73.981022], '033', ('130'), ['NJHT']],
            # 'NNJ' : [[40.32181, 41.04098, -74.48468, -73.94440], '033', ('130'), ['NJHT']],

            'NNJ' : [[40.425, 40.9995, -74.48468, -73.971], '033', ('130'), ['NJHT']],
            'NJ' : [[38.3667007, 41.8881242, -75.4336235, -71.626235], '033', ('8621', '8622', '28221'), ['']],

            ## probably cant get baltimore so should just do this small and baltimore / DC on separate
            'PA+' : [[31.879823, 32.262992, -81.380640, -80.836137], '106', np.arange(121, 132, 1), ['?']],
            # 'DC' : [, '004', '120, 125', ['?']],

            # [40.3471, 40.9906, -74.4847, -73.9442]
            'Charleston':  [[32.5425, 32.9361, -80.2, -79.74], '150', ('102', '105'), ['SCHA']],
            # 'Savannah': [[31.879823, 32.262992, -81.380640, -80.836137], '150', ('102', '105'), ['SAVA']],
            'Savannah': [[31.82, 32.265, -81.380640, -80.65], '150', ('102', '105'), ['SAVA']],
            # 'Miami'   : [[25.01, 27.86, -80.55, -79.94], '048', np.arange(74, 89, 1), ['ZMA1']],
            'Miami1'  : [[25.375, 27.0, -80.5625, -80.0], '048', np.arange(74, 85, 1), ['ZMA1']], #frames prob off
            # 'Philly' : [[39.8438, 40.09006,  -75.3555, -75.00301], '106', ('123', '128'), ['FTML']],
            'Philly' : [[39.25, 40.8, -76.0, -74.35], '106', ('123', '128'), ['FTML']],
            'LongBeach' : [[33.49, 34.13, -118.59, -117.73], '106', ('064', '071'), ['LBCH']], # ARIA, Esther
            'Kiribati' : [[1.6, 2.06, -157.51, -157.0], '095', ('003', '1187'), ['KRTM']], # Pacific Islands
            'SF' : [[37.2244, 38.2736, -123.109, -121.6385], '035', ('117', '122'), ['SVIN']],
            'Houston': [[29.225,  30.06908, -95.8574, -94.725740], '034', ('90', '95'), ['NASA']],
            'Kennedy': [[28.282194898283, 28.745660961871, -80.927309541379, -80.48570961912], '048', ('89'), ['CCV6']],
            'Orlando': [[28.222410, 28.988438, -81.998648, -80.294788], '048', ('89'), ['FLOL']],
            'DC': [[38.77, 39.075, -77.22, -76.801], '004', [120, 121, 124, 125, 126], ['USN8']],
            }

DCT_FRINGE['NYC2'] = DCT_FRINGE['NYC']
DCT_FRINGE['Charleston2'] = DCT_FRINGE['Charleston']
DCT_FRINGE['Savannah2'] = DCT_FRINGE['Savannah']
DCT_FRINGE['Miami'] = DCT_FRINGE['Miami1']

# subset NYC during fringe/ts processing
DCT_REG = DCT_FRINGE

# should replace this with a function to fetch from the UNR files
# lat / lon
DCT_GPS = {'VAHP': '37.062178 -76.403436',
           'LOY2': '36.764016798 -76.237800713',
           'LOYZ' : '36.863597 -76.573562', # -76.6257, move out of water for ARIA
           'LS03' : '36.78868 -75.95924', # one pixel E to coherent
           'SPVA' : '36.94225 -76.32873', # one pixel S to put on land
           'DRV5' : '36.9586563477 -76.5566424163',
           'SCHA' : '32.779477 -79.925365',
           'NYBP' : '40.7021 -74.0143231199', #40.7010691908; actual loc on TG is in water
           'NJHT' : '40.731134 -74.0372', # -74.037577 one pixel right for coh
           'NJI2' : '40.7414064076 -74.1777468016',
           'NYPR' : '40.6378535635 -74.1248108346',
           'NYJM' : '40.6628023005 -73.8070851173',
           'NYBR' : '40.6886607623 -74.0012782594',
           'NYOB' : '40.5515520655 -74.1165786135',
           'NYBK' : '40.7034315045 -73.9789651762',
           'SAVA': '32.078716  -81.089722',
           'ZMA1': '25.8246120479 -80.3191886896',
           'LBCH' : '33.7877726 -118.203350', # Long Beach
           'FTML' : '39.876033 -75.210490', # Philly
           'KRTM' : '2.046560 -157.4476305', # Kirbati
           'SVIN' : '-122.526317813 38.0331786882', # San Rafael
           'NASA' : '29.5519534 -95.0962194', # Houston
           'UHDT' : '29.76596 -95.35944', # Houston
           'CFJV' : '29.88165 -95.55584', # Houston
           'SESG' : '29.98747 -95.42962', # Houston
           'TXTG' : '29.89752 -95.29738', # Houston
           'CSTE': '29.795637 -95.510738', # Houston
           'DMFB': '29.62265 -95.58374', # Houston
           'TXLQ': '29.35796 -94.95285', # Houston
           'TXB6': '29.7569098 -94.9373546', # Houston
           'CCV6': '28.4599573567 -80.5454762056', # Kennedy
           'FLOL': '28.5709739532 -81.4236549764', # Orlando
           'USN7': '38.920566 -77.06615', # DC; moved slightly east to
           'USN8': '38.920566 -77.06615', # DC; moved slightly east to
           'NRL1': '38.8207469933 -77.0244027790',
           }

# GPS station to manually drawn crop
CROP_MAP = {'LOY2': 'HR_South', 'VAHP': 'HR_North', 'LOYZ': 'HR_Southwest',
            'LS03': 'HR_South', 'SPVA': 'HR_South',
            'NJHT': 'NYC_West', 'NYJM': 'NYC_East1', 'NYBK': 'NYC_East1',
            }

DCT_EXCL = {'HR': ['KPVA', 'HB15', 'EVS4', 'DRV6'], # 'DRV5/6'
            'Charleston': [],
            'NYC':'NYBP NYJM ROG1'.split(),
            'NNJ':'NYBP NYJM ROG1 NJMT NJTP NYBK NYBR NJDY'.split(), # last 5 outside of data
            # 'NYC':'NYJM ROG1'.split(),
            'Savannah': [],
            'Houston': 'ANG5 ANG6 COH1 COH6 COTM DISD OKEK DEN2 GAL2 GAL7 TXHU TXKY TXLI TXP5 TXRN TXRO TXRS UH01 UHCO UHKD UHKS'.split(),
            'DC': 'GODE GODN GODS GODZ SA02 S071'.split() # outside InSAR; S071 = old and bad
            }


## ifg utils
def get_annual_pairs(ifgs, months=None):
    """ Takes only the annual pairs from list of igs. Called from get_pairs.

    Specify a string of months in the exp dictionar e.g. "12 1 2" to use just those
    """
    ann_pairs = []
    for ifg in ifgs:
        path, basen = op.dirname(ifg), op.basename(ifg)
        dt          = basen.split('_')[:2]
        ref, sec    = [datetime.strptime(dti, '%Y%m%d') for dti in dt]
        elap        = (sec-ref).days

        if not elap > 364.25:
            continue

        if isinstance(months, str):
            lst_months = [int(m) for m in months.split('-')]
            if ref.month in lst_months:
                ann_pairs.append(op.join(path, basen))
        else:
            ann_pairs.append(op.join(path, basen))

    return ann_pairs


def get_sbas_pairs(ifgs, nconn, max_tbase=36):
    """ Take nconnections of ifgs with temporal baseline <= max_tbase (days) """
    import pandas as pd
    refs, secs, pairs = [], [], []
    for ifg in ifgs:
        path, basen = op.dirname(ifg), op.basename(ifg)
        dt          = basen.split('_')[:2]
        ref, sec    = [datetime.strptime(dti, '%Y%m%d') for dti in dt]
        elap        = (sec-ref).days

        if elap > max_tbase:
            continue

        refs.append(ref)
        secs.append(sec)
        pairs.append(op.join(path, basen))

    df = pd.DataFrame({'ref': refs, 'sec': secs, 'pairs': pairs})

    ## take first nconn occurances; note that max_tbase takes care of some of this
    df = df.sort_values('ref sec'.split()).groupby('ref').head(nconn)
    return df.pairs.tolist()


def get_ifgs(path_ifgs, dct_exp, ext='int', max_tbase=36):
    """ Get all the ifg pairs in a directory """
    ## take only gap filled or not gap filled
    naz, nrg = dct_exp['naz'], dct_exp['nrg']
    ext      = f'gf.{ext}' if dct_exp['gap_fill'] else ext

    search = f'20*{naz}alks_{nrg}rlks.{ext}' if (naz > 1 or nrg > 1) else f'20*[0-9].{ext}'
    log.info ('Search:', search)
    ifgs_all = sorted(list(path_ifgs.glob(search)))

    ## take annual and/or sbas pairs if specified
    ifgs_ann = get_annual_pairs(ifgs_all, dct_exp['ann']) if \
                                                dct_exp['ann'] else []

    ifgs_sbas = get_sbas_pairs(ifgs_all, dct_exp['sbas'], max_tbase) if \
                                                dct_exp['sbas'] else []
    if not ifgs_ann and not ifgs_sbas:
        ifgs = sorted(ifgs_all)
    else:
        ifgs = sorted(np.unique(ifgs_ann + ifgs_sbas))

    return ifgs


def edit_runfiles(path_exp, n=2, nprocs=16):
    """
    Since you specified --connections all in stackSentinel need to reduce n ifgs
    """
    root = op.join(path_exp, 'run_files')
    for f in 'run_13_generate_burst_igram run_14_merge_burst_igram run_15_filter_coherence'.split():
        orig = op.join(root, f)
        backup = f'{orig}_all'
        if not op.exists(backup):
            shutil.copy(orig, backup)

        dst  = f'{orig}_new'
        with open(dst, 'w') as fhw:
            with open(backup, 'r') as fhr:
                ref_prev, i, lineno = '', 0, 0
                for line in fhr:
                    if not line.strip() or line.startswith('wait'):
                        continue

                    # get the new ref
                    ref = line.split()[2][-17:-9]
                    # if the ref exists
                    if ref == ref_prev:
                        # if connections not satisfied
                        if i < n:
                            fhw.write(line)
                            i+=1
                            lineno+=1
                        else:
                            continue

                    else:
                        fhw.write(line)
                        lineno+=1
                        i=1

                    ref_prev = ref

                    # write 'wait' and a blank line so that group of nproc finishes
                    if lineno > 0 and (lineno%nprocs)==0:
                        fhw.write('wait\n\n')
        os.rename(dst, orig)
        print ('Overwrote:', orig)
    return


def get_unw_dirs(path_ps_ds):
    """ Get all unwrap directories (dirs that have 'unwrap' in them) """
    unwrap_dirs = []
    for root, dirs, files in os.walk(path_ps_ds):
        for d in dirs:
            if 'unwrap' in d:
                unwrap_dirs.append(op.join(root, d))
    return unwrap_dirs


## mask utils
def combine_masks(Exp, lst_path_masks, apply_osm=False, mask_all=True, show=False):
    """ Used to combine multiple custom Vup masks / ARIA temporal coherence

    If mask all, pixel is masked if masked in any mask. (FRInGE)
    Else, if pixel is not masked in any, keeps it. (ARIA)
    """
    from mintpy.utils import readfile, writefile
    import h5py, shutil, xarray as xr
    lst_masks = []
    for mask in lst_path_masks:
        arr, meta = readfile.read(mask)
        lst_masks.append(arr)


    ## reapply the watermask to account for some geocoding error
    if apply_osm:
        # open the netcdf watermask for the correct size/shape
        da0 = xr.open_dataset(op.join(Exp.path_mp_exp, 'geo', 'geo_waterMask.nc'))['Band1']
        # open the full watermask
        path_osm = op.join(Exp.path_crops, 'OSM_wmask.tif')
        da = xr.open_dataset(path_osm)['band_data'].squeeze().rename(x='lon', y='lat')

        # crop the full mask
        S, N, W, E = da0.lat.min(), da0.lat.max(), da0.lon.min(), da0.lon.max()
        da_crop    = da.sel(lat=slice(N, S), lon=slice(W, E))
        da_crop    = da_crop.where(~np.isclose(da_crop, 0), np.nan)
        da_mask1   = da_crop.interp_like(da0)
        da_mask1   = da_mask1.where(~np.isnan(da_mask1), 0)
        lst_masks.append(np.flipud(da_mask1.data))

    print ('Using', len(lst_masks), 'masks')

    mask = np.stack(lst_masks)
    if mask_all:
        mask = mask.prod(0, dtype=bool)
    else:
        mask = mask.sum(0, dtype=bool)

    # copy any mask h5 to where to write; use stitched if stitch in any of the masks
    # stitch = True if any(['stitch' in p for p in lst_path_masks]) else False
    stitch   = True if 'stitch' in Exp.mp_exp.lower() else False
    dst = Exp.path_mask_vup_stitch if stitch else Exp.path_mask_vup
    dst = shutil.copy(lst_path_masks[0], dst)

    try:
        writefile.write_hdf5_block(dst, mask, 'waterMask')
    except:
        writefile.write_hdf5_block(dst, mask, 'mask')

    print ('Wrote combined mask to:', dst)

    if show:
        from mintpy.cli import view
        cmd = f'{dst} --noverbose -c binary'
        obj = view.main(cmd.split())

    return dst


## GPS utils; probably to a class
def prep_gps(path_midas, reg, units='m', horiz=False, exclude=True):
    """ Get GPS stas within SNWE. Exclude stations in DCT_EXCL by default """

    df     = pd.read_csv(path_midas)

    keep = 'sta lat lon u_vel u_sig ts_start ts_end'.split()

    if 'u_vel_orig' in df.columns:
        log.warning('Using updated MIDAS rate/unc')
        keep += 'u_vel_lsq u_sig_lsq u_vel_orig u_sig_orig ts_start_orig ts_end_orig'.split()
    if horiz:
        keep += 'e_vel e_sig n_vel n_sig'.split()

    df = df[keep]

    if exclude:
        try:
            df = df[~df.sta.isin(DCT_EXCL[reg])]
        except:
            df = df

    # native units are meters
    cols  = [col for col in df.columns if 'vel' in col or 'sig' in col]
    df    = df.copy()
    if units == 'mm':
        df[cols] *= 1000
    elif units == 'cm':
        df[cols] *= 100

    return df


def gps_rmse_npix_nc(Exp, npix=0, exclude=[], df_gps=None, mask=True):
    """ Calculate the RMSE between each station in insar rates and df_gps within radius (meters)
        - this is in the CompareExps Notebook
        - rate_nc is an opened dataset (xr.open_dataset(Exp.path_rate_msk_nc)['Band1'])
        - df_gps is gotten with prep_gps(Exp.path_gps, Exp.reg)
        used to be called gps_rmse_local
    """
    from BZ import bbGIS
    print ('Dropping stations:', ','.join(exclude)) if exclude else ''
    df_gps = prep_gps(Exp.path_gps, Exp.reg) if df_gps is None else df_gps
    df_gps = df_gps[~df_gps.sta.isin(exclude)]
    df_gps = bbGIS.in_dfll(df_gps, SNWE=Exp.SNWE)
    dct_cmp = {}
    rate_nc, std_nc = Exp.get_rate_unc_nc(mask=mask)
    for ix, row in df_gps.iterrows():
        lat, lon = row['lat lon'.split()]
        da_vel = bbGIS.get_surrounding_pixels(rate_nc, lat, lon, npix=npix)
        da_unc = bbGIS.get_surrounding_pixels(std_nc, lat, lon, npix=npix)
        arr_vel, arr_unc = np.nanmean(da_vel), np.nanmean(da_unc)
        dct_cmp[row.sta] = [row.u_vel*1000, row.u_sig*1000, arr_vel, arr_unc]
        # log.info(f'{row.sta}: InSAR / GPS: {arr_vel:.2f} /' \
        #           f'{1000*row.u_vel:.2f} mm/yr')

    df = pd.DataFrame(dct_cmp).T / 1000 # units are now m
    df.columns='gps_vel gps_unc ins_vel ins_unc'.split()
    df['resid'] = df['gps_vel'] - df['ins_vel']
    n0 = df.shape[0]
    df = df.dropna(subset=['resid'])
    log.warning (f'Dropped {n0-df.shape[0]} of {n0} GPS stations without nearby InSAR.')
    return df


def gps_rmse_npix(Exp, npix=0, exclude=[], df_gps=None, mask=True, verbose=True):
    """ Same as GPS RMSE local except use the pixel coordsinates """
    from BZ import bbGIS
    from mintpy.utils import readfile, utils as ut
    log.setLevel(10) if verbose else log.setLevel(20)
    exclude = [exclude] if isinstance(exclude, str) else exclude
    log.info ('Dropping stations: %s', ','.join(exclude)) if any(exclude) else ''
    df_gps = prep_gps(Exp.path_gps, Exp.reg) if df_gps is None else df_gps
    df_gps = df_gps[~df_gps.sta.isin(exclude)]
    df_gps = bbGIS.in_dfll(df_gps, SNWE=Exp.SNWE)

    arr, arr_unc, mask, meta = Exp.get_rate_unc_h5(mask)
    # stitched  = True if ('stitch' in Exp.mp_exp.lower() and
    #                      op.exists(Exp.path_mask_vup_stitch)) else False
    # path_vup  = Exp.path_vup_geo if not stitched else Exp.path_vup_geo_stitch
    # path_mask = Exp.path_mask_vup if not stitched else Exp.path_mask_vup_stitch

    # arr, meta = readfile.read(path_vup)
    # arr_unc   = readfile.read(path_vup, datasetName='velocityStd')[0]

    # should still support stitched
    coord = ut.coordinate(meta, lookup_file=Exp.path_geom_mp_geo)
    dct_cmp   = {}
    for ix, row in df_gps.iterrows():
        y, x     = coord.geo2radar(row.lat, row.lon)[0:2]
        try:
            if npix:
                arr_vel = np.nanmean(arr[y-npix:y+npix, x-npix:x+npix])
                arr_unc1 = np.nanmean(arr_unc[y-npix:y+npix, x-npix:x+npix])
                if np.isnan(arr_vel):
                    log.warning('%s is masked despite averaging, skipping...', row.sta)
                    continue
            else:
                if isinstance(mask, np.ndarray):
                    if not mask[y,x]:
                        log.warning('%s is masked, skipping...', row.sta)
                        continue
                arr_vel = arr[y,x]
                arr_unc1 = arr_unc[y,x]

        except:
            log.warning('%s is not covered by InSAR, skipping...', row.sta)
            continue

        # log.debug(f'{row.sta}: InSAR / GPS: {1000*arr_vel:.2f} /' \
        #           f'{1000*row.u_vel:.2f} mm/yr')
        dct_cmp[row.sta] = [row.u_vel, row.u_sig, arr_vel, arr_unc1]


    df = pd.DataFrame(dct_cmp).T
    df.columns='gps_vel gps_unc ins_vel ins_unc'.split()
    df['resid'] = df['gps_vel'] - df['ins_vel']
    n0 = df.shape[0]
    df = df.dropna(subset=['resid'])
    log.warning (f'Dropped {n0-df.shape[0]} of {n0} GPS stations without nearby InSAR.')
    return df


def calc_r2(actual_values, predicted_values):
    # Calculate the mean of actual values
    mean_actual = np.mean(actual_values)

    # Calculate the total sum of squares
    total_sum_squares = np.sum((actual_values - mean_actual)**2)

    # Calculate the residual sum of squares
    residual_sum_squares = np.sum((actual_values - predicted_values)**2)

    # Calculate R²
    r_squared = 1 - (residual_sum_squares / total_sum_squares)

    return r_squared


def get_insar_uncertainty(Exp1, Exp2=None):
    """ Get the InSAR only vertical uncertainty """
    from mintpy.utils import readfile
    if Exp2 is None:
        arr_unc, meta = readfile.read(Exp1.path_vlos_geo, datasetName='velocityStd')
        ## get the InSAR only uncertainty
        inc       = readfile.read(Exp1.path_geom_mp_geo, datasetName='incidenceAngle')[0]
        arr_unc  /= np.cos(np.deg2rad(inc)) # bumps it up slightly
        return arr_unc

    ## for stitched; untested
    dct_ref  = {'NJHT': 'NYC_West', 'NYJM': 'NYC_East', 'NYBK': 'NYC_East',
                'NYBR': 'NYC_West'}
    epsg     = 'epsg:4326'

    lst_incs = []
    for Expi in [Exp1, Exp2]:
        da_unc  = h5_to_xr(Expi.path_vlos, 'velocityStd')
        da_unc.rio.write_crs(epsg, inplace=True)

        path_crop    = op.join(Expi.path_crops, f'{dct_ref[Expi.ref_sta]}.GeoJSON')
        gdf_crop     = gpd.read_file(path_crop)
        da_unc_crop  = da_unc.rio.clip(gdf_crop.geometry, epsg, drop=False)
        lst_uncs.append(da_unc_crop)

    unc_stitch  = np.flipud(np.where(np.isnan(lst_uncs[0].data),
                                     lst_uncs[1].data, lst_uncs[0].data))

    return unc_stitch


## GIS utils
def get_extent(path_ds, shrink=None):
    """ Get the WESN extent of some geocoded file (needs to be opened with gdal) """
    from osgeo import gdal
    if isinstance(path_ds, str):
        ds   = gdal.Open(path_ds, gdal.GA_ReadOnly)
    else:
        ds   = path_ds
    trans    = ds.GetGeoTransform()
    # W E S N
    extent   = [trans[0], trans[0] + ds.RasterXSize * trans[1],
                trans[3] + ds.RasterYSize*trans[5], trans[3]]
    if shrink is not None:
        delW, delE, delS, delN = shrink
        extent = [extent[0] + delW, extent[1] - delE, extent[2] + delS, extent[3] - delN]

    del ds
    return extent


def avg_ts(ts0, lalo, rad, path_vel):
    """
    Create a ts of InSAR values 'rad' around lalo (len 2 list)

    ts0 is the data cube with time in first dimension
    """
    from BZ import bbGIS
    assert ts0.ndim in [2, 3], 'Input array must be two or 3 dimensions'

    # get the geocoded velocity for use in circle
    ds_circle = bbGIS.make_circle_rast(lalo, rad, path_vel)
    arr  = ds_circle.ReadAsArray()
    idx  = np.where(arr.reshape(-1) < 9990)[0]

    len_ts = ts0.shape[0]

    ## get disp pixels within radius
    ts_2D         = ts0.reshape(len_ts, -1)[:, idx] if ts0.ndim == 3 else ts0[:, idx]
    masked_pixels = np.isnan(ts_2D).all(axis=0) # mask time series of all nan
    ts_2D         = ts_2D[:, ~masked_pixels]

    if not ts_2D.shape[1] > 0:
        log.info('No valid pixels surrounding %s', lalo)
        return np.full(len_ts, np.nan)

    log.info('Averaging %s pixels', ts_2D.shape[1])
    ts = np.nanmean(ts_2D, axis=1)
    return ts


def buffer_point(lat, lon, rad):
    """ Roughly buffer a lalo point given a radius of meters

    Return S, N, W, E bounding box in WGS84
    """
    from shapely.geometry import Point
    import geopandas as gpd
    gser = gpd.GeoSeries(Point(lon, lat), crs='EPSG:4326')
    gdf  = gser.to_crs('EPSG:4087').buffer(rad).to_crs('EPSG:4326')
    W, S, E, N = gdf.squeeze().bounds
    return S, N, W, E


def vlm_at_points(Exp, df_pts0=None, npix=0, mask=True, verbose=True):
    """ Get the InSAR values in an npix box each point location in 'df_pts' """
    from BZ import bbGIS
    from mintpy.utils import readfile, utils as ut

    log.setLevel(10) if verbose else log.setLevel(20)
    df_pts = bbGIS.in_dfll(df_pts0, SNWE=Exp.SNWE)
    n0, n1 = df_pts0.shape[0], df_pts.shape[0]
    if n0 != n1:
        log.warning(f'Dropped {n0-n1} of {n0} locations outside of SNWE bbox.')

    arr, arr_unc, mask, meta = Exp.get_rate_unc_h5(mask)
    # stitched  = True if ('stitch' in Exp.mp_exp.lower() and
    #                      op.exists(Exp.path_mask_vup_stitch)) else False
    # path_vup  = Exp.path_vup_geo if not stitched else Exp.path_vup_geo_stitch
    # path_mask = Exp.path_mask_vup if not stitched else Exp.path_mask_vup_stitch

    # arr, meta = readfile.read(path_vup)
    # arr_unc   = readfile.read(path_vup, datasetName='velocityStd')[0]

    # should still support stitched
    coord = ut.coordinate(meta, lookup_file=Exp.path_geom_mp_geo)
    dct_res = {}
    for ix, row in df_pts.iterrows():
        y, x     = coord.geo2radar(row.lat, row.lon)[0:2]
        try:
            if npix:
                arr_vel = np.nanmean(arr[y-npix:y+npix, x-npix:x+npix])
                arr_unc1 = np.nanmean(arr_unc[y-npix:y+npix, x-npix:x+npix])
                if np.isnan(arr_vel):
                    log.warning('%s is masked despite averaging, skipping...', row.name)
            else:
                if isinstance(mask, np.ndarray):
                    if not mask[y,x]:
                        log.warning('%s is masked, skipping...', row.name)
                        continue
                arr_vel = arr[y,x]
                arr_unc1 = arr_unc[y,x]

        except:
            log.warning('%s is not covered by InSAR, skipping...', row.name)
            continue

        dct_res[row.name] = arr_vel*1000, arr_unc1*1000, row.lat, row.lon

    df = pd.DataFrame(dct_res).T
    df.columns='ins_vel ins_unc lat lon'.split()
    n0 = df.shape[0]
    df = df.dropna()
    log.warning (f'Dropped {n0-df.shape[0]} of {n0} sites stations without InSAR within npix={npix}.')
    return df


def crop_utm(path_h5, layer, SNWE, dst=None, epsg_src=32615):
    """
    Crops an HDF5 dataset to a specified SNWE bounding box.

    Parameters:
        hdf5_file (str): Path to the input HDF5 file.
        dataset_name (str): Name of the dataset to crop.
        output_file (str): Path to save the cropped HDF5 file.
        snwe_bbox (tuple): Bounding box (S, N, W, E) in dataset coordinates (e.g., UTM meters).

    Returns:
        None
    """
    from mintpy.utils import readfile
    from shapely.geometry import box
    import geopandas as gpd

    S, N, W, E = SNWE
    gser_bbox = gpd.GeoSeries(box(W, S, E, N), crs=4326).to_crs(epsg_src)


    arr, atr = readfile.read(path_h5, datasetName=layer)
    x_first = float(atr['X_FIRST'])
    x_step = float(atr['X_STEP'])
    y_first = float(atr['Y_FIRST'])
    y_step = float(atr['Y_STEP'])

    # x_start_idx = int((W - x_first) / x_step) + x_first
    # x_end_idx = int((E - x_first) / x_step) + x_first + 1
    # y_start_idx = int((N - y_first) / abs(y_step)) + y_first
    # y_end_idx = int((S - y_first) / abs(y_step)) + y_first + 1

    # Calculate column indices
    col_start = int(np.floor((W - x_first) / x_step))
    col_end = int(np.ceil((E - x_first) / x_step))

    # Calculate row indices
    row_start = int(np.floor((y_first - N) / abs(y_step)))
    row_end = int(np.ceil((y_first - S) / abs(y_step)))

    # Clip indices to dataset bounds
    col_start = max(col_start, 0)
    col_end = min(col_end, arr.shape[2])
    row_start = max(row_start, 0)
    row_end = min(row_end, arr.shape[1])

    print(f"Cropping indices:")
    print(f"  Rows: {row_start} to {row_end}")
    print(f"  Columns: {col_start} to {col_end}")

    # Extract the cropped data for all bands
    arr_crop = arr[:, row_start:row_end, col_start:col_end]

    # dst = dst if dst else f'{op.splitext(path_h5)[0]}_crop.h5'
    dst = dst / f'{path_h5.stem}_crop.h5'


    breakpoint()
    # Save the cropped dataset to a new HDF5 file
    with h5py.File(dst, 'w') as h5:
        dset = h5.create_dataset(layer, data=arr_crop, compression="gzip")

        # Update attributes for the cropped dataset
        h5.attrs['X_FIRST'] = x_first + x_start_idx * x_step
        h5.attrs['Y_FIRST'] = y_first - y_start_idx * abs(y_step)
        h5.attrs['X_STEP'] = x_step
        h5.attrs['Y_STEP'] = y_step
        for k, v in atr.items():
            if not k in h5.attrs:
                h5.attrs[k] = v

    print(f'Cropped dataset saved to {dst}')
    return



## Load Utils
def load_exp(fr_exp, mp_exp='Base', gps_exp='Default'):
    if gps_exp == 'Default':
        ref_sta = DCT_REG[fr_exp['root'].split('_')[0]][3][0]
    else:
        ref_sta = gps_exp

    Exp     = ExpBase(fr_exp, mp_exp, gps_exp)
    return Exp


def h5_to_xr(file, dset, ref_date=False):
    """ Updated from save_gdal to use rasterio and handle timeseries """
    from mintpy.utils import readfile
    import rasterio as rio
    import xarray as xr
    from rasterio.crs import CRS as rCRS
    from rasterio.transform import Affine
    import time

    ## read data
    ftype = readfile.read_attribute(file)['FILE_TYPE']

    # grab ref_date from dset
    if ftype == 'timeseries' and dset and '_' in dset:
        ref_date, dset = dset.split('_')
    else:
        ref_date = None

    ds_name = dset if dset else 'data'
    data, atr = readfile.read(file, datasetName=dset)

    if ftype == 'timeseries' and ref_date:
        print(f'read {ref_date} from file: {file}')
        data -= readfile.read(file, datasetName=ref_date)[0]

    transform = Affine.translation(float(atr['X_FIRST']),
                                   float(atr['Y_FIRST'])) * \
                                   Affine.scale(float(atr['X_STEP']),
                                                float(atr['Y_STEP']))

    crs = rCRS.from_epsg(atr.get('EPSG', '4326'))

    data = data[np.newaxis, :, :] if data.ndim == 2 else data
    bands, rows, cols = data.shape

    profile = {
        'driver': 'netcdf',
        'height': rows,
        'width': cols,
        'count': bands,
        'dtype': np.float32,
        'crs': crs,
        'transform': transform,
        'compress': 'lzw'  # Optional: compression
    }

    tmp_path = Path(os.getenv('dataroot')) / 'temp.nc'
    os.remove(tmp_path) if Path.exists(tmp_path) else ''

    with rio.open(tmp_path, 'w', **profile) as ds:
        ds.write(data)
        ds.close()

    # while not Path.exists(tmp_path):
        # log.info('Waiting for netcdf file to be created...')
        # time.sleep(60*5)
        # ds.write(data)
        # ds.close()

    ds = xr.open_dataset(tmp_path)
    if bands == 1:
        da = ds['Band1']
    else:
        # get all the variables to concat
        vars = list(ds.drop('crs').keys())
        da = xr.concat([ds[var] for var in vars], dim='time')

    os.remove(tmp_path)
    return da.rename(dset)


def load_ts_nc(exp, rm_ann=True, overwrite=True, path_ts_geo=None):
    """ Convert / load the h5 timeseries into netcdf. Default to apply mask. """
    from mintpy.utils import readfile
    from mintpy.objects import timeseries
    path_ts_geo = path_ts_geo if path_ts_geo else exp.path_ts_geo
    dst = f'{op.splitext(path_ts_geo)[0]}.nc'
    if op.exists(dst) and not overwrite:
        da = xr.open_dataset(dst)['timeseries']
        print ('Got:', dst)
        return da

    if rm_ann:
        log.debug('Using timeseries with annual cycle removed')
        tsObj = timeseries(exp.path_ts_geo_ann) # what we should be doing
    else:
        tsObj = timeseries(path_ts_geo)

    arr_ts = tsObj.read()
    dates = tsObj.get_date_list()

    ## next 3 are temporary cuz not using annual (which already excluded)
    # date_ix = [i for i, dt in enumerate(dates) if not dt in get_excluded_dates(exp).split()]
    # dates = np.array(dates)[date_ix]
    # arr_ts = arr_ts[date_ix]

    with h5py.File(exp.path_mask_mp_geo, 'r') as h5:
        try:
            mask    = h5['waterMask'][:]
        except:
            mask    = h5['mask'][:]
        arr_ts *= np.where(np.isclose(mask, 0), np.nan, 1)

    attrs = readfile.read_attribute(exp.path_ts_geo_ann)
    if exp.reg == 'Kiribati' or attrs['ORBIT_DIRECTION'] == 'ASCENDING':
        arr_ts = np.fliplr(arr_ts)

    ## for matching the geom
    ds = xr.open_dataset(exp.path_mask_mp_nc)
    da = xr.DataArray(arr_ts, name='timeseries', dims=['time', 'lat', 'lon'],
                        coords={'time': dates, 'lat': ds.lat, 'lon': ds.lon}
                        ).assign_attrs(units='m', ann_cycle_removed=int(rm_ann))

    comp = dict(zlib=True, complevel=3)
    try:
        encoding = {var: comp for var in da.data_vars}
    except:
        encoding = {da.name: comp}
    da.to_netcdf(dst, encoding=encoding)
    print ('Wrote timeseries to netcdf:', dst)

    return da


if __name__ == '__main__':
    path_opera_ts = Path(f'{os.getenv("dataroot")}/VLM/Sentinel1/Houston/OPERA/timeseries_cor_ERA5.h5')
    # path_opera_ts = Path(f'{os.getenv("dataroot")}/VLM/Sentinel1/NYC2/MintPy_2alks_5rlks_33_15/NYBK_ERA5_SET_PM_Stitch_ex_Fast/geo/geo_timeseries_SET_ERA5_demErr.h5')
    # path_opera_ts = Path(f'{os.getenv("dataroot")}/VLM/Sentinel1/NYC2/MintPy_2alks_5rlks_33_15/NYBK_ERA5_SET_PM_Stitch_ex_Fast/geo/geo_timeseries_SET_ERA5_demErr.h5')
    # da = h5_to_xr(path_opera_ts, 'timeseries')
    # print(da)
    from BZ import regions
    crop_utm(path_opera_ts, 'timeseries', regions['Houston'], dst=Path(f'{os.getenv("dataroot")}/VLM/Sentinel1/Houston'))
