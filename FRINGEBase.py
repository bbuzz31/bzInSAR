from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *


class ExpBase(object):
    """ base class for setting up experiments """
    def __init__(self, dct_exp, mp_exp='Base', ref_sta=None, neofs='', v=True):
        self.dct_exp = dct_exp
        self.reg     = re.match('[A-Za-z]+', self.dct_exp['root']).group()
        # self.ref_sta = DCT_REG[self.reg][3][0] if ref_sta is None else ref_sta
        self.ref_sta = DCT_FRINGE[self.reg][3][0] if ref_sta is None else ref_sta
        self.mp_exp0 = mp_exp
        self.mp_exp  = f'{self.ref_sta}_{mp_exp}' # e.g., SCHA_Base
        self.SNWE    = DCT_FRINGE[self.reg][0]
        self.WESN    = [*self.SNWE[2:], *self.SNWE[:2]]
        self.SNWEs   = ' '.join([str(i) for i in self.SNWE])
        self.neofs   = neofs
        self.v       = v

        self._setup_track_frame()
        self._set_win_looks()
        self._set_network()
        self._set_lbl(neofs=neofs)
        self._set_unw_corrs()

        self.set_isce_paths()
        self.set_network_paths()
        self.set_fringe_paths()
        self.set_mintpy_paths()
        self.set_mintpy_paths_vup()

        self.set_netcdfs()
        self.make_struct()


        if not self.dct_exp['cust_net']:
            self.df_ifgs = self.symlink_wrapped()
            self.symlink_unwrapped()

        else:
            ## to do; support reading a reference file:
            pass

        return


    def _setup_track_frame(self, reg=None):
        reg        = self.reg if reg is None else reg
        self.track = str(DCT_REG[reg][1])
        self.frame = DCT_REG[reg][2]#[0]

        return


    def _set_win_looks(self, dct_exp=None):
        """ Set the half window size (x,y) and looks (naz, nrg) """
        dct_exp = self.dct_exp if dct_exp is None else dct_exp

        self.nx  = dct_exp['nx']
        self.ny  = dct_exp['ny']
        self.nxy = f'_{self.nx}_{self.ny}'

        self.naz = dct_exp['naz']
        self.nrg = dct_exp['nrg']

        ## set posting size; 10/15/30 m/90 m
        if self.naz == 1 and self.nrg == 1:
            self.post = '-0.00010416675 0.00010416675'

        elif self.naz == 3 and self.nrg == 3:
            self.post = '-0.000138889 0.0001388895'

        elif self.naz == 2 and self.nrg == 5:
            self.post = '-0.000277778 0.000277778'

        elif self.naz == 7 and self.nrg == 19:
            self.post = '-0.00083333 0.00083333'

        # guess? Sawsen
        elif self.naz == 9 and self.nrg == 4:
            self.post = '-0.00083333 0.00083333'

        else:
            print ('Unsupported multilook size for posting')
            self.post = None

        self.gap_fill = '.gf' if dct_exp['gap_fill'] else ''

        # for upsampling
        self.ml   = f'{self.naz}alks_{self.nrg}rlks' if self.naz>1 or self.nrg>1 else ''

        return


    def _set_network(self, dct_exp=None):
        dct_exp  = self.dct_exp if dct_exp is None else dct_exp
        network  = f'{dct_exp["sref"]}_' if dct_exp['sref'] else ''
        network += f'SBAS-{dct_exp["sbas"]}_' if dct_exp['sbas'] else ''
        if dct_exp['ann']:
            network += f'ANN{dct_exp["ann"]}' if isinstance(dct_exp['ann'],
                                            int) else f'ANN-{dct_exp["ann"]}'
        network     += dct_exp['cust_net']
        self.network = network.rstrip('_')
        return self.network


    def _set_lbl(self, dct_exp=None, neofs=''):
        dct_exp  = self.dct_exp if dct_exp is None else dct_exp
        self.lbl = self.network
        self.cust = dct_exp['custom']
        self.lbl = f'{self.lbl}_GACOS' if dct_exp['gacos'] else self.lbl
        self.lbl = f'{self.lbl}_ERA5' if dct_exp['era5'] else self.lbl
        self.lbl = f'{self.lbl}_brdg'  if dct_exp['bridging'] else self.lbl
        self.lbl = f'{self.lbl}_clos'  if dct_exp['closure'] else self.lbl
        self.lbl = f'{self.lbl}{self.cust}'

        self.lbl = f'{self.lbl}_recon{neofs}' if neofs else self.lbl
        self.lbl = f'{self.lbl}'
        self.neofs = neofs
        return self.lbl


    def _set_unw_corrs(self, dct_exp=None):
        dct_exp  = self.dct_exp if dct_exp is None else dct_exp
        bridging = f'_bridging' if dct_exp['bridging'] else ''
        closure  = f'_phaseClosure'  if dct_exp['closure'] else ''
        self.unw_corrs = f'{bridging}{closure}'
        return


    def set_isce_paths(self, dct_exp=None):
        dct_exp = self.dct_exp if dct_exp is None else dct_exp
        datar   = 'datarootd' if self.reg in ['Kennedy'] and not \
            socket.gethostname().startswith('MT-') else 'dataroot'

        self.path_vlm = Path(os.getenv(datar)) / 'VLM'

        self.path_wd  = self.path_vlm / dct_exp['satellite'] / self.dct_exp['root']
        self.path_stack_ifgs = self.path_wd / 'merged' / 'interferograms'

        ML = '_ML' if self.ml else ''
        self.path_isce_mp = self.path_wd / f'MintPy_ISCE{ML}'
        self.path_isce_mp.mkdir(exist_ok=True)
        return


    def set_fringe_paths(self, dct_exp=None):
        dct_exp = self.dct_exp if dct_exp is None else dct_exp
        self.path_crops = Path(os.getenv('dataroot')) / 'Crops' / self.reg

        self.path_corr = self.path_wraps_all / 'tcorr_ds_ps.bin'
        self.path_ps = self.path_wd / f'ampDispersion{self.nxy}' / 'ps_pixels'

        self.path_masks = self.path_wd / 'masks'
        self.path_mask = self.path_masks / 'waterMask.envi'

        self.path_geom = self.path_wd / 'geometry'
        self.path_dems = self.path_wd / 'DEM'

        ## need to be adjusted for 1 1 ; or possibly not used in full res?
        if self.ml:
            self.path_corr_ml = f'{self.path_corr.with_suffix("").as_posix()}_ISCE.{self.ml}.bin'
            # since using MP for this, use _
            self.path_mask_ml = f'{self.path_mask.with_suffix("").as_posix()}_ISCE_{self.ml}.envi'
            self.path_ps_ml = f'{self.path_ps.as_posix()}.{self.ml}'
            self.path_geom_ml = self.path_geom / f'geometry_{self.ml}'
        else:
            self.path_corr_ml = self.path_corr
            self.path_mask_ml = self.path_mask
            self.path_ps_ml = self.path_ps
            self.path_geom_ml = self.path_geom

        self.path_gps = self.path_vlm / 'GNSS' / 'UNR' / f'Midas_GPS-{self.reg}_cln.csv'
        self.path_unwraps_all = self.path_wraps_all / 'unwrap_ALL'

        self.path_shapes = self.path_wd / 'shapes' # GIS files
        return


    def set_network_paths(self, dct_exp=None):
        """ Set the network specific paths for prep_fringe to grab """
        dct_exp  = self.dct_exp.copy() if dct_exp is None else dct_exp

        ## paths to wrapped/unwrapped ifgs
        self.path_ps_ds     = self.path_wd / f'PS_DS{self.nxy}'
        self.path_wraps_all = self.path_ps_ds / 'ALL'
        self.path_wraps   = self.path_ps_ds / self.network
        self.path_unwraps = self.path_wraps / f'unwrap_{self.network}'

        return


    def set_mintpy_paths(self, sten=''):
        ## ------------------------------------------------------------- Base
        ml                = f'_{self.ml}' if self.ml else ''
        self.path_mp      = Path(self.path_wd) / f'MintPy{ml}{self.nxy}'
        self.path_mp_exp  = self.path_mp / self.mp_exp
        self.path_mp_exp_geo = self.path_mp_exp / 'geo'

        # unfortunately dont have unwrapping corrections for masks
        self.path_tc          = self.path_mp_exp / f'temporalCoherence_{self.network}.h5'
        self.path_tcf         = self.path_mp / f'temporalCoherence_{self.network}.h5'
        self.path_geom_mpf    = self.path_mp / 'inputs' / 'geometryRadar.h5'
        self.path_ifgramStackf = self.path_mp / 'inputs' / f'ifgramStack_{self.network}.h5'

        ## ------------------------------------------------------------- Ifgram/TS
        self.path_ts_raw      = self.path_mp / f'timeseries_{self.network}.h5'
        self.path_ts          = self.path_mp_exp / f'timeseries_{self.network}.h5'
        self.path_numInvIfgs  = self.path_mp / f'numInvIfgram_{self.lbl}.h5'
        self.path_avg_phs_vel = self.path_mp / f'avg_phs_vel_{self.lbl}.h5'
        self.path_ifgramStack = self.path_mp_exp / 'inputs' / f'ifgramStack_{self.network}.h5'
        self.path_geom_mp     = self.path_mp_exp / 'inputs' / 'geometryRadar.h5'


        ## Masks
        ## ------------------------------------------------------------- Mask
        self.path_mask_mpf     = self.path_mp / 'waterMask.h5'
        self.path_mask_mp      = self.path_mp_exp / 'waterMask.h5'
        self.path_mask_mp_geof = self.path_mp / 'geo_waterMask.h5'
        self.path_mask_mp_geo  = self.path_mp_exp_geo / 'geo_waterMask.h5'
        self.path_mask_mp_nc   = self.path_mask_mp_geo.with_suffix('.nc')

        # temporal coherence
        self.path_mask_tc     = self.path_mp_exp / f'maskTempCoh.h5'
        self.path_mask_tcf    = self.path_mp / f'maskTempCoh.h5'

        # conncomp
        self.path_mask_ccf    = self.path_mp / f'maskConnComp_{self.lbl}.h5'
        self.path_mask_cc     = self.path_mp_exp / f'maskConnComp_{self.lbl}.h5'
        # self.path_mask_cc_geof = self.path_mp / f'geo_maskConnComp_{self.lbl}.h5'
        self.path_mask_cc_geo  = self.path_mp_exp_geo / f'geo_maskConnComp_{self.lbl}.h5'

        # map of timeseries standard deviation (may not actually acount for sten)
        self.path_tsig_geo    = self.path_mp_exp_geo / f'geo_TSig_{self.lbl}{sten}.h5'

        # map of residue (sum of |data - model|)
        self.path_mask_residue_geo = self.path_mp_exp_geo / f'geo_residue_{self.lbl}{sten}.h5'

        # mask where the gradient in time and space is above a threshold
        self.path_mask_grad = self.path_mp_exp_geo / f'gradientMask_{self.lbl}{sten}.h5'

        return


    def set_mintpy_paths_vup(self, sten=''):
        """ Set the geo velocity and experiment Vup files for template """
        mp_lbl = '_'.join(self.mp_exp.split('_')[1:])
        if 'SET' in mp_lbl and 'ERA5' in mp_lbl:
            lbl1 = '_SET_ERA5'
        else:
            lbl1 = '' if 'Base' in mp_lbl else f'_{mp_lbl.split("_")[0]}' # for SET/ERA5
        lbl1   = lbl1.replace('_Bridging', '') # i just update the timeseries directly
        lbl1   = lbl1.replace('_PM', '') # i just update the timeseries directly
        lbl2   = self.lbl.replace('SR', '') # for getting eofs

        self.path_vlos_geo    = self.path_mp_exp_geo / f'geo_velocity{sten}.h5'
        self.path_vlos_geo_pm = self.path_vlos_geo.parent / f'{self.path_vlos_geo.stem}_ITRF.h5'

        self.path_geom_mp_geo = self.path_mp_exp_geo / 'geo_geometryRadar.h5'
        self.path_tc_geo      = self.path_mp_exp_geo / 'geo_temporalCoherence.h5'

        demErr = '_demErr' if not 'nodem' in self.mp_exp else ''
        self.path_ts_geo  = self.path_mp_exp_geo / f'geo_timeseries{lbl1}{demErr}{lbl2}.h5'
        self.path_ts_geo_stitch  = self.path_ts_geo.parent / f'{self.path_ts_geo.stem}_stitch.h5'

        ## timeseries with annual signal removed
        self.path_ts_geo_ann = self.path_ts_geo.parent / f"{self.path_ts_geo.stem.replace(f'_recon{self.neofs}', '')}_ann.h5"

        ## Vup
        self.path_vup_geo = self.path_mp_exp_geo / f'geo_Vup_{mp_lbl}_{self.lbl}{sten}.h5'
        self.path_vup_geo_na12 = self.path_vup_geo.parent / f'{self.path_vup_geo.stem}_NA12.h5'
        self.path_vup_geo_stitch = self.path_vup_geo.parent / f'{self.path_vup_geo.stem}_stitch.h5'

        ## Masks
        ## PS points
        self.path_mask_ps_geo = self.path_mp_exp_geo / f'geo_maskPS.h5'

        # temporal coherence
        self.path_mask_tc_geo = self.path_mp_exp_geo / f'geo_maskTempCoh.h5'

        # mask where velocity > (Nv * sigma)
        self.path_mask_velS_geo = self.path_mp_exp_geo / f'geo_maskVelS_{mp_lbl}{sten}.h5'
        self.path_mask_velS_geo_stitch = self.path_mask_velS_geo.parent / f'{self.path_mask_velS_geo.stem}_stitch.h5'

        # mask where velocity values are too high or low
        self.path_mask_vel_geo = self.path_mp_exp_geo / f'geo_maskVel_{mp_lbl}{sten}.h5'
        self.path_mask_vel_geo_stitch = self.path_mask_vel_geo.parent / f'{self.path_mask_vel_geo.stem}_stitch.h5'

        # mask where rate sig values are too high
        self.path_mask_vsig_geo = self.path_mp_exp_geo / f'geo_maskVSig_{mp_lbl}{sten}.h5'
        self.path_mask_vsig_geo_stitch = self.path_mask_vsig_geo.parent / f'{self.path_mask_vsig_geo.stem}_stitch.h5'

        # map of timeseries standard deviation
        self.path_tsig_geo = self.path_mp_exp_geo / f'geo_Tsig.h5'

        # mask where TS sig values are too high
        self.path_mask_tsig_geo = self.path_mp_exp_geo / f'geo_maskTsig.h5'

        # mask where the gradient in time and space is above a threshold
        self.path_mask_grad     = self.path_mp_exp / f'gradientMask_{mp_lbl}.h5'
        self.path_mask_grad_geo = self.path_mp_exp_geo / f'geo_gradientMask_{mp_lbl}.h5'

        self.path_cmp = Path(self.path_wd) / 'comparisons'
        os.makedirs(self.path_cmp, exist_ok=True)
        return

    def set_netcdfs(self, dct_exp=None, sten=''):
        """ Set the paths for the netcdf files for GMT """
        dct_exp  = dct_exp if dct_exp else self.dct_exp

        ## combined mask
        self.path_mask_vup = self.path_mp_exp_geo / f'geo_mask_Vup_{self.lbl}{sten}.h5'
        self.path_mask_vup_stitch = self.path_mask_vup.parent / f'{self.path_mask_vup.stem}_stitch.h5'

        ## water/conncomp/combined masks
        self.path_mask_mp_nc = self.path_mask_mp_geo.with_suffix('.nc')
        self.path_mask_tc_nc = self.path_mask_tc_geo.with_suffix('.nc')
        self.path_mask_vup_nc = self.path_mask_vup.with_suffix('.nc')
        self.path_mask_vup_nc_stitch = self.path_mask_vup_stitch.with_suffix('.nc')

        ## rate/unc, unmasked/masked
        self.path_rate_nc = self.path_mp_exp_geo / f'geo_rate_{self.lbl}{sten}.nc'
        self.path_rate_nc_stitch = self.path_rate_nc.parent / f'{self.path_rate_nc.stem}_stitch.nc'
        self.path_std_nc = self.path_mp_exp_geo / f'geo_std_{self.lbl}{sten}.nc'
        self.path_std_nc_stitch = self.path_std_nc.parent / f'{self.path_std_nc.stem}_stitch.nc'
        self.path_resid_nc = self.path_mp_exp_geo / f'geo_resid_{self.lbl}{sten}.nc'
        self.path_resid_nc_stitch = self.path_resid_nc.parent / f'{self.path_resid_nc.stem}_stitch.nc'

        self.path_rate_msk_nc = self.path_rate_nc.parent / f'{self.path_rate_nc.stem}_masked.nc'
        self.path_rate_msk_nc_stitch = self.path_rate_msk_nc.parent / f'{self.path_rate_msk_nc.stem}_stitch.nc'
        self.path_std_msk_nc = self.path_std_nc.parent / f'{self.path_std_nc.stem}_masked.nc'
        self.path_std_msk_nc_stitch = self.path_std_msk_nc.parent / f'{self.path_std_msk_nc.stem}_stitch.nc'
        self.path_resid_msk_nc = self.path_resid_nc.parent / f'{self.path_resid_nc.stem}_masked.nc'
        self.path_resid_msk_nc_stitch = self.path_resid_msk_nc.parent / f'{self.path_resid_msk_nc.stem}_stitch.nc'

        self.path_vel_kmz = self.path_vup_geo.with_suffix('') / '_rate'
        self.path_vel_kmz_stitch = self.path_rate_nc_stitch.with_suffix('')
        self.path_std_kmz = self.path_std_nc.with_suffix('')
        self.path_std_kmz_stitch = self.path_std_nc_stitch.with_suffix('')
        self.path_resid_kmz = self.path_resid_nc.with_suffix('')
        self.path_resid_kmz_stitch = self.path_resid_nc_stitch.with_suffix('')

        self.path_ts_geo_nc = self.path_ts_geo.with_suffix('.nc')
        self.path_ts_geo_nc_stitch = self.path_ts_geo_stitch.with_suffix('.nc')
        return

    def set_sten_paths(self, sten=''):
        """ For analyzing a temporal subset """
        self.set_mintpy_paths(sten=sten)
        self.set_mintpy_paths_vup(sten=sten)
        self.set_netcdfs(sten=sten)
        return

    def mv_ifgramStack(self):
        """ Move the ifgramStack, or just get the path """
        f_orig = 'ifgramStack'
        f_new  = f'ifgramStack_{self.network}'
        src    = self.path_mp / 'inputs' / f'{f_orig}.h5'
        dst    = self.path_mp / 'inputs' / f'{f_new}.h5'
        os.rename(src, dst) if not dst.exists() else ''
        return dst


    def add_wrapped(self, path_stack=None, path_wraps=None, bbox=None):
        """ Add the wrapped files to an existing ifgramStack

        Need to crop with a bbox (SNWE) in radarcoords if using crop with prep_fringe
        """
        import h5py
        path_stack = self.path_ifgramStackf if path_stack is None else path_stack
        path_wraps = self.path_wraps if path_wraps is None else path_wraps

        path_wraps = list(sorted(path_wraps.glob('*.int')))

        if not path_wraps:
            print ('Couldnt find any matching wrapped ifgs')
            return


        with h5py.File(path_stack, 'a') as h5:
            if 'wrapPhase' in h5.keys():
                # print ('wrapPhase already exists! Not overwriting.')
                del h5['wrapPhase']
                # return

            shp  = h5['unwrapPhase'][:].shape
            dum  = np.zeros(shp)
            dset = h5.create_dataset('wrapPhase', shp)
            print ('Output shape:', shp)

            for i, wrap_file in enumerate(path_wraps):
                ds   = gdal.Open(wrap_file, gdal.GA_ReadOnly)
                data = np.angle(ds.ReadAsArray())

                # note coords might need to be flipped
                data = data[bbox[0]:bbox[1], bbox[2]:bbox[3]] if bbox else data

                dset[i] = np.array(data, dtype=np.float32)
        print ('Wrote wrapPhase to:', path_stack)
        return


    def get_rate_unc_h5(self, mask=True, path_vup_geo=None, units='m'):
        """ Get the rate, uncertainty and metadata from mintpy files """
        from mintpy.utils import readfile
        dct_units = {'m': 1, 'dm': 10, 'cm': 100, 'mm': 1000}
        path_vup_geo = self.path_vup_geo if path_vup_geo is None else self.path_vup_geo
        arr_rate, meta = readfile.read(path_vup_geo)
        arr_unc = readfile.read(path_vup_geo, datasetName='velocityStd')[0]
        arr_mask = readfile.read(self.path_mask_vup)[0] if mask else 1.0
        arr_mask = np.where(arr_mask, 1, np.nan)

        uf = dct_units[units]
        arr_rate = arr_rate * uf * arr_mask
        arr_unc = arr_unc * uf * arr_mask

        return arr_rate, arr_unc, arr_mask, meta


    def get_rate_unc_nc(self, mask=True, units='mm'):
        """ Get the rate/velocity maps in mm  """
        import xarray as xr
        dct_units = {'m': 1, 'dm': 10, 'cm': 100, 'mm': 1000}
        if mask:
            rate = xr.open_dataset(self.path_rate_msk_nc)
            unc = xr.open_dataset(self.path_std_msk_nc)
        else:
            rate = xr.open_dataset(self.path_rate_nc)
            unc = xr.open_dataset(self.path_std_nc)

        # hack for marin
        # log.critical('Using marins data')

        # S, N, W, E = self.SNWE
        # path_disp = self.path_wd / 'OPERA' / 'Houston_3d_disp.h5'
        # da_utm = h5_to_xr(path_disp, 'vlm').rio.write_crs(32615)
        # da_vup = da_utm.rio.reproject(4326).rename(x='lon', y='lat')
        # da_vup = da_vup.sel(lat=slice(N, S), lon=slice(W, E))
        # da_mask = rate.interp_like(da_vup) # tcoh and water from fringe
        # da_mask = xr.where(da_mask.notnull(), 1, np.nan)

        # rate = da_vup.where(np.abs(da_vup)>0)/1000 * da_mask

        # da_utm = h5_to_xr(path_disp, 'vlm_std').rio.write_crs(32615)
        # da_std = da_utm.rio.reproject(4326).rename(x='lon', y='lat')
        # da_std = da_std.sel(lat=slice(N, S), lon=slice(W, E))
        # da_std = da_std.where(np.abs(da_std)>0)/1000
        # unc = da_std * da_mask

        uf = dct_units[units]
        return (rate['Band1']*uf).rio.write_crs(4326), (unc['Band1']*uf).rio.write_crs(4326)


    def make_vel_path(self, sta):
        """ Make the LOS velocities from the timeseries  """
        path_los = self.path_mp / f'geo_vel_LOS_{sta.upper()}_{self.lbl}.h5'
        return path_los


    def _ref_dict(self):
        """ x/y pixel coordinates for each station, for each multilook """
        dct = {'2alks_5rlks': {'LOY2': [439, 7953],
                               'VAHP': [2961, 5375],
                               'LOYZ': [1584, 69]
                               }
               }
        return dct


    def make_struct(self):
        """ Make directories """
        lst_dirs = [self.path_wd, self.path_dems, self.path_masks, self.path_crops,
                    self.path_mp, self.path_unwraps_all, self.path_unwraps]
        [os.makedirs(ld, exist_ok=True) for ld in lst_dirs]
        return


    ## not used?
    def get_wifgs(self, dct_exp=None, hdrs=True):
        """ Get the paths to the wrapped ifgs (.int), optionally also hdrs """
        dct_exp    = self.dct_exp if dct_exp is None else dct_exp
        path_wraps = get_ifgs(self.path_wraps, dct_exp, ext='int')

        if hdrs:
            path_wraps += get_ifgs(self.path_wraps, dct_exp, ext='hdr')

        return path_wraps


    def get_ifg_dates(self, kind='sr'):
        """ Get the ifgs for a specific network from the coreg files """
        path_coreg = self.path_wd / 'coreg_stack'
        if kind.lower() == 'sr':
            vrt_name = 'slcs_base.vrt'

        elif kind.lower() == 'nn':
            vrt_name = 'slcs_nn*.vrt'

        elif kind.lower() == 'ann':
            vrt_name = 'slcs_ann*.vrt'

        path_vrts = sorted(list(path_coreg.glob(vrt_name)))

        ifgs      = []
        ## put this in for NYC which is missing the actual slcs
        if len(path_vrts) == 1 and not self.reg == 'NYC':
            for vrt in path_vrts:
                # print (vrt)
                ds    = gdal.Open(vrt, gdal.GA_ReadOnly)
                dates = []
                files = ds.GetFileList()
                for f in files:
                    dt = op.splitext(op.basename(f))[0]
                    if dt.startswith('2'):
                        dates.append(dt)
                # now make the ifgs
                ref   = dates.pop(0)
                ifgs.extend([f'{ref}_{d}' for d in dates])

        ## in case the actual slcs are missing
        if not ifgs and self.reg == 'NYC' and socket.gethostname().startswith('leffe') \
            and op.exists(path_vrt):

            dates   = []
            with open(path_vrt, 'r') as fh:
                txt = fh.read()
            # iterate over the source filenames in the VRT
            for t in txt.split():
                if t.startswith('<SourceFilename'):
                    dt = op.splitext(op.basename(t[16:-17]))[0]
                    if dt.startswith('2'):
                        dates.append(dt)
            ref   = dates.pop(0)
            ifgs.extend([f'{ref}_{d}' for d in dates])

        return ifgs


    def symlink_wrapped(self, dct_exp=None):
        """ Symlink a customish network from 'ALL' to the correct dir.

        Leverages the slcs_stack folder, multilooking, and gapfilling
        """
        dct_exp  = self.dct_exp if dct_exp is None else dct_exp

        ifg_dates = []

        if dct_exp['sref']:
            ifg_dates.extend(self.get_ifg_dates('sr'))

        if dct_exp['sbas']:
            ifg_dates.extend(self.get_ifg_dates('nn'))

        if dct_exp['ann']:
            ifg_dates.extend(self.get_ifg_dates('ann'))

        naz, nrg = dct_exp['naz'], dct_exp['nrg']
        ext      = f'gf.*' if dct_exp['gap_fill'] else '.*'

        ## get the wrapped
        search = f'20*{naz}alks_{nrg}rlks.{ext}' if (naz > 1 or nrg > 1) else f'20*[0-9].{ext}'

        path_wraps_all = self.path_wraps_all.glob(search)
        path_wraps = [pw for pw in path_wraps_all if pw.name[:17] in ifg_dates]
        path_wraps = sorted(list(set(path_wraps)))
        ## unlink all
        [f.unlink() for f in self.path_wraps.glob('*') if f.is_symlink()]


        for pw in path_wraps:
            (self.path_wraps / pw.name).symlink_to(pw)

        df = pd.DataFrame(path_wraps, columns=['path'])
        df['ifg'] = df['path'].apply(lambda x: op.basename(x)[:17])
        df['ref'] = df['ifg'].apply(lambda x: x.split('_')[0])
        df['sec'] = df['ifg'].apply(lambda x: x.split('_')[1])
        df.drop_duplicates('ifg', inplace=True)
        if df.empty:
            print (f'No {self.network} wrapped ifgs at: {search}.')

        else:
            print (f'Got {df.shape[0]} wrapped ifgs ({self.network})')

        return df['ref sec ifg path'.split()]


    def symlink_unwrapped(self):
        """ Symlink unwrapped ifgs (if they exist) from unwrap_ALL to unwrap_NET

        Gets only those for the custom network by using those in the wrapped dir
        """
        [f.unlink() for f in self.path_unwraps.glob('*') if f.is_symlink()]

        i = 0
        for f in os.listdir(self.path_wraps):
            search    = f.replace('int', 'unw*')
            path_unws = sorted(list(self.path_unwraps_all.glob(search)))
            if not path_unws:
                continue
            for puw in path_unws:
                (self.path_unwraps / puw.name).symlink_to(puw)
            i+=1

        # print (f'Got {i} unwrapped ifgs')
        return


    def get_ref_stas(self):
        """ Mainly for getting stitched reference stations from vel files """
        from mintpy.utils import readfile
        ### grab the reference directly from the velocity files
        use_stitched = True if 'stitch' in self.mp_exp0.lower() else False
        path_vel = self.path_vup_geo_stitch if use_stitched else self.path_vup_geo
        ref_stas = readfile.read_attribute(path_vel).get('ref_stas')
        ref_stas = [ref_stas] if isinstance(ref_stas, str) else ref_stas

        ### below is deprecated
        # stas      = df_gps.sta.tolist()
        # stas_good = [sta in self.mp_exp0 for sta in stas]
        # try:
        #     sta2 = np.array(stas)[stas_good].item()
        # except:
        #     sta2 = False

        # use_stitched = True if 'stitch' in self.mp_exp0.lower() else False
        # if use_stitched:
        #     log.warning(f'Using stitched result (stitch sta: {stitch_sta})')
        #     ref_stas = [self.ref_sta, stitch_sta]

        # elif sta2:
        #     # fancy weighted stitching for NYC
        #     ref_stas = [self.ref_sta, sta2]

        # else:
        #     #  single reference
        #     ref_stas = [self.ref_sta]
        return ref_stas


    def get_vel_unc_near(self, lalo, npix=0, mask=True, lbl=None):
        from mintpy import utils
        arr_rate, arr_unc, arr_mask, meta = self.get_rate_unc_h5(mask)
        arr_rate *= 1000
        arr_unc *= 1000
        coord = utils.utils.coordinate(meta, lookup_file=self.path_geom_mp_geo)
        y, x = coord.geo2radar(lalo[0], lalo[1])[0:2]
        if npix:
            arr_rate1 = np.nanmean(arr_rate[y-npix:y+npix, x-npix:x+npix])
            arr_unc1 = np.nanmean(arr_unc[y-npix:y+npix, x-npix:x+npix])
            if np.isnan(arr_rate1):
                log.warning(f'y={y} x={x} is all nan despite averaging!')
        else:
            if isinstance(arr_mask, np.ndarray):
                if not arr_mask[y,x]:
                    log.warning(f'y={y} x={x} is masked!')
            arr_rate1 = arr_rate[y,x]
            arr_unc1 = arr_unc[y,x]

        lbl = f'y={y} x={x}' if lbl is None else f'{lbl}'
        log.info(f'{lbl}: {arr_rate1:.2f} +/- {arr_unc1:.2f} (mm/yr)')

        return arr_rate1, arr_unc1


class DolphinBase(ExpBase):
    def __init__(self, dct_exp, mp_exp='Base', ref_sta=None, neofs='', v=True):
        super().__init__(dct_exp, mp_exp, ref_sta, neofs, v)
        self.set_dolphin_paths()

    def set_dolphin_paths(self):
        # self.path_ps      = self.path_wd / f'ampDispersion{self.nxy}' / 'ps_pixels'
        # self.path_ps_ds     = self.path_wd / f'PS_DS{self.nxy}'
        self.path_wraps_all     = self.path_wd / 'interferograms'
        self.path_wraps         = self.path_wraps_all
        self.path_unwraps_all   = self.path_wd / 'unwrap_ALL'
        self.path_unwraps       = self.path_unwraps_all / f'unwrap_{self.network}'
        self.path_corr          = self.path_wraps / 'temporal_coherence.tif'
        self.symlink_unwrapped()

        self.path_geom_mp     = self.path_mp / 'geometryGeo.h5'
        self.path_vlos_geo    = self.path_mp_exp / 'velocity.h5'

# for testing
if __name__ == '__main__':
    # Obj = ExpBase(HR_90c).add_wrapped(path_wraps='/u/leffe-data2/buzzanga/data/VLM/Sentinel1/HR2/PS_DS_33_15/Custom')
    # prep_GPS_HR(Obj)
    from BZ.bbPlot import plot_bbox
    Obj = ExpBase(NYC_SR, 'ERA5_Bridging', 'NJHT')
    print (Obj.SNWEs)
    # plot_bbox(Obj.SNWE)
    # plt.show()

    # edit_runfiles(Obj.path_wd)
