import h5py
import xarray as xr
import shutil
from mintpy.utils import readfile, writefile
import geopandas as gpd

from BZ import bbGIS

from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase

from BZ.bbLogger import logger

# for stitching time series
# needs to be updated to weight with either standard deviation of timeseries
# or rate uncertainty
def stitch_ts_NYC(Exp1, Exp2, excl_dates=[]):
    """ Project to East and West and then stitch together

    Writes a _stitch.h5 timeseries and _stitch netcdf
    """
    # get the crops

    epsg     = 'epsg:4326'

    lst_ts  = []
    lst_unc = []
    log.warning ('Need to implement mean in overlapping regions')
    for Expi in [Exp1, Exp2]:
        path_ts_nc  = Expi.path_ts_geo.replace('h5', 'nc')
        assert op.exists(path_ts_nc), \
            'Need to convert timeseries file to netcdf (Eof notebook)'

        ds_ts = xr.open_dataset(path_ts_nc)
        key   = 'timeseries_ann' if 'timeseries_ann' in ds_ts.keys() else 'timeseries'
        print ('Using ts key:', key)
        da_ts = ds_ts[key]
        da_ts.rio.write_crs(epsg, inplace=True)

        ## clip the rates to polygon
        path_crop    = op.join(Expi.path_crops, f'{CROP_MAP[Expi.ref_sta]}.GeoJSON')
        gdf_crop     = gpd.read_file(path_crop)
        da_ts_crop   = da_ts.rio.clip(gdf_crop.geometry, epsg, drop=False)

        lst_ts.append(da_ts_crop)

    ## get only overlapping dates
    set1 = set(lst_ts[0].time.data)
    set2 = set(lst_ts[1].time.data)
    common = list(set1.intersection(set2))
    common.sort()

    ## drop excluded manually excluded dates
    common = [c for c in common if not c in excl_dates]

    ix0 = lst_ts[0].time.isin(common)
    ix1 = lst_ts[1].time.isin(common)

    lst_ts[0] = lst_ts[0][ix0]
    lst_ts[1] = lst_ts[1][ix1]

    overlap = np.sum(lst_ts, 0)

    # nan out the overlap values in the western portion
    if not np.isnan(overlap).all():
        log.warning('There is overlap between the polygons, averaging the two.')

    # nan out the overlap values in the western portion
    # if overlap.any():
    #     overlap_nan = np.where(np.isnan(overlap), 0, np.nan)
    #     lst_ts[0] += overlap_nan

    ts_stitch = np.nanmean(lst_ts, 0)
    # ts_stitch = np.where(np.isnan(lst_ts[0].data),
    #                                  lst_ts[1].data, lst_ts[0].data)

    da_ts_stitch = da_ts.copy()
    da_ts_stitch = da_ts_stitch[da_ts_stitch.time.isin(common)]
    da_ts_stitch.data = ts_stitch

    ## just make both; they'll be the same
    for Expi in [Exp1, Exp2]:
        shutil.copy(Expi.path_ts_geo, Expi.path_ts_geo_stitch)
        with h5py.File(Expi.path_ts_geo_stitch, 'r+') as h5:
            del h5['timeseries']
            del h5['date']
            h5.create_dataset('timeseries', data=np.fliplr(ts_stitch))
            h5.create_dataset('date', data=np.array(common, dtype=np.string_))#[c.encod('utf-u') for c in common])##, data=mask, dtype=bool)

            logger.info('Stitched timeseries shape: %s', h5['timeseries'].shape)

        # writefile.write_hdf5_block(Expi.path_ts_geo_stitch, np.fliplr(ts_stitch), 'timeseries')
        # writefile.write_hdf5_block(Expi.path_ts_geo_stitch, [c.encod('utf-u') for c in common], 'date')
        logger.info ('Wrote: %s', Expi.path_ts_geo_stitch)

        ## update the netcdf
        da_ts_stitch.to_netcdf(Expi.path_ts_geo_nc_stitch)
        logger.info ('Wrote: %s', Expi.path_ts_geo_nc_stitch)

    return


## TO DO
def stitch_ts_mult(lst_exps, excl_dates=[]):
    """ Project timeseries to each station and then stitch

    Writes a _stitch.h5 timeseries and _stitch netcdf
    """
    epsg     = 'epsg:4326'

    lst_ts  = []
    lst_unc = []
    log.warning ('Need to implement mean in overlapping regions')
    for Expi in lst_exps:
        path_ts_nc  = Expi.path_ts_geo.replace('h5', 'nc')
        assert op.exists(path_ts_nc), \
            'Need to convert timeseries file to netcdf (Eof notebook)'

        ds_ts = xr.open_dataset(path_ts_nc)
        key   = 'timeseries_ann' if 'timeseries_ann' in ds_ts.keys() else 'timeseries'
        print ('Using ts key:', key)
        da_ts = ds_ts[key]
        da_ts.rio.write_crs(epsg, inplace=True)

        ## clip the rates to polygon
        path_crop    = op.join(Expi.path_crops, f'{CROP_MAP[Expi.ref_sta]}.GeoJSON')
        gdf_crop     = gpd.read_file(path_crop)
        da_ts_crop   = da_ts.rio.clip(gdf_crop.geometry, epsg, drop=False)

        lst_ts.append(da_ts_crop)

    ## get only overlapping dates
    set1 = set(lst_ts[0].time.data)
    set2 = set(lst_ts[1].time.data)
    common = list(set1.intersection(set2))
    common.sort()

    ## drop excluded manually excluded dates
    common = [c for c in common if not c in excl_dates]

    ix0 = lst_ts[0].time.isin(common)
    ix1 = lst_ts[1].time.isin(common)

    lst_ts[0] = lst_ts[0][ix0]
    lst_ts[1] = lst_ts[1][ix1]

    overlap = np.sum(lst_ts, 0)

    # nan out the overlap values in the western portion
    if not np.isnan(overlap).all():
        log.warning('There is overlap between the polygons, averaging the two.')

    # nan out the overlap values in the western portion
    # if overlap.any():
    #     overlap_nan = np.where(np.isnan(overlap), 0, np.nan)
    #     lst_ts[0] += overlap_nan

    ts_stitch = np.nanmean(lst_ts, 0)
    # ts_stitch = np.where(np.isnan(lst_ts[0].data),
    #                                  lst_ts[1].data, lst_ts[0].data)

    da_ts_stitch = da_ts.copy()
    da_ts_stitch = da_ts_stitch[da_ts_stitch.time.isin(common)]
    da_ts_stitch.data = ts_stitch

    ## just make both; they'll be the same
    for Expi in [Exp1, Exp2]:
        shutil.copy(Expi.path_ts_geo, Expi.path_ts_geo_stitch)
        with h5py.File(Expi.path_ts_geo_stitch, 'r+') as h5:
            del h5['timeseries']
            del h5['date']
            h5.create_dataset('timeseries', data=np.fliplr(ts_stitch))
            h5.create_dataset('date', data=np.array(common, dtype=np.string_))#[c.encod('utf-u') for c in common])##, data=mask, dtype=bool)

            logger.info('Stitched timeseries shape: %s', h5['timeseries'].shape)

        # writefile.write_hdf5_block(Expi.path_ts_geo_stitch, np.fliplr(ts_stitch), 'timeseries')
        # writefile.write_hdf5_block(Expi.path_ts_geo_stitch, [c.encod('utf-u') for c in common], 'date')
        logger.info ('Wrote: %s', Expi.path_ts_geo_stitch)

        ## update the netcdf
        da_ts_stitch.to_netcdf(Expi.path_ts_geo_nc_stitch)
        logger.info ('Wrote: %s', Expi.path_ts_geo_nc_stitch)

    return



if __name__ == '__main__':
    exp       = NYC_SRc
    mp_exp    = 'ERA5_SET_PM_Stitch_ex_Fast2017'

    ## Final Experiments
    ref_sta = 'NJHT'
    ExpW    = ExpBase(exp, mp_exp, ref_sta, neofs=20)

    ref_sta = 'NYBK'
    ExpE    = ExpBase(exp, mp_exp, ref_sta, neofs=20)

    stitch_ts_NYC(ExpW, ExpE)

