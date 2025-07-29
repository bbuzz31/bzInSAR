import geopandas as gpd
import xarray as xr

from BZ import bbGIS

from VLM.bzFRInGE import *
from VLM.bzFRInGE.experiments import *
from VLM.bzFRInGE.FRINGEBase import ExpBase


"""
InSAR at Bridges
"""
path_sent = f'{os.getenv("dataroot")}/VLM/Sentinel1'
ml_fr     = 'MintPy_2alks_5rlks_33_15'

dct_maps = dict(
    NYC    = f'{path_sent}/NYC2/{ml_fr}/NYBK_ERA5_SET_PM_Stitch_ex_Fast/geo/geo_rate_SR_recon20_masked_stitch.nc',
    NNJ    = f'{path_sent}NNJ/{ml_fr}/NJHT_ERA5_SET_PM_ex_Fast/geo/geo_rate_SR_recon15_masked.nc',
    Philly = f'{path_sent}/Philly/{ml_fr}/FTML_ERA5_PM_SET_ex_Fast/geo/geo_rate_SR_recon20_masked.nc',
    HR     = f'{path_sent}/HR/{ml_fr}/SPVA_ERA5_SET_PM_ex_Fast/geo/geo_rate_SR_recon25_masked.nc',
    Charleston = f'{path_sent}/{ml_fr}/MintPy_2alks_5rlks_33_15/SCHA_Base/geo/geo_rate_SR_recon15_masked.nc',
    Savannah   = f'{path_sent}/{ml_fr}/MintPy_2alks_5rlks_33_15/SAVA_ERA5_SET_PM_ex_Fast/geo/geo_rate_SR_masked.nc',
    FL     = f'{path_sent}/Miami/{ml_fr}/ZMA1_ERA5_PM_ex_Fast/geo/geo_rate_SR_masked.nc'
)



def main():
    src = Path(os.getenv('dataroot')) / 'GIS' / 'Bridges' / 'National_Bridge_Inventory_EC.GeoJSON'
    gdf = gpd.read_file(src)

    for k, v in dct_maps.items():
        ds = xr.open_dataset(v)
        da = ds['Band1'] * 1000 # convert to mm

        ## get the bridges in the region
        SNWE = da.lat.min().item(), da.lat.max().item(), da.lon.min().item(), da.lon.max().item()
        gdf_crop = bbGIS.in_dfll(gdf, SNWE=SNWE)

        ## iterate over each bridge to find its nearest location
        ## need to chagne this eventually to get closest within distance

        lst_dat = []
        for pt in gdf_crop.geometry:
            dat = da.sel(lat=pt.y, lon=pt.x, method='nearest').item()
            lst_dat.append(dat)

        gdf_crop[k] = lst_dat
        gdf_good = gdf_crop.dropna(subset=k)
        # gdf_good[gdf_good['deck_width_mt_052'] > 1]

        bbPlot.plot_dfll(gdf_good)
        breakpoint()






if __name__ == '__main__':
    main()