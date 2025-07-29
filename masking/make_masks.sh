cd masks
# OSM_wmask is created from OpenStreetMaps (see jupyter SE notebook(s))

# mintpy can make a radar coded from ./Mintpy/inputs
# generate_mask.py ./inputs/geometryRadar.h5 -m 37
# save_roipac.py mask.h5 # gdal readable

# ampMask made from 30 ifgs with my script on leffe
# makeAmpMask.py


################ SWBD
#### use ISCE to download and then mintpy geom
## use a big fucker
# wbd.py 36.0 38.0 -77.0 -75.0

## use mintpy to radarcode it 
# geocode.py ./swbdLat_N36_N38_Lon_W077_W075.wbd -o swbd.h5 -l ../MintPy_FR/inputs/geometryRadar.h5 --geo2radar

## then make the mask; -1 is water so mask below 0
# generate_mask.py ./swbd.h5 -o waterMask.h5 -m 0

## no mintpy tools work for this; use my own
# hdf2gdal.py waterMask.h5 --of envi

hdf2gdal.py waterMask_2alks_5rlks.h5 --of envi

## convert to byte for snaphu
gdal_translate waterMask_2alks_5rlks.envi  waterMask_2alks_5rlks_snaphu.msk -ot byte -of envi
