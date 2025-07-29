############################# Set up the multilooked geometry files for MintPy
alks=7
rlks=19
cd /Users/buzzanga/data/VLM/Sentinel1/HR_Exps
cd geometry_multilook2

## must have made ../geometry/los.vrt shadowmask and incLocal
# (copy ../geometry/lat.vrt and change path and copy band)
cp ../geometry/*.vrt . &&
#
mv hgt.vrt hgt_full && # otherwise vrt extension gets propagated
multilookWrapped.py hgt_full --naz ${alks} --nrg ${rlks} &&
mv hgt_full hgt_full.vrt && # since its actually a vrt
#
mv lat.vrt lat_full && # otherwise vrt extension gets propagated
multilookWrapped.py lat_full --naz ${alks} --nrg ${rlks} &&
mv lat_full lat_full.vrt && # since its actually a vrt

mv lon.vrt lon_full && # otherwise vrt extension gets propagated
multilookWrapped.py lon_full --naz ${alks} --nrg ${rlks} &&
mv lon_full lon_full.vrt && # since its actually a vrt

mv los.vrt los_full  &&
multilookWrapped.py los_full --naz ${alks} --nrg ${rlks} &&
mv los_full los_full.vrt && # since its actually a vrt

mv shadowMask.vrt shadowMask_full  &&
multilookWrapped.py shadowMask_full --naz ${alks} --nrg ${rlks} &&
mv shadowMask_full shadowMask_full.vrt && # since its actually a vrt

mv incLocal.vrt incLocal_full  &&
multilookWrapped.py incLocal_full --naz ${alks} --nrg ${rlks} &&
mv incLocal_full incLocal_full.vrt && # since its actually a vrt

############################################ make the full files that get read
ln -s hgt_full_ISCE.${alks}alks_${rlks}rlks hgt.rdr.full &&
ln -s lon_full_ISCE.${alks}alks_${rlks}rlks lon.rdr.full &&
ln -s lat_full_ISCE.${alks}alks_${rlks}rlks lat.rdr.full &&
ln -s los_full_ISCE.${alks}alks_${rlks}rlks los.rdr.full &&
ln -s shadowMask_full_ISCE.${alks}alks_${rlks}rlks shadowMask.rdr.full &&
ln -s incLocal_full_ISCE.${alks}alks_${rlks}rlks incLocal.rdr.full &&

ln -s hgt_full_ISCE.${alks}alks_${rlks}rlks.vrt hgt.rdr.full.vrt &&
ln -s lon_full_ISCE.${alks}alks_${rlks}rlks.vrt lon.rdr.full.vrt &&
ln -s lat_full_ISCE.${alks}alks_${rlks}rlks.vrt lat.rdr.full.vrt &&
ln -s los_full_ISCE.${alks}alks_${rlks}rlks.vrt los.rdr.full.vrt &&
ln -s shadowMask_full_ISCE.${alks}alks_${rlks}rlks.vrt shadowMask.rdr.full.vrt &&
ln -s incLocal_full_ISCE.${alks}alks_${rlks}rlks.vrt incLocal.rdr.full.vrt &&

################################################ multilookWrapped.py ps_pixels
# multilook the ps_pixels
cd "$(dirname "$0")" # change to location of script
cd ../ampDispersion


multilookWrapped.py ps_pixels --naz ${alks} --nrg ${rlks} &&

mv ps_pixels ps_pixels_FR &&
mv ps_pixels.xml ps_pixels_FR.xml &&
mv ps_pixels.vrt ps_pixels_FR.vrt &&

ln -s ps_pixels_ISCE.${alks}alks_${rlks}rlks  ps_pixels &&
ln -s ps_pixels_ISCE.${alks}alks_${rlks}rlks.vrt  ps_pixels.vrt &&
ln -s ps_pixels_ISCE.${alks}alks_${rlks}rlks.xml  ps_pixels.xml &&


#################################### make the original lat.vrt which gets read
cd "$(dirname "$0")" # change to location of script
cd geometry_multilook
cp lat.rdr.full.vrt lat.vrt &&
## then you have to write the multilooked src rect physically into the lat.vrt
echo "Manually edit the src rect in lat.vrt"
