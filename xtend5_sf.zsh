#!/bin/zsh

# This program extracts the new event file, light curve and image file filtering with GTI file and point source region file.
# and execute CORTIME filtering.

xselect <<< 'xtd
read eve xa'"$OBSID"'xtd_'"$SUB"'_cl.evt.gz analysis_xtd
set image DET
filter region exclude_calsources.reg
filter region exclude_anomaly.reg
filter region exclude_point.reg
filter pha_cutoff 83 1667
filter time file xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt2.gti
select mkf "CORTIME.gt.6" mkf_name=xa'"$OBSID"'.ehk.gz mkf_dir=analysis_xtd/
extr all
save curve xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt2.lc
save image xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt2_det.img
clear pha_cutoff
extr eve
save eve xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt2.fits
yes
exit
no' 2>&1 | tee xtend5.log

xselect <<< 'xtd
read eve xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt2.fits .
filter pha_cutoff 83 1667
extr image
save image xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt2_sky.img
exit
no' 2>&1 | tee xtend5_sky.log

gzip xa"$OBSID"xtd_"$SUB"_asf_clnevt2.fits
gzip xa"$OBSID"xtd_"$SUB"_asf_clnevt2.lc
mv xtend5.log log_xtd
mv xtend5_sky.log log_xtd
mv xa"$OBSID"xtd_"$SUB"_asf_clnevt2.fits.gz analysis_xtd
mv xa"$OBSID"xtd_"$SUB"_asf_clnevt2.lc.gz analysis_xtd
mv xa"$OBSID"xtd_"$SUB"_asf_clnevt2_det.img analysis_xtd
mv xa"$OBSID"xtd_"$SUB"_asf_clnevt2_sky.img analysis_xtd
rm xsel_timefile.asc
rm xselect.log