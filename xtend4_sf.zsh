#!/bin/zsh

# This program makes the new light curve and image file after removing the anomalous pixels.
# and execute CORTIME filtering.

xselect <<< 'xtd
read eve xa'"$OBSID"'xtd_'"$SUB"'_cl.evt.gz analysis_xtd
set image DET
filter region exclude_calsources.reg
filter region exclude_anomaly.reg
filter pha_cutoff 83 1667
select mkf "CORTIME.gt.6" mkf_name=xa'"$OBSID"'.ehk.gz mkf_dir=analysis_xtd/
extr all
save curve xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt.lc
save image xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt.img
exit
no' 2>&1 | tee xtend4.log

gzip xa"$OBSID"xtd_"$SUB"_asf_clnevt.lc
mv xtend4.log log_xtd
mv xa"$OBSID"xtd_"$SUB"_asf_clnevt.lc.gz analysis_xtd
mv xa"$OBSID"xtd_"$SUB"_asf_clnevt.img analysis_xtd
rm xsel_timefile.asc
rm xselect.log