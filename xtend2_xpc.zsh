#!/bin/zsh

# This program makes the light curve and image file to select the region.

xselect <<< 'xtd
read eve xa'"$OBSID"'xtd_'"$SUB"'_cl.evt.gz analysis_xtd
set image DET
set binsize 128.0
extr all
save curve xa'"$OBSID"'xtd_'"$SUB"'_xpc.lc
save image xa'"$OBSID"'xtd_'"$SUB"'_xpc.img
exit
no' 2>&1 | tee xtend2.log

gzip xa"$OBSID"xtd_"$SUB"_xpc.lc
mv xtend2.log log_xtd
mv xa"$OBSID"xtd_"$SUB"_xpc.lc.gz analysis_xtd
mv xa"$OBSID"xtd_"$SUB"_xpc.img analysis_xtd
rm xselect.log