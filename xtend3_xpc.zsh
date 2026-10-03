#!/bin/zsh

# Remove the unnecessary files.

cd analysis_xtd
[ -e xa"$OBSID"xtd_"$SUB"_xpc_hist.fits ] && rm xa"$OBSID"xtd_"$SUB"_xpc_hist.fits
[ -e xa"$OBSID"xtd_"$SUB"_xpc_reg01_hist.gif ] && rm xa"$OBSID"xtd_"$SUB"_xpc_reg01_hist.gif   # This sentence is changeable (reg02, reg03, ...)
[ -e xa"$OBSID"xtd_"$SUB"_xpc_clnevt.fits.gz ] && rm xa"$OBSID"xtd_"$SUB"_xpc_clnevt.fits.gz
[ -e xa"$OBSID"xtd_"$SUB"_xpc_rmvpix.fpix ] && rm xa"$OBSID"xtd_"$SUB"_xpc_rmvpix.fpix
[ -e xtdpixclip.log ] && rm xtdpixclip.log



# Execute the xtdpixclip with histo mode.
# Region files are changeable.

punlearn xtdpixclip
xtdpixclip evtinfile=xa"$OBSID"xtd_"$SUB"_cl.evt.gz outroot=xa"$OBSID"xtd_"$SUB" pmode=histo incregionfiles="xa${OBSID}xtd_${SUB}_xpc_inc01.reg" excregionfiles="xa${OBSID}xtd_${SUB}_xpc_exc01.reg" thresholds="50.0"

mv xtdpixclip.log xtdpixclip_histo.log
[ -e ../log_xtd/xtdpixclip_histo.log ] && rm ../log_xtd/xtdpixclip_histo.log
mv xtdpixclip_histo.log ../log_xtd

echo "\n histo mode has been completed."



# Execute the xtdpixclip with apply mode.
# Region files and the thresholds are changeable.

punlearn xtdpixclip
xtdpixclip evtinfile=xa"$OBSID"xtd_"$SUB"_cl.evt.gz outroot=xa"$OBSID"xtd_"$SUB" pmode=apply incregionfiles="xa${OBSID}xtd_${SUB}_xpc_inc01.reg" excregionfiles="xa${OBSID}xtd_${SUB}_xpc_exc01.reg" thresholds="50.0"
gzip xa"$OBSID"xtd_"$SUB"_xpc_clnevt.fits

mv xtdpixclip.log xtdpixclip_apply.log
[ -e ../log_xtd/xtdpixclip_apply.log ] && rm ../log_xtd/xtdpixclip_apply.log
mv xtdpixclip_apply.log ../log_xtd

echo "\n apply mode has been completed.\n"