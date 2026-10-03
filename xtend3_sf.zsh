#!/bin/zsh

# This program removes anomalous pixels.

python -u search_flicker.py <<< 'analysis_xtd
xa'"$OBSID"'xtd_'"$SUB"'_bsf.img
xa'"$OBSID"'xtd_'"$SUB"'_sf.expo
xa'"$OBSID"'xtd_'"$SUB"'_asf.img
exclude_anomaly.reg
300
7
5' 2>&1 | tee search_flicker.log   # Threshold for simple counts cut, Poisson clean cell size, and Poisson probability threshold are changeable.

mv search_flicker.log log_xtd
mv analysis_xtd/exclude_anomaly.reg .