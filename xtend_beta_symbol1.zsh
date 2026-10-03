#!/bin/zsh

# This program copies some necessary files.

ln -s ../analysis_xtd/xa"$OBSID"xtd_"$SUB"_xpc_clnevt2_sky.img
ln -s ../analysis_xtd/xa"$OBSID"xtd_"$SUB".expo

cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend_beta_image.py .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend_beta_raytrace.zsh .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend_beta_arf.zsh .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend_panda_to_polygon.py .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend_beta_factor_area.py .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend_beta_grp.zsh .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend_det_to_sky.zsh .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend_sky_to_det.zsh .

touch photon.list
touch raytrace.list
touch radius.list