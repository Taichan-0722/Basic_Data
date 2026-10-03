#!/bin/zsh

# This program makes some symbolic links of necessary files.

ln -s /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtd_count_check_sf_det.zsh
ln -s /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/exclude_calsources.reg
ln -s /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/exclude_edge.reg
ln -s /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/expand_panda.py
ln -s /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/split_panda.py
ln -s /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/split_panda_xtd.py

ln -s ../analysis_xtd/xa"$OBSID"xtd_"$SUB"_cl.evt.gz
ln -s ../analysis_xtd/xa"$OBSID"xtd_"$SUB"_asf_clnevt2_det.img
ln -s ../analysis_xtd/xa"$OBSID"xtd_"$SUB"_asf_clnevt2_sky.img
ln -s ../analysis_xtd/xa"$OBSID"xtd_"$SUB".expo

ln -s ../exclude_anomaly.reg
ln -s ../exclude_point.reg

cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend6_box_circle_sf_det.zsh .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend6_else_sf_det.zsh .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend6_reuse_sf_det.zsh .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtend7.zsh .
cp /Users/mitsuhashitaiki/Documents/analysis/Basic_Data/xtd_det_to_sky.zsh .

touch detx.list
touch dety.list
touch radius.list
touch photon.list
touch raytrace.list