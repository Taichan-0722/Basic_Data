#!/bin/zsh

# This program executes "rebin".
# The number of grouping minimum value is changeable.

for reg_number in {?..?}
do

[ -e xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".grp.gz ] && rm xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".grp.gz

grppha <<< 'xa'"$OBSID"'xtd_'"$SUB"'_src_'"$reg_number"'.pi.gz
xa'"$OBSID"'xtd_'"$SUB"'_beta_'"$reg_number"'.grp
chkey RESPFILE xa'"$OBSID"'xtd_'"$SUB"'_src_'"$reg_number"'.rmf.gz
chkey ANCRFILE xa'"$OBSID"'xtd_'"$SUB"'_beta_'"$reg_number"'.arf.gz
group min 1
exit' 2>&1 | tee grp_beta_"$reg_number".log

gzip xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".grp

done



# Move log files.

mv -f *.log log