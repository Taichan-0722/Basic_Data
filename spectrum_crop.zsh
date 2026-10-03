#!/bin/zsh

echo

for qdp_file in *.qdp; do

base_name=${qdp_file%.*}
qdp <<< ''"$qdp_file"'
/xw
la t
time off
plot
ha '"$base_name"'.ps/cps
q' &> /dev/null
ps2pdf "${base_name}.ps" "${base_name}.pdf"
echo " ${base_name} has been done"

done

echo