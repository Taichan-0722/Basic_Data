#!/bin/zsh

# This program can search the number of counts per selected region.

while :
do
    echo -n ' Enter the region number > '
    read reg_number

    if [ "$reg_number" = 'exit' ]
    then
        exit

    elif [ "$reg_number" = 'none' ]
    then
        echo
        xselect <<< 'xtd
read eve xa'"$OBSID"'xtd_'"$SUB"'_cl.evt.gz .
set image det
filter region exclude_calsources.reg
filter region exclude_anomaly.reg
filter region exclude_point.reg
filter time file ../xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt2.gti
select mkf "CORTIME.gt.6" mkf_name=xa'"$OBSID"'.ehk.gz mkf_dir=../analysis_xtd/
extr all
exit
no' | grep -A2 "Grand Total"
        echo
        rm xselect.log

    elif ! [ -e xa"$OBSID"xtd_"$SUB"_src_"$reg_number".reg ]
    then
        echo "  The region file doesn't exist."

    else
        echo
        xselect <<< 'xtd
read eve xa'"$OBSID"'xtd_'"$SUB"'_cl.evt.gz .
set image det
filter region xa'"$OBSID"'xtd_'"$SUB"'_src_'"$reg_number"'.reg
filter region exclude_calsources.reg
filter region exclude_anomaly.reg
filter region exclude_point.reg
filter time file ../xa'"$OBSID"'xtd_'"$SUB"'_asf_clnevt2.gti
select mkf "CORTIME.gt.6" mkf_name=xa'"$OBSID"'.ehk.gz mkf_dir=../analysis_xtd/
extr all
exit
no' | grep -A2 "Grand Total"
        echo
        rm xselect.log
    fi
done