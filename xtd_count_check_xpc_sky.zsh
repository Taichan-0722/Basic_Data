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
read eve xa'"$OBSID"'xtd_'"$SUB"'_xpc_clnevt2.fits.gz .
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
read eve xa'"$OBSID"'xtd_'"$SUB"'_xpc_clnevt2.fits.gz .
filter region xa'"$OBSID"'xtd_'"$SUB"'_src_'"$reg_number"'.reg
extr all
exit
no' | grep -A2 "Grand Total"
        echo
        rm xselect.log
    fi
done