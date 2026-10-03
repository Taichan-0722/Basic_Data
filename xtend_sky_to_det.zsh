#!/bin/zsh

# This program converts region files in DET coordinates into SKY coordinates.

input_file_name="xa${OBSID}xtd_${SUB}_sky_"
output_file_name="xa${OBSID}xtd_${SUB}_det_"

start_number=1
end_number=1

echo

for region_number in {"$start_number".."$end_number"}
do

    if ! [ -e ${input_file_name}${region_number}.reg ]; then
        echo " ${input_file_name}${region_number}.reg doesn't exist.\n"
        exit 1
    fi

    punlearn coordpnt
    coordpnt input="${input_file_name}${region_number}.reg" outfile="${output_file_name}${region_number}.reg" telescop=XRISM instrume=XTEND teldeffile=CALDB startsys=RADEC stopsys=DET ra="$RA_NOM" dec="$DEC_NOM" roll="$PA_NOM" clobber=yes
    echo " '${input_file_name}${region_number}.reg' has been converted into '${output_file_name}${region_number}.reg'."

done

echo