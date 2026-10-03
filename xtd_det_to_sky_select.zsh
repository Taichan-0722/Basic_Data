#!/bin/zsh

# This program converts region files in DET coordinates into SKY coordinates.

echo -n "\n Enter the input file name > "
read infile

if ! [ -e ${infile} ]; then
    echo " '${infile}' doesn't exist.\n"
    exit 1
fi

echo -n " Enter the output file name > "
read outfile

punlearn coordpnt
coordpnt input="${infile}" outfile="${outfile}" telescop=XRISM instrume=XTEND teldeffile=CALDB startsys=DET stopsys=RADEC ra="$RA_NOM" dec="$DEC_NOM" roll="$PA_NOM" clobber=yes
echo " '${infile}' has been converted into '${outfile}'."