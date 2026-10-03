#!/bin/zsh

# This program makes spectrum, RMF, and ARF files respectively.
# Please Make lists of raytracing event file before.
# Please prepare raytracing event files.



# Set some parameters.

minphoton=100
erange="0.3 16.0 2.0 8.0"
raytrace=($(<raytrace.list))



for reg_number in {?..?}   # reg_number is changeable.
do

    # Make spectrum files

[ -e xa"$OBSID"xtd_"$SUB"_src_"$reg_number".pi.gz ] && rm xa"$OBSID"xtd_"$SUB"_src_"$reg_number".pi.gz
xselect <<< 'xtd
read eve xa'"$OBSID"'xtd_'"$SUB"'_xpc_clnevt2.fits.gz .
filter region xa'"$OBSID"'xtd_'"$SUB"'_src_'"$reg_number"'.reg
extr spectrum
save spectrum xa'"$OBSID"'xtd_'"$SUB"'_src_'"$reg_number"'.pi
exit
no' 2>&1 | tee pi_"$reg_number".log
gzip xa"$OBSID"xtd_"$SUB"_src_"$reg_number".pi



    # Make RMF files

[ -e xa"$OBSID"xtd_"$SUB"_src_"$reg_number".rmf.gz ] && rm xa"$OBSID"xtd_"$SUB"_src_"$reg_number".rmf.gz
punlearn xtdrmf
xtdrmf infile=xa"$OBSID"xtd_"$SUB"_src_"$reg_number".pi.gz outfile=xa"$OBSID"xtd_"$SUB"_src_"$reg_number".rmf rmfparam=CALDB eminin=200 dein="2,24" nchanin="5900,500" eminout=0 deout=6 nchanout=4096 2>&1 | tee rmf_"$reg_number".log
gzip xa"$OBSID"xtd_"$SUB"_src_"$reg_number".rmf



    # Make ARF files

[ -e xa"$OBSID"xtd_"$SUB"_src_"$reg_number".arf.gz ] && rm xa"$OBSID"xtd_"$SUB"_src_"$reg_number".arf.gz
punlearn xaarfgen
xaarfgen xrtevtfile=raytrace_xa"$OBSID"xtd_"$SUB"_src_${raytrace[${reg_number}]}.fits.gz source_ra=1 source_dec=1 telescop=XRISM instrume=XTEND emapfile=xa"$OBSID"xtd_"$SUB".expo regmode=RADEC regionfile=xa"$OBSID"xtd_"$SUB"_src_"$reg_number".reg sourcetype=FLATCIRCLE flatradius=1 rmffile=xa"$OBSID"xtd_"$SUB"_src_"$reg_number".rmf.gz erange="$erange" outfile=xa"$OBSID"xtd_"$SUB"_src_"$reg_number".arf numphoton=1 minphoton="$minphoton" teldeffile=CALDB qefile=CALDB contamifile=CALDB obffile=CALDB fwfile=CALDB onaxisffile=CALDB onaxiscfile=CALDB mirrorfile=CALDB obstructfile=CALDB frontreffile=CALDB backreffile=CALDB pcolreffile=CALDB scatterfile=CALDB mode=h clobber=yes seed=7 imgfile=NONE

mv xaarfgen.log arf_"$reg_number".log
mv xaxmaarfgen.log xma_arf_"$reg_number".log

if ! [ -e xa"$OBSID"xtd_"$SUB"_src_"$reg_number".arf ]
then
    echo "\n The xa${OBSID}xtd_${SUB}_src_${reg_number}.arf cannot be created."
    echo " please change radius, numphoton, minphoton.\n"
    exit
fi

gzip xa"$OBSID"xtd_"$SUB"_src_"$reg_number".arf

done



# Remove unneceaary files

rm coordpnt.log
rm xselect.log
rm xtdrmf.log
rm xaarfgen_region.lis
rm arfgencgrid_full.fits
rm arfgencgrid.fits



# Move log files

[ -d log ] || mkdir log
mv -f *.log log