#!/bin/zsh

# This program makes beta ARF respectively.



# Set some parameters.

photon=10000000
minphoton=100
erange="0.3 16.0 0.5 10.0"
ra=192.2035533
dec=-41.3124625
image_file="aaa.img"



for reg_number in {?..?}   # reg_number is changeable.
do

    # Make a raytracing event file.

[ -e raytrace_xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".fits ] && rm raytrace_xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".fits
[ -e raytrace_xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".fits.gz ] && rm raytrace_xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".fits.gz
[ -e xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".arf.gz ] && rm xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".arf.gz
punlearn xaarfgen
xaarfgen xrtevtfile=raytrace_xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".fits source_ra="$ra" source_dec="$dec" telescop=XRISM instrume=XTEND emapfile=xa"$OBSID"xtd_"$SUB".expo regmode=RADEC regionfile=xa"$OBSID"xtd_"$SUB"_src_"$reg_number".reg sourcetype=IMAGE rmffile=xa"$OBSID"xtd_"$SUB"_src_"$reg_number".rmf.gz erange="$erange" outfile=xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".arf numphoton="$photon" minphoton="$minphoton" teldeffile=CALDB qefile=CALDB contamifile=CALDB obffile=CALDB fwfile=CALDB onaxisffile=CALDB onaxiscfile=CALDB mirrorfile=CALDB obstructfile=CALDB frontreffile=CALDB backreffile=CALDB pcolreffile=CALDB scatterfile=CALDB mode=h clobber=yes seed=7 imgfile="$image_file"

mv xaarfgen.log arf_beta_"$reg_number".log
mv xaxmaarfgen.log xma_arf_beta_"$reg_number".log

if ! [ -e xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".arf ]
then
    echo "\n The xa${OBSID}xtd_${SUB}_beta_${reg_number}.arf cannot be created."
    echo " please change numphoton or minphoton etc.\n"
    exit
fi

gzip xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".arf
gzip raytrace_xa"$OBSID"xtd_"$SUB"_beta_"$reg_number".fits

done



# Remove unneceaary files

rm coordpnt.log
rm xaarfgen_region.lis
rm arfgencgrid_full.fits
rm arfgencgrid.fits
rm _tmp_*
rm arf_columns.lis
rm arf_data.dat
rm heasim_events.fits
rm heasim_photonlist.fits
rm raytracehk_coldef.lis
rm raytracehk.dat
rm raytracehk.fits
rm sim_arf.fits
rm srcfile.txt
rm *arfregion0.reg


# Move log files

[ -d log ] || mkdir log
mv -f *.log log