#!/bin/zsh

# Please execute this program using the source command.
# and define OBSID and Sub-string in ~/.zshrc.



# Reprocessing

echo ""
punlearn xapipeline
xapipeline indir="$OBSID" outdir="$OBSID"_xtd_reproc steminputs=xa"$OBSID" stemoutputs=DEFAULT entry_stage=1 exit_stage=3 instrument=xtend verify_input=no create_ehkmkf=no calc_pointing=yes calc_optaxis=yes seed=650504
[ -d log_xtd ] || mkdir log_xtd
mv *.log log_xtd
echo ""
gzip -r -v "$OBSID"_xtd_reproc 2>&1 | tee xapipeline_gzip.log   # The tee command is important.
mv xapipeline_gzip.log log_xtd



# Make directory "analysis_xtd" and Move files to "analysis_xtd"
# This phase is different depending on the observation data.

[ -d analysis_xtd ] || mkdir analysis_xtd
cd analysis_xtd
ln -s ../"$OBSID"_xtd_reproc/xa"$OBSID"xtd_"$SUB"_cl.evt.gz   # The location of the cl.evt.gz file may depend on the observation data.
ln -s ../"$OBSID"_xtd_reproc/xa"$OBSID".ehk.gz   # The location of the ehk.gz file may depend on the observation data.
ln -s ../"$OBSID"_xtd_reproc/xa"$OBSID"xtd_"$SUB".bimg.gz   # The location of the bimg.gz file may depend on the observation data.
cd ..



# Coordinate check and Substitution

[ -e coordinate_xtd.txt ] && rm coordinate_xtd.txt

for par in RA_NOM DEC_NOM PA_NOM
do
    fkeyprint analysis_xtd/xa"$OBSID"xtd_"$SUB"_cl.evt.gz $par outfile=STDOUT | grep -A1 "EXTENSION:    1" | grep "$par" >> coordinate_xtd.txt
done

read RA_NOM DEC_NOM PA_NOM <<< $(python coordinate_xtd.py)   # This sentence is important. (command replacement, here-string, read command, strip method, etc.)
export RA_NOM DEC_NOM PA_NOM
echo "\n RA_NOM = $RA_NOM"
echo " DEC_NOM = $DEC_NOM"
echo " PA_NOM = $PA_NOM"



# Completed this program

echo "\n This \"xtend1.zsh\" has been completed."
echo " Go to next step \"removing anomalous pixels from Xtend data\".\n"