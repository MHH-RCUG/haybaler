#!/bin/bash
## run on server hpc-bc15-07

# Input data: requires tables output by haybaler_taxonomy.R https://github.com/MHH-RCUG/haybaler
# First run Wochenende, then haybaler_taxonomy.sh
# exclude GC, ref length, any host chr etc
# Sophia Poertner, Colin Davenport, 2020-2021
# Usage: bash run_heattrees.sh 


prepare_files () {
  echo "INFO: Preparing files for R heatmap creation"
  for infile in {RPMM,bacteria_per_human_cell}*haybaler_taxa.csv
        do
        echo "Running on " "$infile"

        #exclude mouse, human, mito
        grep -v "^chr" "$infile" | grep -v "^1_1_1" > "${infile%_haybaler_taxa.csv}"_filt1_heattree.csv

        # using tab delimiters, cut a max of 200 columns out excluding cols 2-3.
        cut -f1,4-200 "${infile%_haybaler_taxa.csv}"_filt1_heattree.csv  > "${infile%_haybaler_taxa.csv}"_filt2_heattree.csv
  done
}


create_heattrees () {
  echo "INFO: Starting batch heattree creation"

# check for rscript, exit if unavailable
rscript_bin="/usr/bin/Rscript"
if [[ ! -f $rscript_bin ]]
        then
        echo "INFO: Rscript binary not found, aborting. Could not find this, is R installed? " $rscript_bin
        exit
fi
echo "INFO: Using rscript binary: " $rscript_bin

# create heattree for each heatmap.csv file
for heattreecsv in {RPMM,bacteria_per_human_cell}*filt2_heattree.csv
        do
        echo "INFO: Creating heattree for file: $heattreecsv"
        # run local
        $rscript_bin create_heattrees.R "$heattreecsv"
done
}


echo "INFO: run this script only on hpc-bc15-07!"

prepare_files
create_heattrees

echo "INFO: Cleanup - creating directories and moving files"
# check if directories exist, create them if not
if [[ ! -d "RPMM_background_heattrees" ]]
        then
	mkdir RPMM_background_heattrees
fi

if [[ ! -d "RPMM_no_background_heattrees" ]]
        then
	mkdir RPMM_no_background_heattrees
fi

if [[ ! -d "bphc_background_heattrees" ]]
        then
	mkdir bphc_background_heattrees
fi

if [[ ! -d "bphc_no_background_heattrees" ]]
        then
	mkdir bphc_no_background_heattrees
fi

# move heattrees in directories
RPMM_count_no_background_pdf=$(ls -1 RPMM_*no_background_heattree.pdf 2>/dev/null | wc -l)
if [[ $RPMM_count_no_background_pdf != 0 ]]
    then
    mv RPMM_*no_background_heattree.pdf RPMM_no_background_heattrees
fi

RPMM_count_background_pdf=$(ls -1 RPMM_*background_heattree.pdf 2>/dev/null | wc -l)
if [[ $RPMM_count_background_pdf != 0 ]]
    then
    mv RPMM_*background_heattree.pdf RPMM_background_heattrees
fi


bphc_count_no_background_pdf=$(ls -1 bacteria_per_human_cell*no_background_heattree.pdf 2>/dev/null | wc -l)
if [[ $bphc_count_no_background_pdf != 0 ]]
    then
    mv bacteria_per_human_cell*no_background_heattree.pdf bphc_no_background_heattrees
fi

bphc_count_background_pdf=$(ls -1 bacteria_per_human_cell*background_heattree.pdf 2>/dev/null | wc -)
if [[ $bphc_count_background_pdf != 0 ]]
    then
    mv bacteria_per_human_cell*background_heattree.pdf bphc_background_heattrees
fi

echo "INFO: Heat tree script completed"
