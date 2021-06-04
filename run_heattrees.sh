#!/bin/bash
# Input data: requires tables output by haybaler_taxonomy.R https://github.com/MHH-RCUG/haybaler
# First run Wochenende, then wochenende_postprocess.sh. (or haybaler_taxonomy.sh manually)
# exclude GC, ref length, any host chr etc
# Sophia Poertner, Colin Davenport, 2020-2021
# Installation: First prepare server with R packages, see https://github.com/MHH-RCUG/haybaler
# Usage: bash run_heattrees.sh


# Run only on certain server
server1=hpc-bc15-07
server2=hpc-bc15-12
server3=hpc06


if [[ $(hostname) == $server1 ]]
        then
        echo "INFO: Found hostname $server1. OK. Will attempt to run heat trees here"
elif [[ $(hostname) == $server2 ]]
        then
        echo "INFO: Found hostname $server2. OK. Will attempt to run heat trees here"
else
        echo "INFO: Can only run heat trees on server where heat-trees dependencies are installed, eg. $server1. We can't run heat trees here!"
fi


set_ulimits () {
        echo "INFO: trying to set ulimits higher from 8GB cstack to 32GB cstack"
        ulimit -s 32384
}

prepare_files () {
  echo "INFO: Preparing files for R heatmap creation"
  for infile in {RPMM,bacteria_per_human_cell}*haybaler_taxa.csv
        do
        echo "Running on " "$infile"

        #exclude mouse, human, mito
        grep -v "^chr" "$infile" | grep -v "^1_1_1" > "${infile%_haybaler_taxa.csv}"_filt1_heattree.csv

        # using tab delimiters, cut a max of 2000 columns out excluding cols 2-3.
        cut -f1,4-2000 "${infile%_haybaler_taxa.csv}"_filt1_heattree.csv  > "${infile%_haybaler_taxa.csv}"_filt2_heattree.csv
  done
}


create_heattrees () {
  echo "INFO: Starting batch heat-tree creation"

    # check for rscript, exit if unavailable
    rscript_bin="/usr/bin/Rscript"
    if [[ ! -f $rscript_bin ]]
            then
            echo "INFO: Rscript binary not found, aborting. Could not find this, is R installed? " $rscript_bin
            exit
    fi
    echo "INFO: Using rscript binary: " $rscript_bin


    # create heat-tree for each heatmap.csv file
    # check if correct files present, else do not run the script, which would create empty files.
    count=$(ls -1 *filt2_heattree.csv 2>/dev/null | wc -l)
    if [[ $count != 0 ]]
        then
        for heattreecsv in {RPMM,bacteria_per_human_cell}*filt2_heattree.csv
            do
            echo "INFO: Creating heat-tree for file: $heattreecsv"
            # run local
            $rscript_bin create_heattrees.R "$heattreecsv"
        done
    else
        echo "Could not find heat-tree input files of type filt2_heattree.csv in the directory"
    fi
}


echo "INFO: run this script only on: $server1 or $server2"

# Actually run functions
set_ulimits
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

# move heat-trees in directories
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

bphc_count_background_pdf=$(ls -1 bacteria_per_human_cell*background_heattree.pdf 2>/dev/null | wc -l)
if [[ $bphc_count_background_pdf != 0 ]]
    then
    mv bacteria_per_human_cell*background_heattree.pdf bphc_background_heattrees
fi

echo "INFO: Heat tree script completed"
