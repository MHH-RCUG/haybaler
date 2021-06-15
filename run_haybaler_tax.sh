#!/bin/bash
# Sophia Poertner
# Run haybaler taxonomy to attribute NCBI taxonomic lineage to Haybaler output files *haybaler.csv in the current dir. See https://github.com/MHH-RCUG/haybaler/
# Usage: bash run_haybaler_tax.sh


echo "Starting Haybaler taxonomy"

# get conda_env from conifg_yaml. Run setup.sh and restart session
. $conda_env
conda activate haybaler


# for samples
count=`ls -1 *haybaler.csv 2>/dev/null | wc -l`
if [[ $count != 0 ]]
    then
    for csv in $(ls *haybaler.csv)
    do
      python3 haybaler_taxonomy.py -i $csv -p .
    done
fi


# uncomment the next section for testing references (just for developers!)
count=`ls -1 *fa*.fai 2>/dev/null | wc -l`
if [[ $count != 0 ]]
    then
    for fai in $(ls *fa*.fai)
    do
      python3 haybaler_taxonomy.py -i $fai -p . -t True
    done
fi
##### taxonomy ######
# python3 haybaler_taxonomy.py  -i 2021_02_human_bact_fungi_vir_masked.fa.fai -p /mnt/ngsnfs/seqres/metagenref/

