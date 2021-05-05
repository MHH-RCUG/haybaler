#!/bin/bash
# Sophia Poertner
# Run haybaler https://github.com/MHH-RCUG/haybaler/

echo "Starting Haybaler taxonomy"

# set directory to get the haybaler heatmaps scripts from
# use default directory if no argument ($1) given

# Users: change this to your haybaler path
haybaler_directory="/mnt/ngsnfs/tools/dev/haybaler/" 


# Users: don't modify this section
if [ -z "$1" ]
then
  haybaler_dir=$haybaler_directory
else
  haybaler_dir="$1"
fi

outputDir=haybaler_output
if [ ! -d $outputDir ]
then
    echo "INFO: Creating directory:" $outputDir
    mkdir $outputDir
fi


# for samples
count=`ls -1 *haybaler.csv 2>/dev/null | wc -l`
if [ $count != 0 ]
    then
    for csv in $(ls *haybaler.csv)
    do
      python3 haybaler_taxonomy.py -i $csv -p .
    done
fi

# uncomment the next section for testing references
count=`ls -1 *fa*.fai 2>/dev/null | wc -l`
if [ $count != 0 ]
    then
    for fai in $(ls *fa*.fai)
    do
      python3 haybaler_taxonomy.py -i $fai -p . -t True
    done
fi

##### taxonomy ######
# python3 haybaler_taxonomy.py  -i 2021_02_human_bact_fungi_vir_masked.fa.fai -p /lager2/rcug/seqres/metagenref/

