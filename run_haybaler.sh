#!/bin/bash

echo "Starting Haybaler"

# set directory to get the haybaler heatmaps scripts from
# use default directory if no argument ($1) given
if [ -z "$1" ]
then
  haybaler_dir=/mnt/ngsnfs/tools/dev/haybaler/  # default directory
else
  haybaler_dir="$1"
fi

outputDir=haybaler_output
if [ ! -d $outputDir ]
then
    echo "INFO: Creating directory:" $outputDir
    mkdir $outputDir
fi

cp $haybaler_dir/prepare_for_R_heatmap.sh $outputDir
cp $haybaler_dir/runbatch_create_heatmap.sh $outputDir
cp $haybaler_dir/create_heatmap.R $outputDir

for csv in $(ls *.bam*.csv)
do
  python3 haybaler.py -i "$csv" -p . -op $outputDir  -o haybaler.csv
done

for csv in $(ls *.bam*.txt)
do
  python3 haybaler.py -i "$csv" -p . -op $outputDir  -o haybaler.csv
done
