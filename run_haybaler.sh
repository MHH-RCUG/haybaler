#!/bin/bash

echo "Starting Haybaler"

outputDir=haybaler_output
if [ ! -d $outputDir ]
then
    echo "INFO: Creating directory:" $outputDir
    mkdir $outputDir
fi

for csv in $(ls *.bam*.csv)
do
  python3 haybaler.py -i "$csv" -p . -op $outputDir  -o haybaler.csv
done

for csv in $(ls *.bam*.txt)
do
  python3 haybaler.py -i "$csv" -p . -op $outputDir  -o haybaler.csv
done
