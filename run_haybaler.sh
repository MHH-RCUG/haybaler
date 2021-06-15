#!/bin/bash
# Sophia Poertner
# Run haybaler https://github.com/MHH-RCUG/haybaler/

version="0.21, June 2021"
echo "Starting Haybaler, run_haybaler.sh version" $version

outputDir=haybaler_output
if [[ ! -d $outputDir ]]
then
    echo "INFO: Creating directory:" $outputDir
    mkdir $outputDir
fi

# get conda_env and haybaler_dir from conifg_yaml. Run setup.sh and restart session
. $conda_env
conda activate haybaler

cp $haybaler_dir/*.py $outputDir
cp $haybaler_dir/runbatch_heatmaps.sh $outputDir
cp $haybaler_dir/create_heatmap.R $outputDir
cp $haybaler_dir/run_heattrees.sh $outputDir
cp $haybaler_dir/run_haybaler_tax.sh $outputDir

input_files=""

# Only run for *bam*.csv if files exist in current dir
count=$(ls -1 *.bam*.csv 2>/dev/null | wc -l)
if [[ $count != 0 ]]
    then
    for csv in ls *.bam*.csv
    do
      input_files="$input_files;$csv"
    done
fi


# Only run for *bam*.txt if files exist in current dir
count=$(ls -1 *.bam*.txt 2>/dev/null | wc -l)
if [[ $count != 0 ]]
    then
    for csv in *.bam*.txt
    do
      input_files="$input_files;$csv"
    done
fi

#python3 haybaler.py -i "$input_files" -p . -op $outputDir  -o haybaler.csv
#testing
python3 haybaler.py -i "$input_files" -p . -op $outputDir  -o haybaler.csv --readcount_limit 1 --rpmm_limit 5

