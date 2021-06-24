#!/bin/bash
# Sophia Poertner
# Run haybaler https://github.com/MHH-RCUG/haybaler/

version="0.22, June 2021"
echo "Starting Haybaler, run_haybaler.sh version" $version

outputDir=haybaler_output
if [[ ! -d $outputDir ]]
then
    echo "INFO: Creating directory:" $outputDir
    mkdir $outputDir
fi

# Setup config
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)
# Setup conda and directories
. $CONDA_SH_PATH
conda activate $HAYBALER_CONDA_ENV_NAME

cp $HAYBALER_DIR/*.py $outputDir
cp $HAYBALER_DIR/runbatch_heatmaps.sh $outputDir
cp $HAYBALER_DIR/create_heatmap.R $outputDir
cp $HAYBALER_DIR/run_heattrees.sh $outputDir
cp $HAYBALER_DIR/create_heattrees.R $outputDir
cp $HAYBALER_DIR/run_haybaler_tax.sh $outputDir

input_files=""

# Only run for *bam*.csv if files exist in current dir
count=$(ls -1 *.bam*.csv 2>/dev/null | wc -l)
if [[ $count != 0 ]]
    then
    for csv in *.bam*.csv
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

python3 haybaler.py -i "$input_files" -p . -op $outputDir  -o haybaler.csv
# for pipeline testing only!!
#python3 haybaler.py -i "$input_files" -p . -op $outputDir  -o haybaler.csv --readcount_limit 1 --rpmm_limit 10
