#!/bin/bash
## run Rscript on a server with R, eg hpc04,5,6, hpc-rc08, etc

# Prepare data for R heatmaps
# exclude GC, ref length, any host chr etc (all distort heatmaps)
# Sophia Poertner, Colin Davenport, 2020-2021


prepare_files () {
  echo "INFO: Preparing files for R heatmap creation"
  for infile in *haybaler.csv
        do
        echo "Running on " $infile

        #exclude mouse, human, mito
        grep -v "^chr" $infile | grep -v "^1_1_1" > $infile.filt1.csv

        # using tab delimiters, cut a max of 2000 columns out excluding cols 2-3. Also restrict to number_of_taxa lines
        cut -f1,4-2000 $infile.filt1.csv | head -n $1 > $infile.filt2.csv

        # remove _complete_genome from labels
        sed "s/_complete_genome//g" $infile.filt2.csv > $infile.filt.heatmap.csv

        # cleanup: remove temp filt1.csv and filt2.csv files
        rm $infile.filt1.csv
        rm $infile.filt2.csv


  done
}


create_heatmaps () {
  echo "INFO: Starting batch heatmap creation"

# check for rscript, exit if unavailable
# rscript_bin="/usr/bin/Rscript"
# get rscript_bin form config_yaml. Run setup.sh and restart session
if [[ ! -f $rscript_bin ]]
        then
        echo "INFO: Rscript binary not found, aborting. Could not find this, is R installed? " $rscript_bin
        exit
fi
echo "INFO: Using rscript binary: " $rscript_bin

# create heatmaps for each heatmap.csv file
for heatmapcsv in *.heatmap.csv
        do
        echo "INFO: Creating heatmap for file: $heatmapcsv"
        # run local
        $rscript_bin create_heatmap.R $heatmapcsv
done
}

# Setup conda and directories
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)

# Create heatmaps with 50 taxa
if [[ ! -d "top_50_taxa" ]]
        then
	mkdir top_50_taxa
fi
# Run bash function for 50 taxa
prepare_files 50
create_heatmaps
if [[ $count_pdf != 0 ]]
    then
    mv *heatmap*.pdf top_50_taxa
fi
if [[ $count_html != 0 ]]
    then
    mv *heatmap*.html top_50_taxa
fi


# Create heatmaps with 200 taxa
if [[ ! -d "top_200_taxa" ]]
	then
	mkdir top_200_taxa
fi
# Run bash function for 200 taxa
prepare_files 200
create_heatmaps
count_html=`ls -1 *heatmap*.html 2>/dev/null | wc -l`
count_pdf=`ls -1 *heatmap*.pdf 2>/dev/null | wc -l`
if [[ $count_pdf != 0 ]]

    then
    mkdir top_200_taxa
fi
# Run bash function for 200 taxa
prepare_files 200
create_heatmaps
count_html=`ls -1 *heatmap*.html 2>/dev/null | wc -l`
count_pdf=`ls -1 *heatmap*.pdf 2>/dev/null | wc -l`
if [[ $count_pdf != 0 ]]
    then
    mv *heatmap*.pdf top_200_taxa
fi
if [[ $count_html != 0 ]]
    then
    mv *heatmap*.html top_200_taxa
fi

  echo "INFO: Script completed"
#else - SP TODO - this was throwing a syntax error
#  echo "no input files found for heatmaps creation. Needs *haybaler.csv as input"
#fi
