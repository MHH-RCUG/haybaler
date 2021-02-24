#!/bin/bash
## run Rscript on a server with R, eg hpc04,5,6, hpc-rc08, etc
# Colin Davenport, Jan 2021

echo "INFO: Starting batch heatmap creation using create_heatmap.R"

# check for rscript, exit if unavailable
rscript_bin="/usr/bin/Rscript"
if [[ ! -f $rscript_bin ]]
        then
        echo "INFO: Rscript binary not found, aborting. Could not find this, is R installed? " $rscript_bin
        exit
fi
echo "INFO: Using rscript binary: " $rscript_bin

# create heatmaps for each heatmap.csv file
for heatmapcsv in `ls *.heatmap.csv`
        do
	echo "INFO: Creating heatmap for file: $heatmapcsv"
        # run local
	$rscript_bin create_heatmap.R $heatmapcsv
done


echo "INFO: Completed run"
