#!/bin/bash
## run Rscript on a server with R, eg hpc04,5,6, hpc-rc08, etc
# Colin Davenport, Jan 2021

echo "INFO: Starting batch heatmap creation using create_heatmap.R"

rscript_bin=/usr/bin/Rscript
echo "Using rscript binary: " $rscript_bin

for heatmapcsv in `ls *.heatmap.csv`
        do
	echo "Creating heatmap for file: $heatmapcsv"
        # run local
	$rscript_bin create_heatmap.R $heatmapcsv
done


echo "INFO: Completed run"
