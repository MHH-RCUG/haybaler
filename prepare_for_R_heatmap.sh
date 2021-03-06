#!/bin/bash
# Colin, Nov 2020
# Prepare data for R heatmaps
# exclude GC, ref length, any host chr etc (all distort heatmaps)

echo "Usage: bash prepare_for_R_heatmap.sh"

# change this to the number of taxa you want in the heatmap (just uses head -n X). 
# 50 seems the max number where all labels on the y axis are still visible
number_of_taxa=50

for infile in `ls *haybaler.csv`
        do
	echo "Running on " $infile

        #exclude mouse, human, mito
        grep -v "^chr" $infile | grep -v "^1_1_1" > $infile.filt.csv

	# using tab delimiters, cut a max of 200 columns out excluding cols 2-3. Also restrict to number_of_taxa lines
	cut -f1,4-200 -d"	" $infile.filt.csv | head -n $number_of_taxa > $infile.filt2.csv

	# remove _complete_genome from labels
	sed "s/_complete_genome//g" $infile.filt2.csv > $infile.filt.heatmap.csv

	# shorten TODO names to 20 chars? awk ?
	#$infile.filt.heatmap.csv > $infile.filt.heatmap.csv

done


echo "INFO: Script completed"
