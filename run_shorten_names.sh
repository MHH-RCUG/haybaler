#!/bin/bash


mkdir short_names

dir="$(pwd)"
file=$1


echo "Running names.R for ${file}, output is in /short"
if ! Rscript names_add.R $dir $file 
then
	echo "This script needs to be run on an R server"
else
	echo "Finished"
fi

