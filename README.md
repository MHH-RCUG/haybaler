# Haybaler
### Sophia Poertner, Colin Davenport, 2020-2021

- Combine your Wochenende .bam.txt or reporting output from multiple samples into one easy matrix.
- Create simple heatmaps in R from the results
- Works only for Wochenende https://github.com/MHH-RCUG/Wochenende, not for Kraken (kraken2table projects exist online for that).


### Details
- Wochenende ouputs many files per sample. Every file contains different result types
- Haybaler collates the different files to one file per result type containing all samples
- one file per sample containing all results -> one file per result type containing all samples
- Less files, easy to compare different samples
- The haybaler output is ordered by read count: the organisms/chromosome with the highest read count over all samples will be at the top of the haybaler output file
- per default filters out every chromosome with less than 10 Reads or less than 300 RPMM in every sample
- excluded taxa are saved (with the reason) in a separate file excluded_taxa.csv
- Scripts allow data preparation for heatmaps and further analyses in R.


### Installation via conda
First install miniconda if you have not already done this.
Required libs are listed in the file env.haybaler.yml
```
# first clone haybaler
git clone https://github.com/MHH-RCUG/haybaler
cd haybaler
# Set up a dedicated environment for haybaler
conda env create -f env.haybaler.yml
conda activate haybaler
```


### Usage

```

# Method 1: Copy the scripts haybaler.py, run_haybaler.sh and all the files you want to combine  into one directory. 
# Then run the run_haybaler script 
bash run_haybaler.sh

# Method 2: Configure your haybaler location and run the wochenende_postprocess.sh script from the Wochenende project
# This will run haybaler as one step
```

#### change read_count or RPMM limit
- default: readcount_limit: 10, rpmm_limit = 300
- open file run_haybaler.sh, edit last line
- add –readcount_limit xx
- add –rpmm_limit xx
```
# example with readcount_limit 20 and rpmm_limit 200
python3 haybaler.py -i "$input_files" -p . -op $outputDir  -o haybaler.csv --readcount_limit 20 --rpmm_limit 200
```
 
### Output

The output is in the created output folder, default haybaler_output

Output is a set of CSVs. These combine the results from the original files into once matrix, so you can better compare your samples.


### Refining the output

You can read the output into R for example and do further analyses yourself, or use our heatmap and heattree scripts.
Heatmaps and Heattrees are also generated with the postprocess script.

### Heatmaps
- exclude mouse, human, mito
- heatmap for the top x organisms (default 50 and 200 taxa in version 0.16)
- use both base R heatmap and heatmaply heatmaps
- display raw and square-rooted results

```
# go to the haybaler output file
# Now run the Rscript to create a heatmap, this requires an R installation.
# Because of a bug with heatmaply and conda, no conda enviroment is allowed to be activated
conda deactivate
bash runbatch_heatmaps.sh 
```

### Heattrees
- Get taxonomy
- Needs taxonkit Datasets installed https://bioinf.shenwei.me/taxonkit/#dataset
```
# copy the taxonomy script in the haybaler_output_directory
cd haybaler_output
cp ../haybaler_taxonomy.py ../run_haybaler_tax.sh .

# run the scripts. Needs Haybaler enviroment and R installation 
conda activate haybaler
bash run_haybaler_tax.sh
```
- Create Heattrees for RPMM and bacteria_per_human_cell files
- needs R installation
- Needs Metacoder package installed (makes troubble installing it, better try it before running the script)
- exclude mouse, human and mito
- one heattree for the sums of all sample
- one heattree for each sample with the sums as "background"
- one heattree for each sample without "background"
```
# copy the heattree scripts in the haybaler_output_directory
cd haybaler_output
cp ../create_heattrees.R ../run_heattrees.sh

# run the scripts. Needs R installation 
bash run_heattrees.sh
```
#### change the level (genus or species) for the heattrees

- edit the create_heattrees.R Script

```
# go to line 31
## column with the lineage information. Uncomment one
# wanted_column <- "genus_lineage"
wanted_column <- "species_lineage"
```


