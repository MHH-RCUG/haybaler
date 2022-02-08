# Haybaler
### Sophia Poertner, Colin Davenport, Lisa Hollstein 2020-2022

- Combine your Wochenende .bam.txt or reporting output from multiple samples into one easy matrix.
- Create heatmaps and heat trees in R from the results
- Works only for Wochenende https://github.com/MHH-RCUG/Wochenende, not for Kraken (kraken2table projects exist online for that).


### Details
- Wochenende outputs many files per sample. Every file contains different result types
- Haybaler collates the different files to one file per result type containing all samples
- one file per sample containing all results -> one file per result type containing all samples
- Less files, easy to compare different samples
- The haybaler output is ordered by read count: the organisms/chromosome with the highest read count over all samples will be at the top of the haybaler output file
- per default filters out every chromosome with less than 10 Reads or less than 300 RPMM in every sample
- excluded taxa are saved (with the reason) in a separate file excluded_taxa.csv
- Scripts allow data preparation for heatmaps and further analyses in R.


### Installation via conda
First install miniconda if you have not already done this. Use mamba instead of conda if you like faster installs (follow the mamba install instructions here https://github.com/mamba-org/mamba )
Required libs are listed in the file env.haybaler.yml
```
# first clone haybaler
git clone https://github.com/MHH-RCUG/haybaler
cd haybaler
# Set up a dedicated environment for haybaler
conda env create -f env.haybaler.yml
# OR mamba
mamba env create -f env.haybaler.yml
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

#### Set lower abundance thresholds, change readcount_limit or RPMM_limit from defaults >10 and >300
- We set configurable lower limits to detected species so as to avoid false positive detection of various taxa. Otherwise these show up in tables and heatmaps and appear to be present, despite very weak evidence for them.
- defaults: readcount_limit: 10 (taxa with less than 10 reads attributed are excluded)
- defaults: rpmm_limit = 300 (taxa with an RPMM value of less than 300 are excluded)
- open file run_haybaler.sh, edit last line
- add --readcount_limit xx
- add --rpmm_limit xx
```
# example with readcount_limit 20 and rpmm_limit 200
python3 haybaler.py -i "$input_files" -p . -op $outputDir  -o haybaler.csv --readcount_limit 20 --rpmm_limit 200
```
 
### Output


The output is in the created output folder (default: haybaler_output)

Output is a set of CSVs. These combine the results from the original files into a single matrix, so you can better compare your samples. Furthermore, heatmaps and heattrees are created, provided you have an R installation set up correctly.


### Refining the output

You can read the output into R for example and do further analyses yourself, or use our heatmap and heattree scripts.
Heatmaps and Heattrees are also generated with the `wochenende_postprocess.sh` script.

### Heatmaps
- exclude mouse, human, mitos
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
- saves every chomosome where it is not possible to identify genus/species in no_genus/species_found_for_chromosomes.csv
- more details about the script: haybaler_taxa.md

```
# copy the taxonomy script in the haybaler_output_directory
cd haybaler_output
cp ../haybaler_taxonomy.py ../run_haybaler_tax.sh .

# run the scripts. Needs Haybaler enviroment and R installation 
conda activate haybaler
bash run_haybaler_tax.sh
```
- Create Heattrees for RPMM and bacteria_per_human_cell files
- for more information about heat trees: https://github.com/grunwaldlab/metacoder
- needs R installation
- Needs Metacoder package installed (may cause trouble when installing, so install before haybaler)
- exclude mouse, human and mito
- one heattree for the sums of all samples
- one heattree for each sample with the sums as "background"
- one heattree for each sample without "background"
- do not create heattrees for empty samples. They are saved in the file `empty_samples.txt`
```
# copy the heattree scripts in the haybaler_output_directory
cd haybaler_output
cp ../create_heattrees.R ../run_heattrees.sh .

# run the scripts. Needs R installation 
bash run_heattrees.sh
```
#### change the level (genus or species) for the heattrees

- edit the create_heattrees.R script

```
# go to line 31
## column with the lineage information. Uncomment one
# wanted_column <- "genus_lineage"
wanted_column <- "species_lineage"
```

#### change node/edge/label sizes
- edit the create_heatrees.R script

```
# go to line 44
# Aesthetic arguments for Heat tree creation(size or text), change as you like 
# for more infomation: https://rdrr.io/cran/metacoder/man/heat_tree.html
node_size_range <- c(0.01, 0.05)
edge_size_range <- c(0.005, 0.025)
node_color_axis_label <- "absolute abundance in\n sample (color)"
node_size_axis_label <- "\nabsolute abundance among\n all samples (size)"
title_size <- 0.05 
```


Contributions

@poer-sophia - main author, taxonomy, heatmaps, testing

@colindaven - concept, code review, testing

@irosenboom heat-trees, testing

@LisaHollstein Code improvements, maintainence, testing, code review, docs

