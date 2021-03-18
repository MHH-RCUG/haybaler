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

### Output

The output is in the created output folder, default haybaler_output

Output is a set of CSVs. These combine the results from the original files into once matrix, so you can better compare your samples.


### Refining the output

You can read the output into R for example and do further analyses yourself, or use our heatmap scripts
- exclude mouse, human, mito
- heatmap for the top x organisms (default 50 taxa in version 0.16)
- use both base R heatmap and heatmaply heatmaps
- display raw and square-rooted results

```
# go to the haybaler output file
# First prepare the data for a heatmap
prepare_for_R_heatmap.sh  

# Now run the Rscript to create a heatmap, this requires an R installation.
# Because of a bug with heatmaply and conda, no conda enviroment is allowed to be activated
conda deactivate
runbatch_create_heatmap.sh  
```
