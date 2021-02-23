# Haybaler
### Sophia Poertner, 2020-2021

- Combine your Wochenende .bam.txt or reporting output from multiple samples into one easy matrix.
- Create simple heatmaps in R from the results
- Works only for Wochenende https://github.com/MHH-RCUG/Wochenende, not for Kraken (kraken2table projects exist online for that).


### Installation via conda
First install miniconda if you have not already done this.
Required libs are listed in the file env.haybaler.yml
```
# first clone haybaler
git clone xxx
cd haybaler
# Set up a dedicated environment for haybaler
conda env create -f env.haybaler.yml
conda activate haybaler
```


### Usage

```

# Copy the scripts haybaler.py, run_haybaler.sh and all the files you want to combine  into one directory. 
OR
# Run the wochenende_postprocess.sh script from the Wochenende project

# Then run the run_haybaler script 
bash run_haybaler.sh
```

### Output

The output is in the created output folder, default haybaler_output

Output is a set of CSVs. These combine the results from the original files into once matrix, so you can better compare your samples.


### Refining the output

You can read the output into R for example and do further analyses yourself, or use our heatmap scripts

# First prepare the data for a heatmap
prepare_for_R_heatmap.sh  

# Now run the Rscript to create a heatmap, this requires an R installation.
runbatch_create_heatmap.sh  
