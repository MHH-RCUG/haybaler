# Script to create heattrees out of haybaler output with Metadata. 
# This is just a template, not an automated script, and will need to be modified before use. 
# Author: Sophia Poertner, April - May 2021

## checks if required packages are installed, install if not. Then load all packages

# packages
packages = c("metacoder", "taxa", "dplyr", "tibble", "ggplot2")

# install uninstalled packages
not_installed <- packages[!(packages %in% installed.packages()[ , "Package"])]  # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos="http://cran.rstudio.com/")  # Install not installed packages from cran 

#load packages
invisible(lapply(packages, library, character.only = TRUE))


## read in data

path = "."                             # Path to your data, default current dir "."
directory = "."                        # directory you want the heatmap to be saved in ("." is current dir)
filename <- "RPMM_filt2_heattree.csv"  # name of your file
# change the file and name to the matadata you want
metadata_name <- "/mnt/ngsnfs/gen/rcug_lw/sophias_projekte/haybaler/heat_trees/Metadata.csv"

## column with the lineage information. Uncomment one
# wanted_column <- "genus_lineage"
wanted_column <- "species_lineage"

# check if wanted column exists
input_file <- read.csv(file = paste0(path,"/",filename), sep = "\t")
if(!(wanted_column %in% colnames(input_file))){
  stop("The wanted column for lineage does not exist.")
}


## clean/filter data and create taxmap

# select the wanted column and name it "lineage". Delete unwanted columns and rows
input_file <- cbind(lineage=input_file[[wanted_column]], input_file)
input_file <- input_file %>% select(-(matches("species|chr_length|gc_ref|genus_name|species_name|species_lineage|genus_lineage")))
input_file <- input_file[input_file$lineage != ';;;;;;',]
input_file <- input_file %>% filter_all(all_vars(.!=Inf))

input_taxmap <- parse_tax_data(input_file,
                               class_cols = "lineage", # the column that contains taxonomic information
                               class_sep = ";", # The character used to separate taxa in the classification
                               class_key = "taxon_name"
)


# Grouped with metadata
# !!! make sure that the samples name in "input_file" and "meatadata" are the same !!!

metadata <- read.csv(file = metadata_name, sep = ";")

# shorten names from input_file, just as test (can be done before reading data in or in R)
colnames(input_taxmap$data$tax_data) <- str_extract(colnames(input_taxmap$data$tax_data), "X19_....|Blank_......|lineage|taxon_id")

col_names <- colnames(input_taxmap$data$tax_data)  # sample names in input_file
metadata <- metadata[metadata$X %in% col_names ,]  # filter out samples which are not in the input_file
group <- "delivery"  # column name of group

# calc abundance for subgroups in group
input_taxmap$data$tax_abund_grouped <- calc_taxon_abund(input_taxmap, "tax_data",
                                                        cols = metadata$X,
                                                        groups = metadata[[group]])
# get the names of the subgroups
sub_groups <- colnames(input_taxmap$data$tax_abund_grouped[, -1])

# create heattree for each subgroup
for (sub_group in sub_groups){
  output_pdf = paste0(filename,"_",group,"_",sub_group,"_heattree.pdf")
  plot <- heat_tree(input_taxmap,
                    node_label = input_taxmap$taxon_names(),
                    node_size = input_taxmap$data$tax_abund_grouped[[sub_group]],
                    node_color = input_taxmap$data$tax_abund_grouped[[sub_group]],
                    node_color_axis_label = "RPMM"
  )
  ggsave(output_pdf, plot=plot, device = "pdf")
}
