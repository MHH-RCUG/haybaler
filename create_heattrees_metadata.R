# Script to create heattrees out of haybaler output with Metadata. 
# This is just a template, not an automated script, and will need to be modified before use. 
# Author: Sophia Poertner, April - August 2021

## checks if required packages are installed, install if not. Then load all packages

# packages
packages = c("metacoder", "taxa", "dplyr", "tibble", "ggplot2", "tidyr", "rcompanion", "stringr")

# install uninstalled packages
not_installed <- packages[!(packages %in% installed.packages()[ , "Package"])]  # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos="http://cran.rstudio.com/")  # Install not installed packages from cran 

#load packages
invisible(lapply(packages, library, character.only = TRUE))


## variables

path = "/mnt/ngsnfs/gen/rcug_lw/sophias_projekte/haybaler/test_daten_ilona/haybaler_output" # Path to your data, default current dir "."
directory = "."                        # directory you want the heat-tree to be saved in ("." is current dir)
filename <- "RPMM_filt2_heattree.csv"  # name of your file
metadata_name <- "/mnt/ngsnfs/gen/rcug_lw/sophias_projekte/haybaler/heat_trees/Metadata.csv"  # path + name to your metadata
group = "timepoint_sampling"  # must be a column in metadata. A new column that combines two existing ones can be created in line 37-46


# column with the lineage information (species or genus). Uncomment one
# wanted_column <- "genus_lineage"
wanted_column <- "species_lineage"

## read in data
setwd(path)
input_file <- read.csv(file = filename, sep = "\t", check.names = FALSE)
metadata <- read.csv(file = metadata_name, sep = ";", check.names = FALSE)


### create combined groups. Only execute this section if you want to combine groups, else ignore it

group_1 <- "delivery"
group_2 <- "timepoint_sampling"
new_col = paste(group_1, group_2, sep = "_and_")

metadata <- metadata %>%
  unite(!!new_col, c(!!group_1, !!group_2), remove = FALSE)

### end combine groups


# check if wanted column exists
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

col_names <- colnames(input_taxmap$data$tax_data)  # sample names in input_file
metadata <- metadata[metadata$samples %in% col_names ,]  # filter out samples which are not in the input_file



### Grouped heat trees with metadata ###
# !!! make sure that the sample names in "input_file" and "metadata" are the same and name the sample column in metadata "samples"!!!

data_type <- str_replace(filename, "_filt2_heattree.csv", "")   # RPMM or bacteria_per_human_cell

# calc abundance for subgroups in group
input_taxmap$data$tax_abund_grouped <- calc_taxon_abund(input_taxmap, "tax_data",
                                                        cols = metadata$samples,
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
                    node_color_axis_label = data_type
  )
  ggsave(output_pdf, plot=plot, device = "pdf")
}




### test compare groups and heattree matrix ###

input_taxmap$data$tax_abund <- calc_taxon_abund(input_taxmap, "tax_data",
                                                        cols = metadata$samples)

# input_taxmap$data$tax_abund <- input_taxmap$data$tax_abund 

input_taxmap$data$compared_groups <- compare_groups(input_taxmap, data = "tax_abund",
               cols = metadata$sample,
               groups = metadata$timepoint_sampling,
               function(abund_1, abund_2) {
                 log_ratio <- log2(median(abund_1) / median(abund_2))
                 if (is.nan(log_ratio)) {
                   log_ratio <- 0
                 }
                 list(log2_median_ratio = log_ratio,
                      median_abund_1 = median(abund_1),
                      median_abund_2 = median(abund_2),
                      median_diff = median(abund_1) - median(abund_2),
                      mean_diff = mean(abund_1) - mean(abund_2),
                      wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
               }
               )
## notes from marie pust: 
# rcompanion package, wilcoxonR /effect size -> ci=TRUE
# r (zwischen 0 und 100), CI = X - Y
# wenn mit p-Werten, p value correction wichtig für multiple comparisons
# Benjamini-Hochberg correction

group <- "median_diff"  # set the group you want. Can be everything from list function above (but not everything works)
set.seed(1)
heat_tree_matrix(input_taxmap,
                 data = "compared_groups",
                 node_size = n_obs, # number of observations
                 node_label = taxon_names,
                 node_color = median_diff, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 # node_color_interval = c(-3, 3), # The range of `node_color` to display. Test with this and see how the legend looks
                 # edge_color_interval = c(-3, 3), # The range of `edge_color` to display. Test with this and see how the legend looks
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = input_taxmap$data$compared_groups[[group]],
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "differential_heat_tree.pdf") # Saves the plot as a pdf file

