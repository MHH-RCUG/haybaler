# Script to create heattrees out of haybaler output.
# one heattree for the sums of all sample
# one heattree for each sample with the sums as "background"
# one heattree for each sample without "background"
# Author: Sophia Poertner, April - June 2021


## check if required packages are installed, install if not. Then load all packages

# packages
packages = c("metacoder", "taxa", "dplyr", "tibble", "ggplot2","stringr","RColorBrewer")

# install uninstalled packages
not_installed <- packages[!(packages %in% installed.packages()[ , "Package"])]  # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos="http://cran.rstudio.com/")  # Install not installed packages from cran 

#load packages
invisible(lapply(packages, library, character.only = TRUE))

# Check and set ulimits in R (TEST). Bash setting might be sufficient
#check Cstack size
#Cstack_info()["size"]
#set 32GB
system("ulimit -s 32384")
#check parameter is bigger
#Cstack_info()["size"]


# args
args <- commandArgs()
filename <-args[6]   # name of file of your data

cmd_msg = "File to process: "
cmd_msg
filename

path = "."              # Path to your data, default current dir "."
directory = "."         # directory you want the heatmap to be saved in ("." is current dir)

## column with the lineage information. Uncomment one
# wanted_column <- "genus_lineage"
wanted_column <- "species_lineage"

# Aesthetic arguments for Heat tree creation(size or text), change as you like 
# for more infomation: https://rdrr.io/cran/metacoder/man/heat_tree.html
node_size_range <- c(0.01, 0.05)
edge_size_range <- c(0.005, 0.025)
node_color_axis_label <- "absolute abundance in\n sample (color)"
node_size_axis_label <- "\nabsolute abundance among\n all samples (size)"
title_size <- 0.05


#####
# file preparation 
#####

# check if wanted column exists
input_file <- read.csv(file = filename, sep = "\t")
if(!(wanted_column %in% colnames(input_file))){
  stop("The wanted column for lineage does not exist.")
}

# select the wanted column and name it "lineage". Delete unwanted columns
input_file <- cbind(lineage=input_file[[wanted_column]], input_file)
input_file <- input_file %>% select(-(matches("species|chr_length|gc_ref|genus_name|species_name|species_lineage|genus_lineage")))
input_file <- input_file[input_file$lineage != ';;;;;;',]
input_file <- input_file %>% filter_all(all_vars(.!=Inf))

input_taxmap <- parse_tax_data(input_file,
                       class_cols = "lineage", # the column that contains taxonomic information
                       class_sep = ";", # The character used to separate taxa in the classification
                       class_key = "taxon_name"
)

data_type <- str_replace(filename, "_filt2_heattree.csv", "")   # RPMM or bacteria_per_human_cell


#####
# all samples sums
#####

num_cols_data <- ncol(input_taxmap$data$tax_data)
input_taxmap$data$tax_abund <- calc_taxon_abund(input_taxmap, "tax_data", cols= 3:num_cols_data)
num_cols_abund <- ncol(input_taxmap$data$tax_abund)
input_taxmap$data$tax_abund$sums <- rowSums(input_taxmap$data$tax_abund[, 3:num_cols_abund])

set.seed(1)
plot <- heat_tree(input_taxmap,
                  node_label = taxon_names(input_taxmap),
                  node_size = input_taxmap$data$tax_abund$sums,
                  node_color =  input_taxmap$data$tax_abund$sums,
                  node_size_range = node_size_range,
                  edge_size_range = edge_size_range,
                  node_color_axis_label = "summarized abundance\n from all samples",
                  title = data_type,
                  title_size = title_size
)
output_pdf = paste0(filename,"_all_samples_heattree.pdf")
ggsave(output_pdf, plot=plot, device = "pdf")


#####
# each sample with all as background
#####

empty_samples <- ""

# get max sample sum
max_col_sum <- max(colSums(input_file[,-1]))

samples <- colnames(input_file[, -1])
for(sample in samples) {
  
  one_sample <- input_file %>% 
    select(lineage, all_of(sample))
  one_sample <- one_sample[one_sample[[sample]] > 0,]
  
  #check if dataframe is empty. skip iteration if so, else a Error would occur
  if (dim(one_sample)[1] == 0) {
    print(paste0("empty sample. Skip heattree creation for sample ", sample))
    empty_samples <- paste(empty_samples, sample, sep="\n")
    next
  }
  
  set.seed(1)
  plot <- heat_tree(input_taxmap,
            node_label = taxon_names(input_taxmap),
            node_size = input_taxmap$data$tax_abund[[sample]],
            node_color =  input_taxmap$data$tax_abund[[sample]],
            node_size_range = node_size_range,
            node_size_interval = range(c(0,max_col_sum), na.rm = TRUE, finite = TRUE),
            edge_size_range = edge_size_range,
            edge_size_interval = range(c(0,max_col_sum), na.rm = TRUE, finite = TRUE),
            node_color_axis_label = node_color_axis_label,
            node_size_axis_label = node_size_axis_label,
            title = data_type,
            title_size = title_size
  )
  output_pdf = paste0(filename,"_",sample,"_background_heattree.pdf")
  ggsave(output_pdf, plot=plot, device = "pdf")
}

write(empty_samples, paste0(path, "/empty_samples.txt"))


#####
# each sample without background
#####

empty_samples <- ""

# get max sample sum
max_col_sum <- max(colSums(input_file[,-1]))


samples <- colnames(input_file[, -1])
for(sample in samples) {
  one_sample <- input_file %>% 
  select(lineage, all_of(sample))
  one_sample <- one_sample[one_sample[[sample]] > 0,]
  
  
  #check if dataframe is empty. skip iteration if so, else a Error would occur
  if (dim(one_sample)[1] == 0) {
    print(paste0("empty sample. Skip heattree creation for sample ", sample))
    empty_samples <- paste(empty_samples, sample, sep="\n")
    next
  }
  
  one_sample_taxmap <- parse_tax_data(one_sample,
                                    class_cols = "lineage", # the column that contains taxonomic information
                                    class_sep = ";", # The character used to separate taxa in the classification
                                    class_key = "taxon_name")
  one_sample_taxmap$data$tax_abund <- calc_taxon_abund(one_sample_taxmap, "tax_data", cols=3)
  
  set.seed(1)
  plot <- heat_tree(one_sample_taxmap,
            node_label = one_sample_taxmap$taxon_names(),
            node_size = one_sample_taxmap$data$tax_abund[[sample]],
            node_color = one_sample_taxmap$data$tax_abund[[sample]],
            node_size_range = node_size_range,
            node_size_interval = range(c(0,max_col_sum), na.rm = TRUE, finite = TRUE),
            edge_size_range = edge_size_range,
            edge_size_interval = range(c(0,max_col_sum), na.rm = TRUE, finite = TRUE),
            node_color_axis_label = node_color_axis_label,
            node_size_axis_label = node_size_axis_label,
            title = data_type,
            title_size = title_size
  )
  output_pdf = paste0(filename,"_",sample,"_no_background_heattree.pdf")
  ggsave(output_pdf, plot=plot, device = "pdf")
}

