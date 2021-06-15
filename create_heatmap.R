# R Script to create a heatmap after preparing the data with haybaler and the preparing script
# Sophia Pörtner, Dec. 2020
# heatmap requires base R
# the interactive heatmaply variant requires heatmaply # install.packages(heatmaply)

usage = "Usage: Rscript create_heatmap.R infile.csv"
usage

version = "0.16"
version

# changelog
# 0.16  remove code which generated the unnamed Rplots.pdf
# 0.15  add simple scale legend on base R heatmap
# 0.14  check if required packages are installed, install them if not.
# 0.13  sqrt as default transform. log(0) is -Inf, log(0.6) is negative, also bact_per_hum_cell too close to 0 to simply add 1
# 0.12  change colours cool_warm blue to red default
# 0.11  change colours
# 0.10  add heatmaply interactive heatmap variant


## check if required packages are installed, install if not. Then load all packages

# packages
packages = c("heatmaply", "RColorBrewer")

# install uninstalled packages
not_installed <- packages[!(packages %in% installed.packages()[ , "Package"])]  # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos="http://cran.rstudio.com/")  # Install not installed packages from cran 

#load packages
invisible(lapply(packages, library, character.only = TRUE))


args <- commandArgs()
# args
file<-args[6]   # name of file of your data

cmd_msg = "File to process: "
cmd_msg
file

# Variables
path = "."              # Path to your data, default current dir "."
directory = "."         # directory you want the heatmap to be saved in ("." is current dir)
# Names are further changed below
output_pdf = paste0(file,"_","heatmap.pdf")
output_html = paste0(file,"_","heatmaply.html")
output_png = paste0(file,"_","heatmaply.png")

# Read the Data
data = sprintf("%s/%s", path, file)
your_data = read.csv(data, sep="\t", header=TRUE, row.names = 1)

# Make your Data numeric so the heatmap function can work with it
your_data_2 = as.matrix(your_data)

# create the heatmap with base R
# Change the size of the Rows (Species), and Columns (Samples) if needed
# Try a few sizes and run this script again until the heatmap looks as you like it
size_rows <- 0.6
size_columns <- 0.6

# 50 taxa is the reference. The proportion to 50 taxa is calculate to adjust the size of the pdf 
num_taxa <- nrow(your_data_2)
heatmap_size <- num_taxa / 50
if(heatmap_size < 1){
        heatmap_size <- 1
}

# Set the working directory
setwd(directory)

# raw heatmap
output_pdf = paste0(file,"_","heatmap1_raw.pdf")
pdf(output_pdf, width=16*heatmap_size, height=8*heatmap_size)
heatmap(
        your_data_2, 
        cexRow = size_rows, 
        cexCol = size_columns, 
        col = colorRampPalette(brewer.pal(8, "Oranges"))(25)
) 
#attempt to add scale/ legend
legend(x="bottomright", legend=c("min", "ave", "max"),
       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))

dev.off()

# square root normalized
output_pdf = paste0(file,"_","heatmap2_sqrt.pdf")
pdf(output_pdf, width=16*heatmap_size, height=8*heatmap_size)
heatmap(
        sqrt(your_data_2),
        cexRow = size_rows,
        cexCol = size_columns,
        col= colorRampPalette(brewer.pal(8, "Oranges"))(25)
)
#attempt to add scale/ legend
legend(x="bottomright", legend=c("min", "ave", "max"), 
       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))

dev.off()

# log heatmap
output_pdf = paste0(file,"_","heatmap1_log2.pdf")
pdf(output_pdf, width=16*heatmap_size, height=8*heatmap_size)
heatmap(
        log2(your_data_2 + 1), 
        cexRow = size_rows, 
        cexCol = size_columns, 
        col = colorRampPalette(brewer.pal(8, "Oranges"))(25)
) 
#attempt to add scale/ legend
legend(x="bottomright", legend=c("min", "ave", "max"),
       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))

dev.off()

# create an interactive HTML heatmap with heatmaply
# Other color options see https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html
# We don't want a divergent colour palette which accentuates deviations from the mean


# raw
output_html = paste0(file,"_","heatmaply1_raw.html")
heatmaply(
        your_data_2, 
        #color = YlGn,
        #color = Blues,
        color = cool_warm,
        #file=c(output_html, output_png)  #png needs another missing dependency
        file=c(output_html)
)

# sqrt
output_html = paste0(file,"_","heatmaply2_sqrt.html")
heatmaply(
        sqrt(your_data_2), 
        #color = YlGn,
        #color = Blues,
        color = cool_warm,
        #file=c(output_html, output_png)  #png needs another missing dependency
        file=c(output_html)
)

# log2
output_html = paste0(file,"_","heatmaply1_log2.html")
heatmaply(
        log2(your_data_2 + 1), 
        #color = YlGn,
        #color = Blues,
        color = cool_warm,
        #file=c(output_html, output_png)  #png needs another missing dependency
        file=c(output_html)
)

# percentize (like percentage ranks)
output_html = paste0(file,"_","heatmaply3_percentize.html")
heatmaply(
        percentize(your_data_2), 
        #color = YlGn,
        #color = Blues,
        color = cool_warm,
        #file=c(output_html, output_png)  #png needs another missing dependency
        file=c(output_html)
)

