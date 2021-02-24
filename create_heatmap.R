# R Script to create a heatmap after preparing the data with haybaler and the preparing script
# Sophia Pörtner, Dec. 2020
# heatmap requires base R
# the interactive heatmaply variant requires heatmaply # install.packages(heatmaply)
usage = "Usage: Rscript create_heatmap.R infile.csv"
usage

version = "0.11"
version

# changelog
# 0.11  change colours
# 0.10  add heatmaply interactive heatmap variant

args <- commandArgs()
#args
file<-args[6]   # name of file of your data

cmd_msg = "File to process: "
cmd_msg
file

# Variables
path = "."              # Path to your data, default current dir "."
directory = "."         # directory you want the heatmap to be saved in ("." is current dir)
outputname = paste0(file,"_","heatmap.pdf")
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

heatmap(
        your_data_2, 
        cexRow = size_rows, 
        cexCol = size_columns
        )

# save the heatmap
setwd(directory)

pdf(outputname, width=16, height=8)
heatmap(
        your_data_2, 
        cexRow = size_rows, 
        cexCol = size_columns
        )
dev.off()

# create an interactive HTML heatmap with heatmaply
# Other color options see https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html
# We don't want a divergent colour palette which accentuates deviations from the mean

# install.packages(heatmaply)
library(heatmaply)

heatmaply(
        your_data_2, 
        #color = YlGn,
        color = Blues,
        #color = cool_warm,
        #file=c(output_html, output_png)  #png needs another missing dependency
        file=c(output_html)
)


