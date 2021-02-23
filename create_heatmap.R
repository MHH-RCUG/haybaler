# R Script to create a heatmap after preparing the data with haybaler and the preparing script
# Sophia Pörtner, Dec. 2020
usage = "Usage: Rscript create_heatmap.R infile.csv"
usage

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

# Read the Data
data = sprintf("%s/%s", path, file)
your_data = read.csv(data, sep="\t", header=TRUE, row.names = 1)

# Make your Data numeric so the heatmap function can work with it
your_data_2 = as.matrix(your_data)

# create the heatmap
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

