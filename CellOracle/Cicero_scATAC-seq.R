### Example on how to prepare the input data prior to building a base GRN using Cicero to extract the cis-regulatory connections between scATAC-seq peaks

#Import libraries 
library(cicero)
library(monocle3)

#make variable containing peak matrix
data_folder <- "snATAC"

# Create a folder to save results
output_folder <- "Ncicero_output"
dir.create(output_folder)

## Load data and make cell Data Set (CDS) object
# Read in matrix data using the Matrix package
indata <- Matrix::readMM(paste0(data_folder, "/matrix.mtx"))

# Binarize the matrix
indata@x[indata@x > 0] <- 1

# Format cell info
cellinfo <- read.table(paste0(data_folder, "/barcodes.tsv"))
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

# Add batch prefix to cell barcodes (rownames) 
#First check rownames are cell barcodes
rownames(cellinfo)
rownames(cellinfo) <- paste("ga22", rownames(cellinfo), sep=".")

# Format peak info
peakinfo <- read.table(paste0(data_folder, "/peaks.bed"))
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# Make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
cell_metadata = cellinfo,
gene_metadata = peakinfo))

# Save as RDS here to merge multiple sample datasets
saveRDS(input_cds, paste0(output_folder,"/40y_input_cds.Rds"))

# Read in all of the sample cds objects 
ga22_input_cds <- readRDS("ga22_input_cds.Rds") # Make as many variables as required
ga24_input_cds <- readRDS("ga24_input_cds.Rds")
m1_input_cds <- readRDS("1m_input_cds.Rds")
m3_input_cds <- readRDS("3m_input_cds.Rds")
m6_input_cds <- readRDS("6m_input_cds.Rds")
m10_input_cds <- readRDS("10m_input_cds.Rds")
y1_input_cds <- readRDS("1y_input_cds.Rds")
y2_input_cds <- readRDS("2yr_input_cds.Rds")
y6_input_cds <- readRDS("6y_input_cds.Rds")
y8_input_cds <- readRDS("8y_input_cds.Rds")
y10_input_cds <- readRDS("10y_input_cds.Rds")
y14_input_cds <- readRDS("14y_input_cds.Rds")
y16_input_cds <- readRDS("16y_input_cds.Rds")
y20_input_cds <- readRDS("20y_input_cds.Rds")
y25_input_cds <- readRDS("25y_input_cds.Rds")
y40_input_cds <- readRDS("40y_input_cds.Rds")

# Merge the multiple cds objects. Include as many input cds objects as required 
# this will merge them using a 'stack' method. 
merged_cds <- combine_cds(list(ga22_input_cds,ga24_input_cds,m1_input_cds,m3_input_cds,m6_input_cds,m10_input_cds,y1_input_cds,y2_input_cds,y6_input_cds,y8_input_cds))

input_cds <- monocle3::detect_genes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

##Quality check and filtering 
# Visualize peak_count_per_cell
hist(Matrix::colSums(exprs(input_cds)))

# Filter cells by peak_count
# Please set an appropriate threshold values according to your data
max_count <-  15000
min_count <- 2000
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) >= min_count]
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) <= max_count]

## Process Cicero-CDS object
# Data preprocessing
set.seed(2017)

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Save Cds object (Optional)
saveRDS(cicero_cds, paste0("Ncicero_cds.Rds"))

## Load reference genome information 
# !!Please make sure that the reference genome information below matches your scATAC-seq reference genome.

# If your scATAC-seq was aligned to the mm10 reference genome, you can read the chromosome length file using the following command.
download.file(url = "https://raw.githubusercontent.com/morris-lab/CellOracle/master/docs/demo_data/mm10_chromosome_length.txt",
              destfile = "./mm10_chromosome_length.txt")
chromosome_length <- read.table("./mm10_chromosome_length.txt")

# For mm9 genome, you can use the following command.
#data("mouse.mm9.genome")
#chromosome_length <- mouse.mm9.genome

# For hg19 genome, you can use the following command.
#data("human.hg19.genome")
#chromosome_length <- mhuman.hg19.genome

## Run Cicero 
# Run the main function
conns <- run_cicero(cicero_cds, chromosome_length) # Takes a few minutes to run

# Save results (Optional)
saveRDS(conns, paste0(output_folder, "/Ncicero_connections.Rds"))

# Check results
head(conns)

## Save results for the next step
# creating csv output folders 
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = paste0(output_folder, "/Nall_peaks.csv"))
write.csv(x = conns, file = paste0(output_folder, "/Ncicero_connections.csv"))

