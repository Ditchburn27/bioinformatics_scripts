## Convert Seurat object to anndata 
# Need to have Seurat and SeuratDisk packages installed 

# Load libraries
library(Seurat)
library(SeuratDisk)

# Read in seurat object
L5_6 <- readRDS("/home/lditchburn/working_data_02/snRNA/L5_6.Rds")

# First save object as Seurat hdf5 file 
SaveH5Seurat(allPN, filename = "allPN.h5Seurat")

# Convert Seurat h5 to anndata 
Convert("allPN.h5Seurat", dest = "h5ad")
