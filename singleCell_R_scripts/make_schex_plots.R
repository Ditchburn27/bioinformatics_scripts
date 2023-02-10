# run on the stiletto RStudio 

# run in terminal
# change to your username and key
# ssh -N -Y -f -i "/Users/leightonditchburn/Documents/Server/lditchburn_ngs_id.txt" -L 8787:127.0.0.1:8787 lditchburn@130.95.176.224

# open http://127.0.0.1:8787/ in browser

# to install packages run
# install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("scater")
# BiocManager::install("schex")

library(scater)
library(schex)

sce <- readRDS("/Users/leightonditchburn/Documents/Data/Brain_dev/sce_raw.RDS")

# change the number of nbins to get a better picture
sce <- make_hexbin(sce, nbins = 200, 
                   dimension_reduction = "UMAP", use_dims=c(1,2))

gene_id <-"BCL11A"
plot_hexbin_feature(sce, type="counts", feature=gene_id, 
                    action="prop_0", xlab="UMAP1", ylab="UMAP2", 
                    title=paste0("Prop of ", gene_id))

plot_hexbin_feature(sce, type="logcounts", feature=gene_id, 
                    action="mean", xlab="UMAP1", ylab="UMAP2", 
                    title=paste0("Mean of ", gene_id))

