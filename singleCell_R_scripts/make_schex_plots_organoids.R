# run on the stiletto RStudio 

# run in terminal
# change to your username and key
# ssh -N -Y -f -i ~/Documents/Server/lditchburn_ngs_id.txt -L 8787:127.0.0.1:8787 lditchburn@130.95.176.224

# open http://127.0.0.1:8787/ in browser

# to install packages run
# install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("scater")
# BiocManager::install("ggpubr")
# BiocManager::install("ggplot.multistats")

library(scater)
library(ggplot2)
library(ggplot.multistats)
library(ggpubr)

all_sce <- readRDS("/Users/leightonditchburn/Documents/Data/Daniel_organoids/sce_organoids.RDS")


df <- lapply(all_sce, function(x) reducedDim(x, "UMAP"))
df <- do.call(rbind, df)

img <- png::readPNG("/Users/leightonditchburn/Documents/Data/Daniel_organoids/background.png")

gene.id <- "KCNH1"
gene.ens <- rownames(all_sce[[1]])[match(gene.id, rowData(all_sce[[1]])$Symbol)]
if(is.na(gene.ens)){
    
    stop("Gene cannot be found!")
}

df$Gene_log <- do.call(c, sapply(all_sce, function(x) as.vector(logcounts(x[gene.ens,]))))
df$Gene <- do.call(c, sapply(all_sce, function(x) as.vector(counts(x[gene.ens,]))))
df$Gene <- as.numeric(df$Gene>0)

p1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
    background_image(img) +
    stat_summaries_hex(aes(z = Gene_log, fill = stat(mean)),
    funs = c('mean'), bins = 100) + scale_fill_viridis_c() + 
    theme_classic() + scale_x_continuous(limits=c(-6.063561,21.449993)) +
    scale_y_continuous(limits=c(-9.861524,15.003452)) +
    guides(fill=guide_legend(title=gene.id))

p2 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
    background_image(img) +
    stat_summaries_hex(aes(z = Gene, fill = stat(prop)),
        funs = normalize_function_list(list(prop = function(x) 
            sum(x)/length(x))), bins = 100) + scale_fill_viridis_c() + 
    theme_classic() + scale_x_continuous(limits=c(-6.063561,21.449993)) +
    scale_y_continuous(limits=c(-9.861524,15.003452)) +
    guides(fill=guide_legend(title=gene.id))


ggsave(p1, file=paste0("Proportion_", gene_id, "_organoids.png"))
ggsave(p2, file=paste0("Mean_", gene_id, "_organoids.png"))
