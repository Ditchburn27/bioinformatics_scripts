## Adding basepair info to DESeq Object
## + TPM function 
# Need to read in a RNA-seq counts matrix and sample info after 
# 'library(tidyverse)' step and instantiate a DESeq object 
library(GenomicFeatures)
library(biomaRt)

txdb <- makeTxDbFromGFF("gencodeV19_110_103_hChr2_ERCC.gtf.gz", format="gtf")
# Collect exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, 
# calculate their lengths (widths) and sum
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))

library(tidyverse)

# Convert ENS gene ID to human gene name
# connect to biomart database 
ensembl <- useEnsembl(biomart = 'genes', 
                      version = 75,
                      host = "https://feb2014.archive.ensembl.org")
# choose dataset 
mart <- useDataset(dataset = 'hsapiens_gene_ensembl', mart = ensembl) 

ens_genes <-names(exonic.gene.sizes)
ens_genes <- stringr::str_remove_all(ens_genes, pattern = "\\..*")
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=ens_genes,mart= mart)
# Rename human ENSM gene list
exonic.gene.sizes.df <- data.frame(exonic.gene.sizes) %>% 
  tibble::rownames_to_column(var ="ensembl_gene_id")
exonic.gene.sizes.df$ensembl_gene_id <- stringr::str_remove_all(
  exonic.gene.sizes.df$ensembl_gene_id, pattern = "\\..*")

# Intersect Ensemble gene ids with the gene symbol and total exon sizes
exonic.gene.sizes.symbol.df <- dplyr::left_join(exonic.gene.sizes.df, 
                                                G_list, by = "ensembl_gene_id" )

# Subset the DESeq2 object with the gene symbols I have gene length info
dds.genelength<- dds[rownames(dds) %in% 
                       exonic.gene.sizes.symbol.df$ensembl_gene_id,]
# Only keep unique gene symbol entries
exonic.gene.sizes.symbol.unique.df <- 
  exonic.gene.sizes.symbol.df[! duplicated(
    exonic.gene.sizes.symbol.df$ensembl_gene_id), ]
# Add the exon length to the metadata column `basepairs`
mcols(dds.genelength)$basepairs <- 
  exonic.gene.sizes.symbol.unique.df[
    exonic.gene.sizes.symbol.unique.df$ensembl_gene_id %in% 
      rownames(dds.genelength), "exonic.gene.sizes"]

# TPM function 
tpm <- function(counts, gene_bp_sizes) {
  # Calculate kilobase pair sizes
  gene_kbp_sizes <- gene_bp_sizes / 1000
  # Divide all counts by the rpk size of the gene
  counts_rpk <- apply(counts, 2, function(x) x / gene_kbp_sizes)
  # Scaling factor per million
  scaling <- colSums(counts_rpk) / 1E6
  # Divide each row by the scaling factor
  counts_tpm <- counts_rpk / scaling
  return(counts_tpm)
}

tpm_results <- tpm(counts(dds.genelength), mcols(dds.genelength)$basepairs )
tpm_results <- as.data.frame(tpm_results)
tpm_results <- tibble::rownames_to_column(tpm_results, 'row_names')
names(tpm_results)[1] <- 'ensembl_gene_id'
