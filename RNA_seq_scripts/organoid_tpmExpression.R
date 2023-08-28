# Use this script to have a quick look at the expression of a gene of 
# interest in organoids, using the bulk RNA-seq dataset from my honours project.

library(tidyverse)
library(ggplot2)

# Read in data table with TPM values and hgnc gene symbols
tpm_results <- read.table(
  file = "organoid_bulkRNAseq_tpmExpression_geneSymbol_table.txt", sep = ","
  , header = TRUE)

# Read in metadata to include organoid conditions
metadata <- read.table("sample_number_condition.txt", sep = ",", 
                         header = TRUE)

# Function to plot a histogram of TPM values in all organoid conditions
plot_tpms <- function(gene = "") {
  expression <- subset(tpm_results, tpm_results$hgnc_symbol == gene)
  t.expression <- pivot_longer(expression, cols = 1:24,
                               names_to = 'sample', values_to = 'tpm')
  t.expression <- t.expression[, -1]
  t.expression <- left_join(x = t.expression, y = metadata,"sample")
  
  return(ggplot(t.expression, aes(x = sample, y = tpm, fill = condition)) + 
           geom_col() + ggtitle(paste(
             gene, "expression in organoids", sep = " ")) +
    theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5)) +
    theme(legend.position = "right") + 
    theme(plot.title = element_text(hjust = 0.5)))
}

# Example: plot_tpms(gene = "SMAD3")