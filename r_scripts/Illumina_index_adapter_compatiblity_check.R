#!/usr/local/bin/Rscript

# version 1.1 with command line options

suppressMessages(library(DNABarcodes))
suppressMessages(library(Biostrings))
suppressMessages(library(optparse))
suppressMessages(library(superheat))


### USER PARAMETERS ###

option_list = list(
  make_option(c("-i", "--index_set1"), type="character", default=NULL, 
              help="CSV file name of index set1. First column index names, second column index sequences. Separated by ','", metavar=" <TruSeq_i7_set1.csv>"),
  make_option(c("-j", "--index_set2"), type="character", default=NULL, 
              help="CSV file name of index set2. First column index names, second column index sequences. Separated by ','", metavar=" <Nextera_i7_set1.csv>"),
  make_option(c("-p", "--project"), type="character", default="Barcode Collision Overview", 
              help="Name of the project. Can include spaces if the whole name is wrapped by \". E.g.  \"Barcode Collision Nextera vs TruSeq\" ", metavar="character"),
  make_option(c("-n", "--number_bases_in_common"), type="integer", default=5, 
              help="Maximum number of bases allowed to be in common between provided indexes. Makes it easier to read and find compatible index combos.", metavar="character"),
  make_option(c("-d", "--use_hamming_dist"), type="logical", default=FALSE,
              help="Set to TRUE if you want to see hamming distance (number of bases different) instead of number of bases in common (default)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Check options

if (any(is.null(opt$index_set1), is.null(opt$index_set2))){
  print_help(opt_parser)
  stop("Two index csv files need to be provided. It can be the same index file for both index sets if you want to compare all indexes against each other.n", call.=FALSE)
}


use_hamming_distance=opt$use_hamming_dist # Hamming distance or number bases in common
project=opt$project # Project name

file_name1 <- opt$index_set1
file_name2 <- opt$index_set2
bases_in_common <- opt$number_bases_in_common
input_file1 <- read.csv(file_name1, header=TRUE)
input_file2  <- read.csv(file_name2, header=TRUE)
use_reverse_complements_list1=FALSE
use_reverse_complements_list2=FALSE

length_truncated_barcodes=8 #NULL if you don't want to truncate your barcodes, desired length otherwise


########

list1 <- input_file1[,2]
list2 <- input_file2[,2]

#truncate barcodes
if (!is.null(length_truncated_barcodes)){
	list1 <- substring(list1, 1, length_truncated_barcodes)
	list2 <- substring(list2, 1, length_truncated_barcodes)}

#reverse complement barcodes
if (isTRUE(use_reverse_complements_list1)){
  RC_list1 <- as.character(reverseComplement(DNAStringSet(list1)))
  list1 <- RC_list1}
if (isTRUE(use_reverse_complements_list2)){
  RC_list2 <- as.character(reverseComplement(DNAStringSet(list2)))
  list2 <- RC_list2}

barcode_list <- c(list1, list2)

#calculate distances
dist_matrix <- barcode.set.distances(barcode_list, metric="hamming")
if (isFALSE(use_hamming_distance)){
	dist_matrix <- nchar(barcode_list[1]) - dist_matrix }

new_matrix <- as.matrix(dist_matrix)
rownames(new_matrix) <- paste(c(input_file1[,1], input_file2[,1]), c(list1, list2), sep="-")
colnames(new_matrix) <- paste(c(input_file1[,1], input_file2[,1]), c(input_file1[,2], input_file2[,2]), sep="-")
sub_matrix <- new_matrix[c(1:length(list1)), -c(0:length(list1))]

#create heatmap
col <- colorRampPalette(c("firebrick1","coral1","grey","seagreen3","seagreen4","seagreen"))(6)
if (isFALSE(use_hamming_distance)){
	col <- colorRampPalette(c("seagreen","seagreen4","seagreen3","grey","coral1","firebrick1"))(6) }

type="hamming distance"
my_breaks=c(0,0.1,0.3,0.5,1)
if (isFALSE(use_hamming_distance)){
type="similarity"
my_breaks=c(0,0.5,0.7,0.9,1)
}

png(paste(project, "barcode_heatmap.png", sep="_"), height = 1900, width = 1800)
map <- superheat(sub_matrix, scale=F, pretty.order.rows = F ,  pretty.order.cols = F, 
bottom.label.text.angle = 90, legend = FALSE,
left.label.text.col = "chartreuse4",
bottom.label.text.col = "cornflowerblue",
left.label.text.size = 6, bottom.label.text.size = 6,
grid.hline.col = "white", grid.vline.col = "white", grid.hline.size = 1.3, grid.vline.size = 1.3,
title = paste("Project ", project, " barcode ", type, sep=""),  title.size = 8, #title.alignment = "center",
row.title = paste(file_name1,colnames(input_file1)[1], sep=" - "), row.title.size = 8, column.title = paste(file_name2, colnames(input_file2)[1], sep=" - "), column.title.size = 8,
heat.pal=col, heat.pal.values = my_breaks,  X.text = sub_matrix)
dev.off()


# Subselect the samples with less than 'n' or less bases in common
sub_matrix_filtered <- apply(sub_matrix, 2, function(x) sum(x > bases_in_common))
sub_matrix_filtered <- sub_matrix_filtered == 0

png(paste(project, "barcode_heatmap_filtered.png", sep="_"), height = 1900, width = 1800)
map <- superheat(sub_matrix[, sub_matrix_filtered], scale=F, pretty.order.rows = F ,  pretty.order.cols = F, 
                 bottom.label.text.angle = 90, legend = FALSE,
                 left.label.text.col = "chartreuse4",
                 bottom.label.text.col = "cornflowerblue",
                 left.label.text.size = 6, 
                 bottom.label.text.size = 6,
                 grid.hline.col = "white", grid.vline.col = "white", grid.hline.size = 1.3, grid.vline.size = 1.3,
                 title = paste("Project ", project, " barcode ", type, sep=""),  title.size = 8, #title.alignment = "center",
                 row.title = paste(file_name1,colnames(input_file1)[1], sep=" - "), 
                 row.title.size = 8, 
                 column.title = paste(file_name2, colnames(input_file2)[1], sep=" - ")[sub_matrix_filtered], column.title.size = 5,
                 heat.pal=col, heat.pal.values = my_breaks,  X.text = sub_matrix[, sub_matrix_filtered])
dev.off()



