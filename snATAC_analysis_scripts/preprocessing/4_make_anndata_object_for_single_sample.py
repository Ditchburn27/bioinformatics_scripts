import snapatac2 as snap
import argparse
import os

######################################################
# CLI inputs
arg_description = '''Script to make h5ad file
from cellranger fragments.tsv.gz outputs for 
individual snATAC-seq samples using snapATAC2.'''

parser = argparse.ArgumentParser(description=arg_description)
parser.add_argument('fragment_file', type=str, help='''Path to 
                    combined fragment file
                    ''')
parser.add_argument('output_file', type=str, help='Anndata file name')
parser.add_argument('bin_size', type=int, help='Genomic region bin size')
args = parser.parse_args()

fragment_file = args.fragment_file
output_file = args.output_file
bin_size = args.bin_size

blacklist_file="/group/ll005/reference/bowtie2_hg38/blacklist/hg38-blacklist.v2.bed"

######################################################
# Import fragment file and make anndata object. 
print("Importing fragment file...")
adata = snap.pp.import_data(fragment_file,
                            chrom_sizes=snap.genome.GRCh38,
                            sorted_by_barcode=False)

######################################################
# Calculate metrics & preprocess
print("Calculating TSSE scores...")
snap.metrics.tsse(adata, snap.genome.GRCh38)

print("Filtering cells...")
snap.pp.filter_cells(adata, min_counts=2000, min_tsse=1.5)

print(f'Adding tile matrix with {bin_size} bp bins')
snap.pp.add_tile_matrix(adata, bin_size=bin_size, counting_strategy='paired-insertion')

print("Selecting features...")
snap.pp.select_features(adata, blacklist=blacklist_file)

print("Detecting and removing doublets...")
snap.pp.scrublet(adata)
snap.pp.filter_doublets(adata)

print("Performing spectral embedding...")
snap.tl.spectral(adata)

print("Performing UMAP...")
snap.tl.umap(adata)

print("Running clustering analysis -> knn & leiden")
snap.pp.knn(adata)
snap.tl.leiden(adata)

print(f'Anndata object saved as {output_file}')
adata.write(output_file)
