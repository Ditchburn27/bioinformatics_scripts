import snapatac2 as snap
import argparse
import os

######################################################
# CLI inputs
arg_description = '''Script to make h5ad file
from combined cellranger fragments.tsv.gz outputs with 
library_id prefixed barcodes snATAC-seq using snapATAC2'''

parser = argparse.ArgumentParser(description=arg_description)
parser.add_argument('fragment_file', type=str, help='''Path to 
                    combined fragment file
                    ''')
args = parser.parse_args()

fragment_file = args.fragment_file

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
snap.pp.filter_cells(adata, min_counts=2000, min_tsse=None)

print("Adding Tile Matrix...")
snap.pp.add_tile_matrix(adata, bin_size=5000, counting_strategy='fragment')

print("Selecting features...")
blacklist = '/group/ll005/reference/bowtie2_hg38/blacklist/hg38-blacklist.v2.bed'
snap.pp.select_features(adata, n_features=500000, blacklist=blacklist)

adata.write('braindev_hg38_snATAC_tile_matrix.h5ad')
