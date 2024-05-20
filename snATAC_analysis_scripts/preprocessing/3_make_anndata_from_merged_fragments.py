import snapatac2 as snap
import argparse

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
adata = snap.pp.import_data(fragment_file,
                            chrom_sizes=snap.genome.GRCh38,
                            sorted_by_barcode=False)

adata.write('braindev_hg38_snATAC_raw.h5ad')