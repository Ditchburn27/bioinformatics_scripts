import snapatac2 as snap
import numpy
import os
import argparse
import multiprocess
from tqdm import tqdm

######################################################
# CLI inputs
arg_description = '''Script to make h5ad files
from cellranger fragments.tsv.gz output for 
snATAC-seq using snapATAC2'''

parser = argparse.ArgumentParser(description=arg_description)
parser.add_argument('input_dir', type=str, help='''Path to 
                    parent directory of sample directories
                    that contain fragment files''')
args = parser.parse_args()

input_dir = args.input_dir

######################################################
# Make dir for output files
h5ad_path = f'{input_dir}/h5ad_files'
os.makedirs(h5ad_path, exist_ok=True)

######################################################
# Make lists of input fragment file paths, output file names
# for h5ad files.
fragment_files = [f'{input_dir}/{RL}/outs/fragments.tsv.gz'
             for RL in os.listdir(input_dir) if RL.startswith("RL")]

h5ad_outputs = []
for RL in fragment_files:
    name = RL.split('/')[-3]
    h5ad_outputs.append(f'{h5ad_path}/{name}.h5ad')

######################################################
# Make h5ad files
for fragment_file, h5ad_output in zip(fragment_files,
                                      h5ad_outputs):
    snap.pp.import_data(fragment_file, file=h5ad_output,
                    chrom_sizes=snap.genome.hg38,
                     min_num_fragments=1000,
                     sorted_by_barcode=False)