### Script to create Anndata object from 
### AnnDataSet of Anndata objects in 
### provided input directory.
### Script will add tile matrix, select features,
### remove doublets, and make barcode names unique
### by adding sample_ids to barcodes. 

##################################################
# import libraries
import snapatac2 as snap
import numpy as np
import os
from tqdm import tqdm
import argparse

################################################## 
# CLI inputs
arg_description = '''Script to make Anndata
object from anndata objects in input directory'''

parser = argparse.ArgumentParser(description=arg_description)
parser.add_argument('input_dir', type=str, help='''Path to 
                    input directory containing anndata objects''')
parser.add_argument('output_name', type=str, help='''Name of 
                    output Anndata object file''')
args = parser.parse_args()

input_dir = args.input_dir
output_name = args.output_name

################################################## 
# Read in anndata objects
print('Reading anndata objects...')
adata_files = [f'{input_dir}/{RL}' for RL in os.listdir(input_dir)]
adatas = []
for file in adata_files:
    adata = snap.read(file)
    adatas.append(adata)

################################################## 
# Add cell-by-bin matrix
print('Adding cell-by-bin matrix...')
snap.pp.add_tile_matrix(adatas, bin_size=5000)

################################################## 
# Select features
print('Selecting features...')
snap.pp.select_features(adatas)

################################################## 
# Remove doublets
print('Removing doublets...')
snap.pp.scrublet(adatas)
snap.pp.filter_doublets(adatas)

##################################################
# Creat AnnDataSet object
print('Creating AnnDataSet & Making barcodes unique...')
adataset = snap.AnnDataSet(
    adatas=[RL.split('.h5ad')[0] for RL in adatas],
    filename=f'{input_dir}/{output_name}'
)
# Make the obs_names unique by adding sample names to barcodes
adataset.obs_names = adataset.obs['sample'] + '+' + np.array(adataset.obs_names)

##################################################
# Make AnnDataSet object into anndata
print('''Converting AnnDataSet object to Anndata
      and saving...''')
adata = adataset.to_adata()
adata.write(input_dir+'/'+output_name)
print(f'Anndata object saved here:{input_dir+'/'+output_name}')
