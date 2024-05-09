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
output_path = f'{input_dir}/{output_name}'

################################################## 
# Read in anndata objects
def read_adata(file):
    return snap.read(file)

def read_adatas(adata_files):
    adatas = []
    for file in adata_files:
        adata = read_adata(file)
        adatas.append(adata)
    return adatas

if __name__ == '__main__':
    print('Reading anndata objects...')
    adata_files = [f'{input_dir}/{RL}' for RL in os.listdir(input_dir)]
    adatas = read_adatas(adata_files)
        
    ################################################## 
    # Add cell-by-bin matrix
    print('Adding cell-by-bin matrix...')
    snap.pp.add_tile_matrix(adatas, bin_size=5000, n_jobs=16)

    ##################################################
    # Creat AnnDataSet object
    print('Creating AnnDataSet & Making barcodes unique...')
    adataset = snap.AnnDataSet(
        adatas=[(f.filename.split('/')[-1].split('.h5ad')[0], f) for f in adatas],
        filename=f'{input_dir}/{output_name}'
    )
    # Make the obs_names unique by adding sample names to barcodes
    adataset.obs_names = adataset.obs['sample'] + '+' + np.array(adataset.obs_names)

    ##################################################
    # Make AnnDataSet object into anndata
    print('Converting AnnDataSet object to Anndata')
    adata = adataset.to_adata()

    ################################################## 
    # Select features
    print('Selecting features...')
    snap.pp.select_features(adata, n_jobs=16)

    ################################################## 
    # Remove doublets
    print('Removing doublets...')
    snap.pp.scrublet(adata, n_jobs=16)
    snap.pp.filter_doublets(adata, n_jobs=16)
    ################################################## 
    # Remove doublets
    adata.write(output_path)
    print(f'Anndata object saved here:{output_path}')
