import snapatac2 as snap
import argparse
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import os
import sys
import concurrent.futures

# Import custom functions
sys.path.append('/group/ll005/lditchburn/bioinformatics_scripts/snATAC_analysis_scripts')
from snatac_utils import region_overlaps_bed, select_bin_region_features

######################################################
# Define functions to perform operations on each AnnData object
def select_features(ad):
    snap.pp.select_features(ad, blacklist=blacklist_file)
    return ad

def run_scrublet(ad):
    snap.pp.scrublet(ad)
    if 'doublet_probability' in ad.obs:
        ad.obs['predicted_doublet'] = ad.obs['doublet_probability'] > 0.5
    else:
        raise ValueError("The 'doublet_probability' column does not exist in one of the AnnData objects.")
    return ad

if __name__ == "__main__":
    ######################################################
    # CLI inputs
    arg_description = '''Script to make h5ad file
    from combined cellranger fragments.tsv.gz outputs with 
    library_id prefixed barcodes snATAC-seq using snapATAC2'''

    parser = argparse.ArgumentParser(description=arg_description)
    parser.add_argument('fragment_file', type=str, help='''Path to 
                        combined fragment file
                        ''')
    parser.add_argument('output_file', type=str, help='Anndata file name')
    parser.add_argument('bin_size', type=int, help='Genomic region bin size')
    parser.add_argument('metadata_file', type=str, help='Path to metadata data csv file')
    parser.add_argument('region_bed_files_dir', type=str, help='''Path to directory
                        containing genomic region/ feature bed files (e.g. promoters)''')
    args = parser.parse_args()

    fragment_file = args.fragment_file
    output_file = args.output_file
    bin_size = args.bin_size
    metadata_file = args.metadata_file
    bed_file_dir = args.region_bed_files_dir

    blacklist_file = "/group/ll005/reference/bowtie2_hg38/blacklist/hg38-blacklist.v2.bed"

    # Create bed file variables
    files = os.listdir(bed_file_dir)
    bed_files = {}
    for file in files:
        variable_name = os.path.splitext(file)[0]
        variable_value = os.path.join(bed_file_dir, file)  # Include the full path
        bed_files[variable_name] = variable_value
        exec(f"{variable_name} = '{variable_value}'")
    
    print("Provided bed files are:")
    for variable_name, variable_value in bed_files.items():
        print(f"{variable_name} = {variable_value}")

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
    snap.pp.add_tile_matrix(adata, bin_size=bin_size, counting_strategy='paired-insertion')

    print("Calculating log10 fragments per cell")
    # Function for converting 'n_fragment' to int and calculating log10 'n_fragment' for plotting.
    def n_frag_int(adata, obs_col='n_fragment'):
        # Get obs column series
        n_fragment = adata.obs[obs_col]
        # Ensure it is a pandas series
        if not isinstance(n_fragment, pd.Series):
            n_fragment = pd.Series(n_fragment)
        # Convert data type to integer
        n_fragment = n_fragment.astype(int)
        # Get log10 of n_fragment
        log10_n_fragment = np.log10(n_fragment)
        # Write the modified column back to the AnnData object
        adata.obs['n_fragment'] = n_fragment
        adata.obs['log10_n_fragment'] = log10_n_fragment
        return
    n_frag_int(adata)

    ##########################################################
    # Move adata to memory to add metadata and remove doublets
    print('Adding metadata to anndata object...')
    adata.to_memory()
    adata.write('tmp.h5ad')

    adata = sc.read_h5ad('tmp.h5ad')

    metadata = pd.read_csv(metadata_file, header=[0])
    # Extract library ids from barcode names
    barcodes = adata.obs_names
    adata.obs['barcodes'] = barcodes
    library_ids = [barcode.split('+')[0] for barcode in adata.obs_names]
    # Add new .obs column for library id
    adata.obs['library_id'] = library_ids
    adata.obs['library_id'] = adata.obs['library_id'].astype('category')
    # Add metadata to adata.obs
    adata.obs = adata.obs.set_index('library_id').join(metadata.set_index('library_id')).reset_index()
    # Add batch order
    adata.uns['batch_order'] = ['RL2366', 'RL2207', 'RL2367', 'RL1914', 'RL2208', 'RL2371', 'RL2209', 'RL1784', 'RL2210', 
                                'RL2364', 'RL1994', 'RL2368', 'RL2372', 'RL1785', 'RL2085', 'RL2369', 'RL2373']

    # Set column to categorical 
    adata.obs['stage_id'] = adata.obs['stage_id'].astype('category')

    # Before subsetting save adata.uns['reference_sequences'] to variable
    reference_sequences = adata.uns['reference_sequences']

    ## Identify doublets on a per sample basis
    print("Splitting anndata object samplewise...")
    library_ids = adata.uns['batch_order']
    # Subset anndata for each sample
    adatas = [adata[adata.obs['library_id'] == lib_id].copy() for lib_id in library_ids]
   # Select features and run scrublet for each sample in parallel
    print("Selecting features and identifying doublets...")
    with concurrent.futures.ThreadPoolExecutor() as executor:
        adatas = list(executor.map(select_features, adatas))
        adatas = list(executor.map(run_scrublet, adatas))
    # Concatenate anndata objects back into 1 object
    adata = anndata.concat(adatas, axis=0)

    ## Add adata.uns after concatenating
    # Set order and colour of stages
    adata.uns['stage_order'] = ['Fetal', 'Neonatal', 'Infancy', 'Childhood', 'Adolescence', 'Adult']

    # Add batch order
    adata.uns['batch_order'] = ['RL2366', 'RL2207', 'RL2367', 'RL1914', 'RL2208', 'RL2371', 'RL2209', 'RL1784', 'RL2210', 
                                'RL2364', 'RL1994', 'RL2368', 'RL2372', 'RL1785', 'RL2085', 'RL2369', 'RL2373']

    # Add age order 
    adata.uns['age_order'] = ['ga22', 'ga24', '28d', '86d', '179d', '301d', '627d', '758d', '4yr', 
                              '6yr', '8yr', '10yr', '14yr', '16yr', '20yr', '25yr', '40yr']
    # Add colours and colour dict to be used for plotting
    stage_ordered_colors = ["#512568", "#443682", "#3D6B93", "#20988C", "#98CA43", "#F9E51B"]
    adata.uns['stage_colors_dict'] = dict(zip(adata.uns['stage_order'], stage_ordered_colors))
    adata.uns['stage_id_colors'] = [adata.uns['stage_colors_dict'][ii] for ii in adata.obs['stage_id'].cat.categories]

    # Add reference_sequences
    adata.uns['reference_sequences'] = reference_sequences

    #########################################################
    # Feature selection, dimensionality reduction, clustering
    # Repeat feature selection for concatenated data
    print("Selecting features...")
    snap.pp.select_features(adata, blacklist=blacklist_file)

    print("Identifying bins containing genomic feature regions in provided bed files...")
    adata.var['regions'] = adata.var_names
    region_dict = {key + '_regions': value for key, value in bed_files.items()}
    select_bin_region_features(adata, region_dict)

    print("Calculating fraction of reads in genomic regions of provided bed files...")
    # Create frip_dict from bed_file dict
    frip_dict = {key + '_frac': value for key, value in bed_files.items()}
    snap.metrics.frip(adata, frip_dict)

    print("Performing spectral embedding...")
    snap.tl.spectral(adata)

    print("Performing UMAP...")
    snap.tl.umap(adata)

    print("Running clustering analysis -> knn & leiden")
    snap.pp.knn(adata)
    snap.tl.leiden(adata)

    #########################################################
    # Save anndata object
    print(f'Anndata object saved as {output_file}.h5ad')
    adata.write(f'{output_file}.h5ad')
