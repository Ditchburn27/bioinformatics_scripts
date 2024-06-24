### Functions for working with snATAC-data
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import snapatac2 as snap
from pybedtools import BedTool
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

#####################################
## Plot cell counts of specified obs for 
## different tsse thresholds. 
def plot_tsse_countplots(adata, groupby_observation, log_scale=False):
    """
    Plots cell counts of specified obs for different tsse thresholds.
    
    Parameters:
    - adata: Anndata object containing the observations data.
    - groupby_observation: The key in adata.obs to group by (e.g., 'cell_type').
    - log_scale: Boolean, if True, sets y-axis to log scale.
    """
    # Create figure and 2x2 subplots with smaller figsize
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharey=True)
    # Flatten the axes array for easy iteration
    axes = axes.flatten()
    # Define the conditions and titles for subplots
    conditions = [
        (adata.obs['tsse'] < 5),
        (adata.obs['tsse'] >= 5) & (adata.obs['tsse'] < 10),
        (adata.obs['tsse'] >= 10) & (adata.obs['tsse'] < 20),
        (adata.obs['tsse'] >= 20)
    ]
    titles = [
        'tsse < 5',
        '5 <= tsse < 10',
        '10 <= tsse < 20',
        'tsse >= 20'
    ]
    # Plot histograms for each condition
    for ax, condition, title in zip(axes, conditions, titles):
        subset = adata.obs[condition]
        sns.countplot(data=subset, x=groupby_observation, ax=ax, palette="Set2", color="lightblue")
        ax.set_title(title)
        ax.set_xlabel(groupby_observation.capitalize())
        ax.set_ylabel('Number of Cells')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)  # Rotate x-axis labels for readability
        if log_scale:
            ax.set_yscale('log')  # Set y-axis to log scale
    # Adjust layout for better fit
    plt.tight_layout()
    plt.show()

#####################################
## Function for mapping sub_cluster cell types to major_cluster
# Define the function to map cell types to major clusters
def map_to_major_cluster(cell_type):
    if 'dev' in cell_type:
        if cell_type.startswith('L2/3_') or cell_type.startswith('L2_') or cell_type.startswith('L3_'):
            return 'L2/3_dev'
        elif cell_type.startswith('L4_'):
            return 'L4_dev'
        elif cell_type.startswith('L5/6_'):
            return 'L5/6_dev'
        elif cell_type.startswith('Astro_'):
            return 'Astro_dev'
        elif cell_type.startswith('PV_'):
            return 'PV_dev'
        elif cell_type.startswith('SST_'):
            return 'SST_dev'
        elif cell_type.startswith('VIP_'):
            return 'VIP_dev'
    elif cell_type.startswith('L2/3_') or cell_type.startswith('L2_') or cell_type.startswith('L3_'):
        return 'L2/3'
    elif cell_type.startswith('L4_'):
        return 'L4'
    elif cell_type.startswith('L5/6_'):
        return 'L5/6'
    elif cell_type.startswith('Astro_'):
        return 'Astro'
    elif cell_type.startswith('Oligo'):
        return 'Oligo'
    elif cell_type.startswith('PV_'):
        return 'PV'
    elif cell_type.startswith('SST_'):
        return 'SST'
    elif cell_type.startswith('VIP_'):
        return 'VIP'
    elif cell_type.startswith('Vas_'):
        return 'Vas'
    elif cell_type.startswith('MGE_'):
        return 'MGE'
    elif cell_type.startswith('LAMP5_'):
        return 'LAMP5'
    elif cell_type.startswith('Micro'):
        return 'Micro'
    elif cell_type.startswith('CCK_'):
        return 'CCK'
    elif cell_type.startswith('ID2_'):
        return 'ID2'
    else:
        return cell_type
# example usage:
#adata.obs['major_cluster'] = adata.obs['cell_type'].apply(map_to_major_cluster)

#####################################
## Function to map major_cluster into broader groups
# Define the function to map major clusters to cluster groups
def map_to_cluster_group(major_cluster):
    if major_cluster in ['L2/3_dev', 'L4_dev', 'L5/6_dev']:
        return 'PN_dev'
    elif major_cluster in ['L2/3', 'L4', 'L5/6']:
        return 'PN'
    elif major_cluster in ['PV_dev', 'SST_dev', 'VIP_dev']:
        return 'IN_dev'
    elif major_cluster in ['PV', 'SST', 'VIP', 'LAMP5']:
        return 'IN'
    else:
        return major_cluster
# example usage:
#data.obs['cluster_group'] = adata.obs['major_cluster'].apply(map_to_cluster_group)

#####################################
## Function to plot a grid of subplots for counting cell numbers
## for each label in the provided count key
def plot_grouped_counts(adata, group_key, count_key, num_cols=5):
    """
    Plots subplots of count plots for each unique group in adata.obs.

    Parameters:
    - adata: Anndata object containing the observations data.
    - group_key: The key in adata.obs to group by (e.g., 'batch').
    - count_key: The key in adata.obs to count (e.g., 'major_cluster').
    - num_cols: Number of columns in the subplot grid.
    """
    unique_groups = adata.obs[group_key].unique()
    num_subplots = len(unique_groups)
    num_rows = -(-num_subplots // num_cols)  # Ceiling division to ensure at least num_subplots number of rows
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(5*num_cols, 5*num_rows), sharey=True)
    # Flatten the axes array if it's not already 1-dimensional
    axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    # Define a color palette for each unique count_key
    palette = sns.color_palette("husl", len(adata.obs[count_key].unique()))
    # Iterate over each unique group
    for i, group in enumerate(unique_groups):
        # Filter data for current group
        group_data = adata.obs[adata.obs[group_key] == group]    
        # Count plot for each count_key label
        sns.countplot(x=count_key, data=group_data, ax=axes[i], palette=palette)
        axes[i].set_title(f'{group_key.capitalize()} {group}')      
        # Rotate x-labels for better readability
        axes[i].tick_params(axis='x', rotation=45)        
        # Set labels for better readability
        axes[i].set_xlabel(count_key.capitalize())
        axes[i].set_ylabel('Count')
    # Remove empty subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout()
    plt.show()

#####################################
## Downsampling function - 
## adapted from Chuck's function for snRNA-seq
# standardize library sizes by removing nucs with low libs and downsampling the rest
def standardize_libs(adata, tar_lib_sz=2500, seed=123):
    """
    Standardizes library sizes in the Anndata object by removing low-count cells and downsampling.

    Parameters:
    - adata: Anndata object containing the observations data.
    - tar_lib_sz: Target library size for standardizing cell counts. Default is 2500.
    - seed: Random seed for reproducibility of the downsampling. Default is 123.

    This function performs the following steps:
    1. Identifies and prints the number of cells (nuclei) with total counts less than or equal to the target library size,
       grouped by 'library_id'.
    2. Removes cells with total counts less than the target library size.
    3. Downsamples the total counts of each remaining cell to the target library size.
    4. Resets the data type of the counts matrix to float to address a specific bug in the scanpy version used.

    Modifies:
    - The function modifies the adata object in-place by potentially removing low-count cells and adjusting the counts 
      of the remaining cells to standardize their library sizes.
    """
    # make sure counts for each nuc are current
    #sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)

    # print number of nucs to be removed
    print("Number of low lib nucs removed per run:")
    low_mask = (adata.obs['n_fragment'] <= tar_lib_sz)
    print(pd.value_counts(adata.obs['library_id'].values[low_mask]).to_string())
    
    # remove nucs with low library sizes
    #### comment out line below to keep low UMI count nuclei ####
    sc.pp.filter_cells(adata, min_counts=tar_lib_sz, copy=False)
    
    # downsample nuc total counts over target library size
    sc.pp.downsample_counts(adata, counts_per_cell=tar_lib_sz, copy=False, replace=False, random_state=seed)
    
    # due to scanpy bug in this version have to reset dtype
    adata.X = adata.X.astype(float)
    return

#####################################
## Function for spectral embedding to umap
def spec_to_umap(adata, n_comps=30, distance_metric='cosine', sample_size=1.0, n_neighbors=50, seed=123, 
                 weighted_by_sd=True):
    """
    Performs spectral embedding followed by UMAP for dimensionality reduction on the given Anndata object.

    Parameters:
    - adata: Anndata object containing the observations data.
    - n_comps: Number of principal components to use for spectral embedding. Default is 30.
    - distance_metric: Distance metric to use for spectral embedding. Default is 'cosine'.
    - sample_size: Proportion of the data to sample if using a distance metric other than 'cosine'. 
                   A value between 0.0 and 1.0. Default is 1.0.
    - n_neighbors: Number of neighbors to use for the k-nearest neighbors graph construction in UMAP. Default is 50.
    - seed: Random seed for reproducibility of the spectral embedding and UMAP. Default is 123.
    - weighted_by_sd: If True, weights by the standard deviation when performing spectral embedding. Default is True.

    Modifies:
    - The function modifies the adata object in-place, adding the UMAP coordinates and cluster assignments.
    """
    if distance_metric=='cosine':
        snap.tl.spectral(adata, random_state=seed, n_comps=n_comps, weighted_by_sd=weighted_by_sd)
        snap.pp.knn(adata, n_neighbors=n_neighbors, random_state=seed)
        snap.tl.leiden(adata)
        snap.tl.umap(adata, random_state=seed)
    else:
        snap.tl.spectral(adata, distance_metric=distance_metric,
                         sample_size=sample_size, random_state=seed,
                         n_comps=n_comps, weighted_by_sd=weighted_by_sd)
        snap.pp.knn(adata, n_neighbors=n_neighbors, random_state=seed)
        snap.tl.leiden(adata)
        snap.tl.umap(adata, random_state=seed)
    return

#####################################
# Function for identifying genomic bins
# that overlap with features in provided
# bed file e.g. promoter, enhancers
def region_overlaps_bed(args):
    region, bed_file_bed = args
    chrom, start_end = region.split(':')
    start, end = map(int, start_end.split('-'))
    region_bed = BedTool(f"{chrom} {start} {end}", from_string=True)
    return bool(region_bed.intersect(bed_file_bed, u=True))

def select_bin_region_features(adata, region_dict, region_column='regions'):
    """
    Adds columns to adata.var indicating if each region overlaps with any regions in the provided BED files.
    
    Parameters:
    adata (AnnData): The AnnData object containing the data.
    region_dict (dict): A dictionary where keys are the output column names and values are the paths to BED files.
    region_column (str): The column in adata.var that contains the regions.
    """
    # Convert adata.var regions to a DataFrame
    adata_df = pd.DataFrame({
        'chrom': [region.split(':')[0] for region in adata.var[region_column]],
        'start': [int(region.split(':')[1].split('-')[0]) for region in adata.var[region_column]],
        'end': [int(region.split(':')[1].split('-')[1]) for region in adata.var[region_column]],
        'region': adata.var[region_column]
    })

    # Loop through each item in the dictionary
    for output_column, bed_file_path in region_dict.items():
        # Read the BED file using pybedtools
        bed_file_bed = BedTool(bed_file_path)

        # Prepare arguments for multiprocessing
        args = [(region, bed_file_bed) for region in adata.var[region_column]]

        # Use multiprocessing to process regions in parallel with progress bar
        with Pool(cpu_count()) as pool:
            results = list(tqdm(pool.imap(region_overlaps_bed, args), total=len(args), desc=f"Processing {output_column}"))

        # Add the results to the new column in adata.var
        adata.var[output_column] = results

#####################################
# Function to print counts of regions
# that are identified to overlap with
# features from bed file. 
def count_region_features(adata, region_feature_label='selected'):
    true_count = adata.var[region_feature_label].sum()
    print(f"Number of genomic bins that are {region_feature_label}: {true_count}")