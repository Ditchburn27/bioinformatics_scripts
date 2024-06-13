### Functions for working with snATAC-data
import matplotlib.pyplot as plt
import seaborn as sns

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
