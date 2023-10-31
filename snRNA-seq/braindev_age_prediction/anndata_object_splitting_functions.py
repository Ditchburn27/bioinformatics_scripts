import scanpy as sc
import numpy as np
from sklearn.model_selection import train_test_split
from geosketch import gs
import pandas as pd


# Function to downsample data to a set number of cells based on the density of the data
def density_split (adata, target_variable = 'numerical_age', split = 0.8):
    # Calculate number of cells that accounts for 80% of the dataset
    split = split
    num_cells = int(split*len(adata.obs))
    # Extract distances from adata and store it in udist
    udist = adata.obsp['distances']
    # Compute the sum of distances for each cell and convert to 1D array
    dist_sum = udist.sum(1).A1
    # Calculate the 95th percentile of distance sums to normalize distances
    pert_th = np.percentile(dist_sum, 95)
    # Calculate density values by subtracting distance sums from the 95th percentile
    density = pert_th - dist_sum
    # Ensure density values are not zero by setting values less than 0 to 0.0001
    density[density < 0] = 0.0001
    # Calculate normalized density values by subtracting each density value from the maximum density
    mxdens = max(density) - density
    # Normalize the density values so they sum up to 1
    adata.obs['norm_density'] = mxdens / sum(mxdens)
    # Randomly select num_cells from adata based on the normalized density values
    bcs = np.random.choice(adata.obs_names, size=num_cells, replace=False, p=adata.obs['norm_density'])
    test_bcs = []
    for barcode in adata.obs.index:
                if barcode not in set(bcs):
                        test_bcs.append(barcode)
    # Subset adata to train and test sets
    train_adata = adata[bcs].copy()
    test_adata = adata[test_bcs].copy()
    X_train = train_adata.X
    X_test = test_adata.X
    y_train = train_adata.obs[target_variable]
    y_test = test_adata.obs[target_variable]
    print(f"X_train shape is: {X_train.shape}")
    print(f"X_test shape is: {X_test.shape}")
    print(f"y_train shape is: {y_train.shape}")
    print(f"y_test shape is: {y_test.shape}")

    return X_train, X_test, y_train, y_test

###################################################################################
###################################################################################

# Function to randomly split anndata object into train and test subsets
# 80% for training and 20% for testing
def random_split(adata, split=0.8, target_variable='numerical_age'):
    split = split
    train_barcodes = np.random.choice(adata.obs.index, replace = False, size= int(split * adata.shape[0]))
    test_barcodes = np.asarray([barcode for barcode in adata.obs.index if barcode not in set(train_barcodes)])
    train_adata = adata[train_barcodes]
    test_adata = adata[test_barcodes]
    # generate train and test data by extracting matrices from adata.X
    X_train = train_adata.X
    X_test = test_adata.X
    y_train = train_adata.obs[target_variable]
    y_test = test_adata.obs[target_variable]
    print(f"X_train shape is: {X_train.shape}")
    print(f"X_test shape is: {X_test.shape}")
    print(f"y_train shape is: {y_train.shape}")
    print(f"y_test shape is: {y_test.shape}")

    return X_train, X_test, y_train, y_test

###################################################################################
###################################################################################
# Function to split anndata object into train and test
# 80% of cells from each sample (adata.obs['batch']) will be allocated for training
# 20% of cells from each sample will be allocated for testing
def batch_split (adata, target_variable ='numerical_age', split = 0.8):
    # Identiies the unique batches or samples in adata
    unique_batches = adata.obs['batch'].unique()
    # Calculate the number of cells to be included from each batch for 80% training
    split = split
    n_train_cells_per_batch = {
        batch: int(split * len(adata.obs[adata.obs['batch'] == batch]))
        for batch in unique_batches}
    # Calculate the number of cells to be included from each batch for testing
    n_test_cells_per_batch = {
        batch: len(adata.obs[adata.obs['batch'] == batch]) - n_train_cells_per_batch[batch]
        for batch in unique_batches}
    # Initialize empty arrays to hold train and test indicies
    train_indices = np.array([], dtype=int)
    test_indices = np.array([], dtype=int)
    # Loop through unique batches and randomly sample cells 
    # for both train and test sets
    for batch in unique_batches:
         # Get the indices of cells in the current batch
         batch_indices = np.where(adata.obs['batch'] == batch)[0]

         # Split the batch indices into train and test indices
         train_batch_indices, test_batch_indices = train_test_split(batch_indices, 
                                                                   train_size=n_train_cells_per_batch[batch], 
                                                                   test_size=n_test_cells_per_batch[batch], 
                                                                   random_state=42)

         # Add the batch-specific train and test indices to the overall arrays
         train_indices = np.concatenate((train_indices, train_batch_indices))
         test_indices = np.concatenate((test_indices, test_batch_indices))
    # Subset with the train and test indicies to generate the datasets     
    X_train = adata.X[train_indices]
    X_test = adata.X[test_indices]
    y_train = adata.obs[target_variable][train_indices] 
    y_test = adata.obs[target_variable][test_indices]
    print(f"X_train shape is: {X_train.shape}")
    print(f"X_test shape is: {X_test.shape}")
    print(f"y_train shape is: {y_train.shape}")
    print(f"y_test shape is: {y_test.shape}")
    return X_train, X_test, y_train, y_test

###################################################################################
###################################################################################

# Function to split anndata into train and test data using geosketching
# Need to start with an anndata object that has already had PCA run on it
def geosplit (adata, split = 0.8, target_variable = 'numerical_age'):
    split = split
    # Geosketching to get 80% of the dataset
    N = int(split*len(adata.obs))
    sketch_index = gs(adata.obsm['X_pca'], N, replace = False)
    # Create a list of the train indicies (80% geosketch)
    # and a list of the test indicies (remaining 20%)
    df = adata.obsm['X_pca']
    df = pd.DataFrame.from_records(df)
    df_sketch = df.iloc[sketch_index, :]
    train_indices = df_sketch.index.to_list()
    all_indices = df.index.to_list()
    test_indices = list(set(all_indices).difference(train_indices))
    # Subset with the train and test indicies to generate the datasets
    X_train = adata.X[train_indices]
    X_test = adata.X[test_indices]
    y_train = adata.obs[target_variable][train_indices] 
    y_test = adata.obs[target_variable][test_indices]
    print(f"X_train shape is: {X_train.shape}")
    print(f"X_test shape is: {X_test.shape}")
    print(f"y_train shape is: {y_train.shape}")
    print(f"y_test shape is: {y_test.shape}")
    return X_train, X_test, y_train, y_test




