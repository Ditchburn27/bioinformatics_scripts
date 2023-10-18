import scanpy as sc
import numpy as np



# Function to downsample data to a set number of cells based on the density of the data
def density_ds (adata, num_cells = 10000):
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
    # Subset adata to include only the selected cells
    ds_adata = adata[bcs]
    return ds_adata

# Function to randomly split anndata object into train and test subsets
# 80% for training and 20% for testing
def split_adata(adata, split=0.8):
    split = split
    train_barcodes = np.random.choice(adata.obs.index, replace = False, size= int(split * adata.shape[0]))
    test_barcodes = np.asarray([barcode for barcode in adata.obs.index if barcode not in set(train_barcodes)])
    train_adata = adata[train_barcodes]
    test_adata = adata[test_barcodes]

    print (len(train_adata), len(test_adata))
    return train_adata, test_adata

# Function that returns the numerical_age in numpy arrays for the corresponding train and test data
def num_age_split(train_adata, test_adata):
    y_train = train_adata.obs['numerical_age'].to_numpy().reshape(-1,1)
    y_test = test_adata.obs['numerical_age'].to_numpy().reshape(-1,1)
    print(y_train.shape, y_test.shape)
    return y_train, y_test
