######### WORK IN PROGRESS DOES NOT CURRENTLY WORK #########
### Script to tune ElasticNet Hyperparameters on 
### different methods of splitting brain-dev dataset
## Import libraries
import scanpy as sc
import anndata as ad
import pandas as pd
import cupy as cp
from cupyx.scipy.sparse import csr_matrix as csr_gpu
from cuml import ElasticNet
from cuml.metrics.regression import r2_score
from cuml.metrics.regression import mean_absolute_error
from sklearn.model_selection import RandomizedSearchCV
import anndata_object_splitting_functions as osf 
import cudf
import rapids_singlecell as rsc

sc.logging.print_header()

## read in anndata
adata = sc.read('2023-05-31_braindev_hg38_post-ds1000cts-clustered.h5ad')
rsc.utils.anndata_to_GPU(adata)

# subset adata to principal neurons only
PNs_annots = ['L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev']
PN_adata = adata[adata.obs['major_clust'].isin(PNs_annots)].copy()

## set up parameter list
enet_params = [{'alpha':0.1}, {'alpha':0.2}, {'alpha':0.3}, {'alpha':0.4}, 
                 {'alpha':0.5}, {'alpha':0.6}, {'alpha':0.7}, {'alpha':0.8}, {'alpha':0.9}, {'alpha':1},
                 {'l1_ratio':0.1}, {'l1_ratio':0.2}, {'l1_ratio':0.3}, {'l1_ratio':0.4}, 
                 {'l1_ratio':0.5}, {'l1_ratio':0.6}, {'l1_ratio':0.7}, {'l1_ratio':0.8}, {'l1_ratio':0.9}, {'l1_ratio':1},
                 ]

# create dataframe for storing results
results_df = pd.DataFrame(columns = ['data_split_method',
                                     'params', 'r2_score', 'MAE'])

## Create 2D list of names of data splitting methods,
# splitting functions, and the model to be trained
training_list = [
    ["Density", osf.density_split, ElasticNet, enet_params],
    ["Random", osf.random_split,  ElasticNet, enet_params],
    ["Batch", osf.batch_split,  ElasticNet, enet_params],
    ["Geosketch", osf.geosplit,  ElasticNet, enet_params]
]

## Loop to train model on each of the ways
# to split the data
# results are appended to the results df
for split_method, Split, Model, params_list in training_list:
    X_train, X_test, y_train, y_test = Split(PN_adata)
    for params in params_list:
        model = Model(**params, max_iter=50)
        trained_enet = model.fit(X_train, y_train)
        y_pred = trained_enet.predict(X_test)
        result_df = pd.DataFrame([{'data_split_method':f'{split_method}', 
                                 'params':model.get_params ,
                                 'r2_score': r2_score(y_test, y_pred),
                                 'MAE': mean_absolute_error(y_test, y_pred)}])
    results_df = pd.concat([results_df, result_df]) 

results_df.to_csv('gpu_randomizedsearchcv_elastic_net_training_different_datasplits_results.csv')