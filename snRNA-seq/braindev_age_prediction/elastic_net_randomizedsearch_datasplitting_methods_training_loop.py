### Script to tune ElasticNet Hyperparameters on 
### different methods of splitting brain-dev dataset
## Import libraries
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.linear_model import ElasticNet
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
import anndata_object_splitting_functions as osf 

sc.logging.print_header()

## read in anndata
adata = sc.read('2023-05-31_braindev_hg38_post-ds1000cts-clustered.h5ad')

# subset adata to principal neurons only
PNs_annots = ['L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev']
PN_adata = adata[adata.obs['major_clust'].isin(PNs_annots)].copy()

## initialize model
model = ElasticNet(random_state=1)

## set up parameter grid
grid = dict()
grid['alpha'] = np.arange(0, 1, 0.1)
grid['l1_ratio'] = np.arange(0, 1, 0.1)

# initialize gridsearch object
gm = RandomizedSearchCV(model, param_distributions=grid, cv=5, n_iter= 50, n_jobs=-1)

# create dataframe for storing results
results_df = pd.DataFrame(columns = ['data_split_method',
                                     'best_params', 'r2_score', 'MAE'])

## Create 2D list of names of data splitting methods,
# splitting functions, and the model to be trained
training_list = [
    ["Density", osf.density_split, gm],
    ["Random", osf.random_split, gm],
    ["Batch", osf.batch_split, gm],
    ["Geosketch", osf.geosplit, gm]
]

## Loop to train model on each of the ways
# to split the data
# results are appended to the results df
for split_method, Split, model in training_list:
    X_train, X_test, y_train, y_test = Split(PN_adata)
    gm.fit(X_train, y_train)
    y_pred = gm.predict(X_test)
    result_df = pd.DataFrame([{'data_split_method':f'{split_method}', 
                                 'best_params':gm.best_params_ ,
                                 'r2_score': r2_score(y_test, y_pred),
                                 'MAE': mean_absolute_error(y_test, y_pred)}])
    results_df = pd.concat([results_df, result_df]) 

results_df.to_csv('randomizedsearchcv_elastic_net_training_different_datasplits_results.csv')