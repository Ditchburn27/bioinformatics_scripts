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
grid['alpha'] = [arange(0, 1, 0.01), 1, 5, 100]
grid['l1_ratio'] = arange[0, 1, 0.01]

# initialize gridsearch object
gm = GridSearchCV(model, param_grid=grid, cv=5)

## Density split data
X_train, X_test, y_train, y_test = osf.density_split(PN_adata)
# Train on density split data
gm.fit(X_train, y_train)
# predict
y_pred = gm.predict(X_test)
# scoring
print(f'Best Paramaters:{gm.best_params}')
print(f'r2 score:{r2_score(y_test, y_pred)}')
print(f'MAE:{mean_absolute_error(y_test, y_pred)}')

###############################
## Random split data
X_train, X_test, y_train, y_test = osf.random_split(PN_adata)
# Train on density split data
gm.fit(X_train, y_train)
# predict
y_pred = gm.predict(X_test)
# scoring
print(f'Best Paramaters:{gm.best_params}')
print(f'r2 score:{r2_score(y_test, y_pred)}')
print(f'MAE:{mean_absolute_error(y_test, y_pred)}')

###############################
## Batch split data
X_train, X_test, y_train, y_test = osf.batch_split(PN_adata)
# Train on density split data
gm.fit(X_train, y_train)
# predict
y_pred = gm.predict(X_test)
# scoring
print(f'Best Paramaters:{gm.best_params}')
print(f'r2 score:{r2_score(y_test, y_pred)}')
print(f'MAE:{mean_absolute_error(y_test, y_pred)}')

###############################
## Geosketch split data
X_train, X_test, y_train, y_test = osf.geosplit(PN_adata)
# Train on density split data
gm.fit(X_train, y_train)
# predict
y_pred = gm.predict(X_test)
# scoring
print(f'Best Paramaters:{gm.best_params}')
print(f'r2 score:{r2_score(y_test, y_pred)}')
print(f'MAE:{mean_absolute_error(y_test, y_pred)}')