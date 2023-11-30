### Script to tune ElasticNet hyperparameters
## on brain_dev dataset with genes intersecting 
## with Velmeshev
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNet
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import RandomizedSearchCV
import anndata_object_splitting_functions as osf 

## read in anndata
adata = sc.read('2023-11-30_braindev_hg38_aligned_with_velnm_intersecting_genes.h5ad')

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
gm = RandomizedSearchCV(model, param_distributions=grid, cv=5, n_iter= 50, 
                        n_jobs=-1, scoring='r2')

# split data into train/test
X_train, X_test, y_train, y_test = osf.geosplit(PN_adata)

# Run model training
gm.fit(X_train, y_train)
y_pred = gm.predict(X_test)
result_df = pd.DataFrame([{'best_params':gm.best_params_ ,
                                'r2_score': r2_score(y_test, y_pred),
                                'MAE': mean_absolute_error(y_test, y_pred)}]) 

result_df.to_csv('2023-11-30_elasticnet_tuning_on_intersecting_genes.csv', sep=',')