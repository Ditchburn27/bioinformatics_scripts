### Script to tune hyperparamters of logistic regression elasticnet
## on brain-dev dataset downsampled to 2500 counts, and counts binarized
## Velmeshev 2019 dataset will be used for evaluation/scoring. 
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import top_k_accuracy_score
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import matthews_corrcoef
from matplotlib.pyplot import rc_context
from sklearn.metrics import ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import anndata_object_splitting_functions as osf
import pickle

# set seed for random state
seed = 123

# import anndata objects
brain_dev = sc.read_h5ad("2023-12-11_brain_dev_ds2500cts.h5ad")
velm = sc.read_h5ad("2023-12-11_velm_ds2500cts.h5ad")

# split data into X & y for train and test
X_train, y_train = osf.xy_split(brain_dev, target_variable='stage_id',
                                binarize=True)
X_test, y_test = osf.xy_split(velm, target_variable='stage_id',
                              binarize=True)

# initialize model
eln = LogisticRegression(penalty='elasticnet', solver='saga',
                          multi_class='multinomial', n_jobs=-1, random_state=seed)

# set up parameter grid
grid = dict()
grid['C'] = np.arange(0.01, 1, 0.01)
grid['l1_ratio'] = np.arange(0.01, 1, 0.01)
grid['max_iter'] = range(100, 500)


# Initialize randomsearch object
gm = RandomizedSearchCV(eln, param_distributions=grid,
                        cv=5, random_state=seed, scoring='roc_auc_ovo', n_jobs=-1)

# Model training
gm.fit(X_train, y_train)

print(gm.best_params_)
print(gm.best_score_)

# Using best performing model
# Evaluating performance
best_model = gm.best_estimator_
best_model.fit(X_train, y_train)

y_pred = best_model.predict(X_test)
y_proba = best_model.predict_proba(X_test)


class_labels = ['Adolescence', 'Adult', 'Childhood', 'Fetal', 'Infancy', 'Neonatal']

result_df = pd.DataFrame([{'mean_accuracy': best_model.score(X_test, y_test),
                                'roc_auc': roc_auc_score(y_test, y_proba, multi_class='ovo', labels=class_labels),
                                                    'balanced_accuracy': balanced_accuracy_score(y_test, y_pred),
                                                    'matthews_corrcoef':matthews_corrcoef(y_test, y_pred),
                                                      'top_k_accuracy': top_k_accuracy_score(y_test, y_proba, labels=class_labels),
                                                      'cohen_kappa':cohen_kappa_score(y_test, y_pred, labels=class_labels)}])

result_df

for param, value in best_model.get_params(deep=True).items():
    print(f"{param} -> {value}")

# Save trained model
filename = '2023-12-18_PN_logisticRegression_ELN_stage_id_classifier'
pickle.dump(best_model, open(filename, 'wb'))


