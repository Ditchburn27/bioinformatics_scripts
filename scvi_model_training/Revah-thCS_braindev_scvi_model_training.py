# Train SCVI model on the brain development dataset using the scArches optimisation
# Integrate snRNA-seq data from transplanted hCS (Revah, 2022)

import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown

sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)

# Read in anndata objects
adata_ref = sc.read("/scratch/ll005/lditchburn/ds-1000cts_braindev-5k-hvgs_ref_adata.h5ad")
adata_query = sc.read("/scratch/ll005/lditchburn/braindev-matching-hvgs_thCS_query-adata.h5ad")

# Creat SCVI model and train it on reference dataset 
sca.models.SCVI.setup_anndata(adata_ref)

vae = sca.models.SCVI(
    adata_ref,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)

vae.train()

# Create anndata file of latent representation and compute UMAP
reference_latent = sc.AnnData(vae.get_latent_representation())
reference_latent.obs["major_clust"] = adata_ref.obs["major_clust"]
reference_latent.obs["stage_id"] = adata_ref.obs["stage_id"] 

sc.pp.neighbors(reference_latent, n_neighbors=8)
sc.tl.leiden(reference_latent)
sc.tl.umap(reference_latent)
sc.pl.umap(reference_latent,
           color=['major_clust', 'stage_id'],
           frameon=False,
           wspace=0.6, save='reference_latent_UMAP.png'
           )

ref_path = '/scratch/ll005/lditchburn/braindev_ref_model'
vae.save(ref_path, overwrite=True)

# Perform surgery on reference model and train on query dataset 
model = sca.models.SCVI.load_query_data(
    adata_query,
    ref_path,
    freeze_dropout = True,
)

model.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))

query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['batch'] = adata_query.obs['batch']

sc.pp.neighbors(query_latent)
sc.tl.leiden(query_latent)
sc.tl.umap(query_latent)
plt.figure()
sc.pl.umap(
    query_latent,
    color=["batch"],
    frameon=False,
    wspace=0.6, save='query_latent_UMAP.png'
)

surgery_path = '/scratch/ll005/lditchburn/thCS-braindev_surgery_model'
model.save(surgery_path, overwrite=True)

adata_full = adata_ref.concatenate(adata_query)
full_latent = sc.AnnData(model.get_latent_representation(adata=adata_full))
full_latent.obs['major_clust'] = adata_full.obs['major_clust']
full_latent.obs['stage_id'] = adata_full.obs['stage_id']

sc.pp.neighbors(full_latent)
sc.tl.leiden(full_latent)
sc.tl.umap(full_latent)
plt.figure()
sc.pl.umap(
    full_latent,
    color=["major_clust", "stage_id"],
    frameon=False,
    wspace=0.6, save='ref_and_query_latent_representation_UMAP.png'
)