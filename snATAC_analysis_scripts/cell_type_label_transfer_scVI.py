### Script to transfer cell type labels from snRNA-seq to snATAC-seq
#######################################
## Import libraries
import snapatac2 as snap
import scanpy as sc
import numpy as np 
import pandas as pd
import scvi
import anndata as ad
import argparse
import matplotlib.pyplot as plt

######################################
## CLI inputs
arg_description = '''Script to transfer cell type labels
 from snRNA-seq to snATAC-seq, using scVI'''

parser = argparse.ArgumentParser(description=arg_description)
parser.add_argument('snRNA_h5ad', type=str, help='''Path to
                    annotated snRNA-seq anndata object''')
parser.add_argument('snATAC_h5ad', type=str, help='''Path to 
                    snATAC anndata object''')
args = parser.parse_args()

reference_file = args.snRNA_h5ad
atac_file = args.snATAC_h5ad

######################################
## Get objects ready to be trained
# Make gene activity matrix
print('Reading in datasets and making gene activity matrix...')
reference = sc.read_h5ad(reference_file)
reference.obs['cell_type'] = reference.obs['major_clust']
atac = sc.read_h5ad(atac_file)
query = snap.pp.make_gene_matrix(atac, gene_anno=snap.genome.GRCh38)
query.obs['cell_type'] = pd.NA

# Merge RNA & gene activity matrix
data = ad.concat([reference, query], join = 'inner', label = 'dataset', 
                 keys=['reference', 'query'], index_unique='_')

# Filter cells & genes
sc.pp.filter_genes(data, min_cells=5)
sc.pp.highly_variable_genes(data, n_top_genes=5000, flavor='seurat_v3', 
                            batch_key='dataset', subset=True)

######################################
# Pre-train model
print("Pre-training the model...")
scvi.model.SCVI.setup_anndata(data, batch_key='dataset')
vae = scvi.model.SCVI(data, n_layers=2, n_latent=365, gene_likelihood='nb', 
                      dispersion='gene-batch')
vae.train(max_epochs=1000, early_stopping=True)

# Plot elbo validation
ax = vae.history['elbo_train'][1:].plot()
vae.history['elbo_validation'].plot(ax=ax)
plt.savefig('elbo_validation_plot.png')

# Train model
print("Training the model...")
# Add observation labels
data.obs['celltype_scanvi'] = 'Unknown'
ref_idx = data.obs['batch'] == 'reference'
data.obs['celltype_scanvi'][ref_idx] = data.obs['cell_type'][ref_idx]

lvae = scvi.model.SCANVI.from_scvi_model(vae, adata=data, labels_key='celltype_scanvi', 
                                         unlabeled_category='Unknown')
lvae.train(max_epochs=1000, n_samples_per_label=100)

# Plot elbo train
lvae.history['elbo_train'][1:].plot()
plt.savefig('elbo_train_plot.png')

######################################
# Label Transfer
print("Transferring labels...")
data.obs['C_scANVI'] = lvae.predict(data)
data.obsm['X_scANVI'] = lvae.get_latent_representation(data)

# Calculate neighbors and umap for joint embedding
print("Calculating neighbors and umap for joint embedding...")
sc.pp.neighbors(data, use_rep="X_scANVI")
sc.tl.umap(data)

######################################
# Save objects
# Add predicted labels to original atac object
atac.obs['cell_type'] = data.obs.loc[atac.obs_names + '_query']['C_scANVI'].to_numpy()
# Save jointly embedded object
data.write('braindev_hg38_snRNA_snATAC_jointly_embedded.h5ad')
# Save annotated snATAC
atac.write('braindev_hg38_snATAC_annotated.h5ad')
