### Script to preprocess snRNA-seq data

############################################################################################
###  Import CLI libraries
import os
import argparse

############################################################################################
### Deal with command line arguments
arg_description = '''Preprocess snRNA-seq Data:
                This script cleans snRNA-seq matrices by,
                   removing doublets, and filtering for
                high ribosomal and mitochondral % nuclei'''
parser = argparse.ArgumentParser(description=arg_description)
parser.add_argument('file', type=str, help='Path to the input 10X h5 matrix file')
parser.add_argument('-o', '--output_filename', type=str, required=True,
                    help='Provide a name for the output files')
args = parser.parse_args()

data = args.file
filename = args.output_filename

#############################################################################################
### Import libraries snRNA-seq data processing
import scanpy as sc
import anndata as ad
import pandas as pd 
import numpy as np
import scrublet as scr
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import rc_context
sc.logging.print_header()

############################################################################################
### Make and set paths to directories for data and figures
data_path = "./out_data/"
fig_path = "./out_figs/"

os.mkdir(data_path)
os.mkdir(fig_path)

############################################################################################
### Read in count matrix
adata = sc.read_10x_h5(data)
adata.var_names_make_unique()

#############################################################################################
### QC nuclei
# Get scanpy quality metrics
sc.pp.calculate_qc_metrics( adata, 
                           percent_top=None, inplace=True)
# Remove genes with fewer than 5 total counts 
sc.pp.filter_genes( adata, min_counts=5, inplace=True)

#############################################################################################
### Remove doublets with scrublet
# Remove doublets on a per sample basis
dub_mk = np.array( [False]*adata.shape[0], 
                  dtype=bool)
dub_sc = np.array( [False]*adata.shape[0],
                   dtype=float)
scrub = scr.Scrublet(adata.X.tocsc(), 
                     expected_doublet_rate=0.1)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3,
                                                         min_gene_variability_pctl=85,
                                                         n_prin_comps=10,
                                                         verbose=False)
dub_mk = predicted_doublets
dub_sc = doublet_scores
# add scrublet doublet scores to adata
adata.obs['doublet_score'] = dub_sc
# remove predicted doublets
adata = adata[~dub_mk]
# recalculate qc metrics
sc.pp.calculate_qc_metrics( adata, percent_top=None, inplace=True)

###############################################################################################
### Library sizes and feature counts per sample 
adata.obs['log10_gene_counts'] = np.log10( adata.obs['n_genes_by_counts'])
adata.obs['log10_UMI_counts']  = np.log10( adata.obs['total_counts'])

with rc_context({'figure.figsize': (5,5)}):
    ax = sc.pl.violin( adata, ['log10_gene_counts'], stripplot=True, 
                  inner='box', rotation=70, show=False)
    plt.savefig(f"{fig_path}{filename}_violin_feature-cts.png", format='png', bbox_inches='tight')

with rc_context({'figure.figsize': (5,5)}):
    ax = sc.pl.violin( adata, ['log10_UMI_counts'], stripplot=True, 
                  inner='box', rotation=70, show=False)
    plt.savefig(f"{fig_path}{filename}_violin_library-sizes.png", format='png', bbox_inches='tight')

###############################################################################################
### Remove nuclei with high ribosomal and mitochondrial count percentages
# add mitochondrial gene list to object
mito_genes = adata.var_names.str.startswith('MT-')
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# download list of ribosomal genes
ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
ribo_genes = pd.read_table( ribo_url, header=[1])
# create mask for all ribo genes
ribo_mk = np.in1d( adata.var_names.values.astype(str), ribo_genes)
adata.obs['percent_ribo'] = np.sum( adata[:,ribo_mk].X, axis=1).A1 / np.sum( adata.X, axis=1).A1
# violin plots of percentages
with rc_context({'figure.figsize':(5,5)}):
    ax = sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'percent_ribo'],
             jitter=0.4, multi_panel=True)
    plt.savefig(f"{fig_path}{filename}_mito_ribo_pct_violin.png", format='png', bbox_inches='tight')
# Outlier function for MADs calculation
from scipy.stats import median_abs_deviation
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M)
    return outlier
# find gene count outliers 
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5))
adata.obs.outlier.value_counts()
# find mitochondrial pct outliers 
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 20)
adata.obs.mt_outlier.value_counts()
# find ribosomal pct outliers
adata.obs["ribo_outlier"] = is_outlier(adata, "percent_ribo", 3) | (
    adata.obs["percent_ribo"] > 20)
adata.obs.ribo_outlier.value_counts()
# Number of cells filtered based on ribo and mito genes
print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier) & (~adata.obs.ribo_outlier)].copy()
print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
# Remove highly expressed genes as they can have outstized library effects
adata = adata[:,~adata.var_names.str.endswith('MALAT1')]
adata = adata[:, ~np.in1d(adata.var_names, np.append(mito_genes, ribo_genes))]

######################################################################################################
### Plot highest expressed genes and save anndata object
with rc_context({'figure.figsize':(5,5)}):
    ax = sc.pl.highest_expr_genes(adata, n_top=20)
    plt.savefig(f"{fig_path}{filename}_top20_highest_expressed_genes.png", format='png', bbox_inches='tight')

adata.write(f"{data_path}{filename}_cleaned_filtered.h5ad")