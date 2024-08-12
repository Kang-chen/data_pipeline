#!/usr/bin/env python
# coding: utf-8


# import packages
import re
import scanpy as sc
import numpy as np
import pandas as pd
import scrublet as scr
import scipy.sparse as sp
import lamindb as ln
import mygene
import os
import argparse

# 创建解析器对象
parser = argparse.ArgumentParser(description='Process some integers.')

# 添加参数
parser.add_argument('--source_id', type=str, required=True, help='The source ID to process')

# 解析参数
args = parser.parse_args()

# 获取参数值
source_id = args.source_id

# 打印参数值
print(f'Source ID: {source_id}')


# import adata
adata = sc.read_h5ad(f'./{source_id}/outs/processed.h5ad')


# Defines a function to check whether it is an integer
def is_integer_matrix(matrix):
    if sp.issparse(matrix):
        return np.allclose(matrix.data.astype(int), matrix.data.astype(float))
    else:
        return np.allclose(matrix.astype(int), matrix.astype(float))

# Check whether a count matrix of integer type is included
count_layer_key = None

# First check if adata.X is an integer
if is_integer_matrix(adata.X):
    count_layer_key = 'X'
else:
    for layer_key in adata.layers.keys():
        if is_integer_matrix(adata.layers[layer_key]):
            count_layer_key = layer_key
            break

if count_layer_key:
    if count_layer_key != 'X':
        adata.X = adata.layers[count_layer_key].copy()
    print(f"Count matrix found in layer '{count_layer_key}' and moved to 'X'. Proceeding with analysis.")
else:
    print("Error: No integer count matrix found in the provided .h5ad file.")
    raise ValueError("No integer count matrix found in the provided .h5ad file.")


# Check if the expression matrix contains NaN or Inf values
contains_nan = np.isnan(adata.X.data).any() if sp.issparse(adata.X) else np.isnan(adata.X).any()
contains_inf = np.isinf(adata.X.data).any() if sp.issparse(adata.X) else np.isinf(adata.X).any()

# Check if there are empty genes or cells
empty_genes = (adata.X.sum(axis=0) == 0).any()
empty_cells = (adata.X.sum(axis=1) == 0).any()


# List to store processing measures
measures = []

# Handle NaN values
if contains_nan:
    adata = adata[~np.isnan(adata.X).any(axis=1), :]
    adata = adata[:, ~np.isnan(adata.X).any(axis=0)]
    measures.append("Removed cells and genes containing NaN values")

# Handle Inf values
if contains_inf:
    adata = adata[~np.isinf(adata.X).any(axis=1), :]
    adata = adata[:, ~np.isinf(adata.X).any(axis=0)]
    measures.append("Removed cells and genes containing Inf values")

# Handle empty genes
if empty_genes:
    adata = adata[:, adata.X.sum(axis=0) != 0]
    measures.append("Removed empty genes")

# Handle empty cells
if empty_cells:
    adata = adata[adata.X.sum(axis=1) != 0, :]
    measures.append("Removed empty cells")

# Print processing measures
if measures:
    print("Processing measures:")
    for measure in measures:
        print("-", measure)
else:
    print("Data is normal, no processing needed")


# Print basic information about the data
num_genes = adata.shape[1]
num_cells = adata.shape[0]
is_expression_matrix = isinstance(adata.X, (np.ndarray, sp.csr_matrix, sp.csc_matrix))
mean_expression = adata.X.mean()
median_expression = np.median(adata.X.data) if sp.issparse(adata.X) else np.median(adata.X)

print(f"Number of genes: {num_genes}")
print(f"Number of cells: {num_cells}")
print(f"Is expression matrix: {is_expression_matrix}")
print(f"Mean expression: {mean_expression}")
print(f"Median expression: {median_expression}")


# Identify species and gene forms
gene_names = adata.var_names
is_symbol = not gene_names[0].startswith("ENSG") and not gene_names[0].startswith("ENSMUS")

# Species are identified by the genetic symbol or the Ensemble ID prefix
if is_symbol:
    human_genes = sum(gene.startswith(("A1BG", "A2M", "NAT1")) for gene in gene_names) > 0  # example human genes
    mouse_genes = sum(gene.startswith(("Atp5g1", "Actb", "Gapdh")) for gene in gene_names) > 0  # example mouse genes
else:
    human_genes = sum(gene.startswith("ENSG") for gene in gene_names) > 0
    mouse_genes = sum(gene.startswith("ENSMUS") for gene in gene_names) > 0

species = None
if human_genes > mouse_genes:
    species = "human"
elif mouse_genes > human_genes:
    species = "mouse"
else:
    species = "unknown"

print(f"The species of the sequencing data is: {species}")


# Use mygene for gene name conversion
mg = mygene.MyGeneInfo()

if is_symbol:
    if species == "human":
        query_result = mg.querymany(adata.var.index.tolist(), scopes='symbol', fields='ensembl.gene', as_dataframe=True, species='human', build='grch37')
    elif species == "mouse":
        query_result = mg.querymany(adata.var.index.tolist(), scopes='symbol', fields='ensembl.gene', as_dataframe=True, species='mouse', build='grcm38')
else:
    # If it is already an Ensemble ID, you do not need to convert it
    query_result = [{"query": gene, "ensembl": {"gene": gene}} for gene in gene_names]


# Reset index, remove weight, fill in missing values
query_result = query_result.reset_index().drop_duplicates(subset='query')
query_result['ensembl.gene'] = query_result['ensembl.gene'].fillna('0')


# Create a dictionary for mapping gene symbols to Ensembl gene ids
gene_map = query_result.set_index('query')['ensembl.gene'].to_dict()


# Update adata.var.index
adata.var.index = adata.var.index.map(gene_map).fillna('0')


# Filter out genes that fail to annotate
adata = adata[:, adata.var.index != '0']
adata.var.index = adata.var.index.map(gene_map).fillna('0')


# Remove duplicate genes
# var_names are indexes, so reset_index and drop_duplicates are required
adata.var['ensembl_id'] = adata.var.index
adata.var.reset_index(inplace=True)
adata.var.drop_duplicates(subset='ensembl_id', keep='first', inplace=True)
adata.var.set_index('index', inplace=True)

# Handle duplicates in var_names (if any)
# Remove duplicates from var_names or optionally use the corresponding Ensembl ID as the new var_names
adata = adata[:, ~adata.var_names.duplicated()]


# Make sure var_names are synchronized with var.index
adata.var_names = adata.var.index

# Assign the processed data to filtered_adata
filtered_adata = adata


# Test the quality of the data after the gene conversion
num_genes = filtered_adata.shape[1]
is_primary = num_genes > 10000
filtered_adata.obs['is_primary1'] = is_primary



# Scanpy preprocessing steps
# Filter cells and genes
sc.pp.filter_cells(filtered_adata, min_genes=200)
sc.pp.filter_genes(filtered_adata, min_cells=3)


# Deal with mitochondrial genes, according to species
if species == "human":
    mt_genes = [
    "ENSG00000198888", "ENSG00000198763", "ENSG00000198804", "ENSG00000198712",
    "ENSG00000198886", "ENSG00000212907", "ENSG00000198899", "ENSG00000198938",
    "ENSG00000198727", "ENSG00000198840", "ENSG00000210049", "ENSG00000198786",
    "ENSG00000198695", "ENSG00000198711", "ENSG00000198848", "ENSG00000198886",
    "ENSG00000198888", "ENSG00000198888", "ENSG00000210100", "ENSG00000211459",
    "ENSG00000210082", "ENSG00000210107", "ENSG00000198920", "ENSG00000198899",
    "ENSG00000198911", "ENSG00000198804", "ENSG00000210082", "ENSG00000198786",
    "ENSG00000198727", "ENSG00000198695", "ENSG00000198886", "ENSG00000198888",
    "ENSG00000198888", "ENSG00000198888", "ENSG00000198886", "ENSG00000198786",
    "ENSG00000198712", "ENSG00000198899", "ENSG00000210107", "ENSG00000198804",
    "ENSG00000198840"
    ]
elif species == "mouse":
    mt_genes = [
    "ENSMUSG00000064341", "ENSMUSG00000064345", "ENSMUSG00000064351", "ENSMUSG00000064354",
    "ENSMUSG00000064356", "ENSMUSG00000064358", "ENSMUSG00000064363", "ENSMUSG00000064365",
    "ENSMUSG00000064367", "ENSMUSG00000064370", "ENSMUSG00000064373", "ENSMUSG00000064376",
    "ENSMUSG00000064377", "ENSMUSG00000064378", "ENSMUSG00000064379", "ENSMUSG00000064380",
    "ENSMUSG00000064381", "ENSMUSG00000064382", "ENSMUSG00000064383", "ENSMUSG00000064384",
    "ENSMUSG00000064385", "ENSMUSG00000064386", "ENSMUSG00000064387", "ENSMUSG00000064388",
    "ENSMUSG00000064389", "ENSMUSG00000064390", "ENSMUSG00000064391", "ENSMUSG00000064392",
    "ENSMUSG00000064393", "ENSMUSG00000064394", "ENSMUSG00000064395", "ENSMUSG00000064396",
    "ENSMUSG00000064397"
    ]

filtered_adata.var['mt'] = filtered_adata.var_names.isin(mt_genes)
sc.pp.calculate_qc_metrics(filtered_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# Filter out cells with a high proportion of mitochondrial genes
filtered_adata = filtered_adata[filtered_adata.obs.pct_counts_mt < 5, :]


# Remove doublets using scrublet
import scrublet as scr
scrub = scr.Scrublet(filtered_adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# Add the doublet information to the adata object
filtered_adata.obs['doublet_scores'] = doublet_scores
filtered_adata.obs['predicted_doublets'] = predicted_doublets

# Remove predicted doublets
filtered_adata = filtered_adata[~filtered_adata.obs['predicted_doublets'], :]

# Print some statistics about the results of two-cell recognition
print(f"Detected doublets: {np.sum(predicted_doublets)}")
print(f"Doublet rate: {np.mean(predicted_doublets) * 100:.2f}%")
print(f"After removing doublets, data contains {filtered_adata.n_obs} cells and {filtered_adata.n_vars} features.")


# The number of features per cell was calculated
filtered_adata.obs['n_features'] = (filtered_adata.X > 0).sum(axis=1)


# Count the number of cells in each sample
sample_cell_counts = filtered_adata.obs['orig.ident'].value_counts()

# Add the number of sample cells to the obs
filtered_adata.obs['sample_cell_count'] = filtered_adata.obs['orig.ident'].map(sample_cell_counts)
filtered_adata.obs['sample_cell_count'] = filtered_adata.obs['sample_cell_count'].astype(int)

# Add the is_primary column
filtered_adata.obs['is_primary2'] = (
    (filtered_adata.obs['n_features'] >= 400) & 
    (filtered_adata.obs['sample_cell_count'] >= 1000)
)

# Output result check
# print(filtered_adata.obs[['n_features', 'sample_cell_count', 'is_primary2']])


# Comprehensively judge the data quality
filtered_adata.obs['is_primary'] = (
    (filtered_adata.obs['is_primary1']) & 
    (filtered_adata.obs['is_primary2'])
)


# Update adata
adata = filtered_adata


# Print basic information after preprocessing
print(f"Post-processing number of genes: {adata.shape[1]}")
print(f"Post-processing number of cells: {adata.shape[0]}")
print(f"Post-processing mean expression: {adata.X.mean()}")


# Print loci
# Subsites of 2.5% and 97.5% of cell counts and features were counted
# Convert a sparse matrix to a dense matrix
X_dense = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
# Subsites of 2.5% and 97.5% of cell counts and features were counted
cell_counts = np.sum(X_dense, axis=1)
cell_features = np.sum(X_dense > 0, axis=1)

counts_2_5_percentile = np.percentile(cell_counts, 2.5)
counts_97_5_percentile = np.percentile(cell_counts, 97.5)
features_2_5_percentile = np.percentile(cell_features, 2.5)
features_97_5_percentile = np.percentile(cell_features, 97.5)

print(f"Cell counts 2.5th percentile: {counts_2_5_percentile}")
print(f"Cell counts 97.5th percentile: {counts_97_5_percentile}")
print(f"Cell features 2.5th percentile: {features_2_5_percentile}")
print(f"Cell features 97.5th percentile: {features_97_5_percentile}")


path = './dataforload/kang_processing.h5ad'

dir_path = os.path.dirname(path)

if not os.path.exists(dir_path):
    os.makedirs(dir_path)
    print(f"The directory {dir_path} was created.")

# save data
adata.write('./dataforload/kang_processing.h5ad')





