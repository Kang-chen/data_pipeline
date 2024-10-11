#!/usr/bin/env python
# coding: utf-8


# import packages
import re
import scanpy as sc
import numpy as np
import pandas as pd
import scrublet as scr
import scipy.sparse as sp
import mygene
import os
import argparse
import shutil
import warnings

# Ignore the warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Observation names are not unique")
warnings.filterwarnings("ignore", category=UserWarning, message="Variable names are not unique")
warnings.filterwarnings("ignore", category=UserWarning, message="Trying to modify attribute `.obs` of view, initializing view as actual")

# define rguments
def parse_args():
    parser = argparse.ArgumentParser(description="Quality control of adata.")
    parser.add_argument('--source_id', type=str, required=True, help="GEO series ID (e.g., GSE156793)")
    return parser.parse_args()


# define a function to find and move integer matrix
def find_and_move_integer_matrix(adata):
   
    def is_integer_matrix(matrix):
        if sp.issparse(matrix):
            return np.allclose(matrix.data.astype(int), matrix.data.astype(float))
        else:
            return np.allclose(matrix.astype(int), matrix.astype(float))

    # Initialize the count matrix
    count_layer_key = None

    # Check whether adata.X is an integer matrix
    if is_integer_matrix(adata.X):
        count_layer_key = 'X'
    else:
        # Check whether the matrix in adata.layers is an integer matrix
        for layer_key in adata.layers.keys():
            if is_integer_matrix(adata.layers[layer_key]):
                count_layer_key = layer_key
                break

    # If the integer count matrix is found, move it
    if count_layer_key:
        if count_layer_key != 'X':
            adata.X = adata.layers[count_layer_key].copy()
        print(f"Count matrix found in layer '{count_layer_key}' and moved to 'X'. Proceeding with analysis.")
    else:
        raise ValueError("No integer count matrix found in the provided .h5ad file.")
    
    return adata 


# define a function to clean expression matrix
def clean_expression_matrix(adata):
    
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
        if sp.issparse(adata.X):
            nan_rows = np.isnan(adata.X.data).any(axis=1)
            nan_cols = np.isnan(adata.X.data).any(axis=0)
            adata = adata[~nan_rows, :]
            adata = adata[:, ~nan_cols]
        else:
            adata = adata[~np.isnan(adata.X).any(axis=1), :]
            adata = adata[:, ~np.isnan(adata.X).any(axis=0)]
        measures.append("Removed cells and genes containing NaN values")

    # Handle Inf values
    if contains_inf:
        if sp.issparse(adata.X):
            inf_rows = np.isinf(adata.X.data).any(axis=1)
            inf_cols = np.isinf(adata.X.data).any(axis=0)
            adata = adata[~inf_rows, :]
            adata = adata[:, ~inf_cols]
        else:
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
    return adata


# define a function toidentify species
def identify_species(adata):
    gene_names = adata.var_names
    is_symbol = not gene_names[0].startswith("ENSG") and not gene_names[0].startswith("ENSMUS")

    if is_symbol:
        human_genes = sum(gene.startswith(("A1BG", "A2M", "NAT1")) for gene in gene_names)
        mouse_genes = sum(gene.startswith(("Atp5g1", "Actb", "Gapdh")) for gene in gene_names)
    else:
        human_genes = sum(gene.startswith("ENSG") for gene in gene_names)
        mouse_genes = sum(gene.startswith("ENSMUS") for gene in gene_names)
    
    return 'human' if human_genes > mouse_genes else 'mouse'


# define a function to preprocess data
def preprocess_data(adata, species):

    # filter cells and genes
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # deal with mitochondrial genes, according to species
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
    
    # labeled mitochondrial gene
    adata.var['mt'] = adata.var_names.isin(mt_genes)

    # calculate mitochondrial gene ratio
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # filter cells according to mitochondrial gene ratio
    adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

    return adata


# define a function to ensure batch column
def ensure_batch_column(adata):
    if 'batch' not in adata.obs.columns:
        adata.obs['batch'] = 1
    return adata


# define a function to remove doublets and process
def remove_doublets_and_process(adata):

    # Remove doublets using scrublet
    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # Add the doublet information to the adata object
    adata.obs['doublet_scores'] = doublet_scores
    adata.obs['predicted_doublets'] = predicted_doublets

    # Print some statistics about the results of two-cell recognition
    # print(f"Doublet rate: {np.mean(predicted_doublets) * 100:.2f}%")

    # Remove predicted doublets
    adata = adata[~adata.obs['predicted_doublets'], :].copy()

    # Print some statistics about the results of two-cell recognition
    print(f"After removing doublets, data contains {adata.n_obs} cells and {adata.n_vars} features.")

    # The number of features per cell was calculated
    adata.obs['n_features'] = (adata.X > 0).sum(axis=1)

    # Count the number of cells in each sample
    adata = ensure_batch_column(adata)
    sample_cell_counts = adata.obs['batch'].value_counts()

    # The number of cells in the sample is mapped back into the obs
    adata.obs['sample_cell_count'] = adata.obs['batch'].map(sample_cell_counts)
    adata.obs['sample_cell_count'] = adata.obs['sample_cell_count'].astype(int)

    return adata


# define a function to evaluate data quality
def evaluate_data_quality(adata, gene_threshold=10000, min_features=400, min_sample_cells=1000):
    
    # test the quality of the data after the gene conversion
    num_genes = adata.shape[1]
    criterion1 = num_genes > gene_threshold
    adata.obs['criterion1'] = criterion1

    # Add the is_primary column
    adata.obs['criterion2'] = (
        (adata.obs['n_features'] >= min_features) &
        (adata.obs['sample_cell_count'] >= min_sample_cells)
    )

    # Comprehensively judge the data quality
    adata.obs['is_primary'] = (
        (adata.obs['criterion1']) &
        (adata.obs['criterion2'])
    )
    
    return adata


# define a function to convert csc to csr
def convert_csc_to_csr(adata):

    if sp.issparse(adata.X):
        if isinstance(adata.X, sp.csc_matrix):
            adata.X = adata.X.tocsr()
            print("adata.X has been converted from CSC to CSR format.")
    
    return adata


# save adata
def save_adata(adata):
    
    save_dir = './dataforload'
    save_path = os.path.join(save_dir, 'qc_adata.h5ad')

    # create folders
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # save adata
    adata.write(save_path)


if __name__ == "__main__":
    
    # get arguments
    args = parse_args() 

    # import adata
    adata = sc.read_h5ad('./process/merged_adata.h5ad')

    # process adata
    adata = find_and_move_integer_matrix(adata)
    adata = clean_expression_matrix(adata)

    # process adata
    species = identify_species(adata)
    adata = preprocess_data(adata, species)
    adata = remove_doublets_and_process(adata)
    adata = evaluate_data_quality(adata)
    adata = convert_csc_to_csr(adata)

    # save adata
    adata.obs['dataset_id'] = args.source_id
    adata.obs['organisms'] = species
    save_adata(adata)

    # delete folders
    shutil.rmtree('./process')


