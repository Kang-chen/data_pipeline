#!/usr/bin/env python
# coding: utf-8

# import packages
import scanpy as sc
import pandas as pd
import anndata as ad
import os
import ftplib
import gzip
import shutil
import tarfile
import argparse
import warnings

# Ignore the warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Observation names are not unique")
warnings.filterwarnings("ignore", category=UserWarning, message="Variable names are not unique")
warnings.filterwarnings("ignore", category=UserWarning, message="Trying to modify attribute `.obs` of view, initializing view as actual")


# define rguments
def parse_args():
    parser = argparse.ArgumentParser(description="Process GEO data.")
    parser.add_argument('--source_id', type=str, required=True, help="GEO series ID (e.g., GSE156793)")
    parser.add_argument('--output_dir', type=str, default='./geo_supp_files', help="Output directory for GEO files")
    return parser.parse_args()


# Download supplementary files in GEO
def download_geo_supp_files(geo_id, output_dir):
    geo_series = geo_id[:-3] + 'nnn'
    ftp_url = f"ftp.ncbi.nlm.nih.gov"
    ftp_dir = f"/geo/series/{geo_series}/{geo_id}/suppl/"
    
    # create output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # connect ftp
    ftp = ftplib.FTP(ftp_url)
    ftp.login()
    
    # cd output_dir
    ftp.cwd(ftp_dir)

    # file lists
    files = ftp.nlst()

    # download
    for file in files:
        local_file_path = os.path.join(output_dir, file)
        with open(local_file_path, 'wb') as local_file:
            ftp.retrbinary(f"RETR {file}", local_file.write)

    # shut down ftp
    print(f"Downloaded: All {geo_id} files.")
    ftp.quit()


# Detect tar files
def check_and_extract_raw_tar(output_dir):
    
    files = os.listdir(output_dir)
    tar_files = [f for f in files if 'raw' in f.lower() and f.endswith('.tar')]

    if not tar_files:
        return

    for tar_file in tar_files:
        tar_path = os.path.join(output_dir, tar_file)
        
        # Open tar
        if tarfile.is_tarfile(tar_path):
            with tarfile.open(tar_path, "r") as tar:
                tar.extractall(path=output_dir)
            print(f"Extraction complete for {tar_file}.")


# Detect tar.gz files for sub files
def check_and_extract_tar_gz(output_dir):
    files = os.listdir(output_dir)
    
    tar_gz_files = [f for f in files if f.endswith('.tar.gz')]
    
    if not tar_gz_files:
        return
    
    for tar_gz_file in tar_gz_files:
        tar_gz_path = os.path.join(output_dir, tar_gz_file)
        
        # Open tar.gz
        if tarfile.is_tarfile(tar_gz_path):
            with tarfile.open(tar_gz_path, "r:gz") as tar:
                tar.extractall(path=output_dir)
            print(f"Extraction complete for {tar_gz_file}.")


# remove header
def check_and_remove_header(file_path, file_type):
    lines = []
    with gzip.open(file_path, 'rt') as f:
        for i, line in enumerate(f):
            lines.append(line)
            if i == 0:
                first_line = line.strip().split('\t')

                # for features
                if file_type == 'features':
                    header_keywords = ['gene', 'id', 'x']
                    if any(keyword in first_line[0].lower() for keyword in header_keywords):
                        lines = lines[1:]  # remove header
                
                # for barcode
                elif file_type == 'barcodes':
                    header_keywords = ['barcode', 'id', 'x']
                    if any(keyword in first_line[0].lower() for keyword in header_keywords):
                        lines = lines[1:]  # remove header
    
    # overwrite
    with gzip.open(file_path, 'wt') as f_out:
        for line in lines:
            f_out.write(line)


# check and update features
def check_and_update_features(file_path):
    
    with gzip.open(file_path, 'rt') as f:
        features = pd.read_csv(f, sep='\t', header=None)

    # Determine if there are only one columns
    if features.shape[1] == 1:

        # Add a first and third column
        features[1] = range(1, len(features) + 1)
        features[2] = 'Gene Expression'
        new_order = [1, 0, 2]
        features = features.iloc[:, new_order]
        with gzip.open(file_path, 'wt') as f_out:
            features.to_csv(f_out, sep='\t', index=False, header=False)

    
    # Determine if there are only two columns
    if features.shape[1] == 2:

        # Add a third column with the value 'Gene Expression'
        features[2] = 'Gene Expression'
        # Overwrite source file
        with gzip.open(file_path, 'wt') as f_out:
            features.to_csv(f_out, sep='\t', index=False, header=False)


# identify species
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


# filter adata by species
def filter_adata_by_species(adata_list):
    human_count = 0
    mouse_count = 0
    human_adata = []
    mouse_adata = []
    
    # identify species
    for adata in adata_list:
        species = identify_species(adata)
        if species == 'human':
            human_count += 1
            human_adata.append(adata)
        else:
            mouse_count += 1
            mouse_adata.append(adata)
    
    # Retaining adata according to counts
    if human_count > mouse_count:
        print("Retaining human samples only.")
        return human_adata
    else:
        print("Retaining mouse samples only.")
        return mouse_adata


# convert gene symbols to ensembl
def convert_symbol_to_ensembl(adata, human_mapping="../python/geneid/human_grch38_symbol_ensembl.csv", mouse_mapping="../python/geneid/mouse_grcm39_symbol_ensembl.csv"):
    
    # Determine the mapping file based on the species
    species = identify_species(adata)
    if species.lower() == "human":
        mapping_file = human_mapping
    elif species.lower() == "mouse":
        mapping_file = mouse_mapping
    else:
        raise ValueError("Species must be either 'human' or 'mouse'.")
    
    # Load the appropriate mapping file
    ensembl_mapping = pd.read_csv(mapping_file)
    ensembl_mapping.columns = ['EnsemblID', 'GeneSymbol']
    ensembl_mapping.set_index('GeneSymbol', inplace=True)
    
    # Check if adata.var.index is in EnsemblID format
    if adata.var.index.str.startswith('ENSG').all() or adata.var.index.str.startswith('ENSMUSG').all():
        return adata
    else:       
        # Extract gene symbols from adata.var.index
        gene_symbols = adata.var.index
        matched_ensembl_ids = ensembl_mapping.reindex(gene_symbols)['EnsemblID']
        
        # Filter out genes that could not be matched to Ensembl IDs
        valid_ensembl_ids = matched_ensembl_ids.dropna()

        # Remove duplicates
        valid_ensembl_ids = valid_ensembl_ids[~valid_ensembl_ids.index.duplicated(keep='first')]
        adata = adata[:, ~adata.var.index.duplicated(keep='first')]
        
        # Retain only the genes that could be matched successfully
        adata = adata[:, valid_ensembl_ids.index]
        adata.var.index = valid_ensembl_ids.values

        print(f"Successfully converted Gene Symbols to Ensembl IDs, retained {len(valid_ensembl_ids)} genes.")
        
        return adata



# detect adata_list by gene number
def filter_adata_by_gene_count(adata_list, gene_threshold=3000):

    filtered_adata_list = []
    
    for adata in adata_list:
        gene_count = adata.var.shape[0] 
        
        # save object more than 3000
        if gene_count > gene_threshold:
            filtered_adata_list.append(adata)
    
    return filtered_adata_list


# check and transpose adata
def check_and_transpose_adata(adata):
    gene_prefix_human = "ENSG"
    gene_prefix_mouse = "ENSMUS"
    human_gene_symbols = ["A1BG", "A2M", "NAT1"]  # Extendable list of common human gene symbols
    mouse_gene_symbols = ["Atp5g1", "Actb", "Gapdh"]  # Extendable list of common mouse gene symbols

    # Check if var_names or obs_names contain gene names (symbols or Ensembl IDs)
    def contains_gene_names(name_list):
        # Check if Ensembl ID or common gene symbols are present
        return (any(name.startswith(gene_prefix_human) or name.startswith(gene_prefix_mouse) for name in name_list) or
                any(name in human_gene_symbols for name in name_list) or
                any(name in mouse_gene_symbols for name in name_list))

    # Check if var_names and obs_names contain gene names
    obs_has_genes = contains_gene_names(adata.obs_names)

    # If obs_names contain gene names, transpose the data
    if obs_has_genes:
        adata = adata.transpose()

    return adata


# filter single adata by gene count
def filter_single_adata_by_gene_count(adata, gene_threshold=3000):
    gene_count = adata.var.shape[0]
    if gene_count > gene_threshold:
        return adata
    else:
        return None 


# read and merge 10x files
def read_and_merge_mtx_files(output_dir):
    
    files_and_dirs = os.listdir(output_dir)
    files = [f for f in files_and_dirs if os.path.isfile(os.path.join(output_dir, f))]
    dirs = [d for d in files_and_dirs if os.path.isdir(os.path.join(output_dir, d))]
    
    features_files = [f for f in files if 'features.tsv.gz' in f or 'genes.tsv.gz' in f]
    barcodes_files = [f for f in files if 'barcodes.tsv.gz' in f or 'identities.tsv.gz' in f]
    matrix_files = [f for f in files if 'mtx.gz' in f]
    
    if not (features_files and barcodes_files and matrix_files and not dirs):
        if dirs:
            print("Detected folders, processing them with read_10x_mtx...")
        else:
            print("Required feature, barcode, and matrix files not found.")
            return None

    # get prefixes
    prefixes = set([f.rsplit('_', 1)[0] for f in matrix_files])
    
    adata_list = []
    
    for prefix in prefixes:
        feature_file = [f for f in features_files if f.startswith(prefix)][0]
        barcode_file = [f for f in barcodes_files if f.startswith(prefix)][0]
        matrix_file = [f for f in matrix_files if f.startswith(prefix)][0]

        # create folders
        folder_name = os.path.join(output_dir, prefix.strip('_'))
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        # move flies and rename
        shutil.move(os.path.join(output_dir, feature_file), os.path.join(folder_name, 'features.tsv.gz'))
        shutil.move(os.path.join(output_dir, barcode_file), os.path.join(folder_name, 'barcodes.tsv.gz'))
        shutil.move(os.path.join(output_dir, matrix_file), os.path.join(folder_name, 'matrix.mtx.gz'))
        
        # check features
        check_and_remove_header(os.path.join(folder_name, 'features.tsv.gz'), file_type='features')
        check_and_remove_header(os.path.join(folder_name, 'barcodes.tsv.gz'), file_type='barcodes')
        check_and_update_features(os.path.join(folder_name, 'features.tsv.gz'))

        # read 10x files
        print(f"Reading {folder_name}...")
        adata = sc.read_10x_mtx(folder_name, var_names='gene_symbols')
        adata = filter_single_adata_by_gene_count(adata)
        if adata is None:
            continue
        adata = convert_symbol_to_ensembl(adata)
        adata_list.append(adata)

    # directly read folders
    for folder in dirs:
        folder_path = os.path.join(output_dir, folder)
        print(f"Directly reading from folder {folder_path}...")
        adata = sc.read_10x_mtx(folder_path, var_names='gene_symbols')
        adata = filter_single_adata_by_gene_count(adata)
        if adata is None:
            continue
        adata = convert_symbol_to_ensembl(adata)
        adata_list.append(adata)
        
    # merge adata
    adata_list = filter_adata_by_species(adata_list)
    adata_list = filter_adata_by_gene_count(adata_list)

    if len(adata_list) > 1:
        print("Merging AnnData objects...")
        adata_merged = ad.concat(adata_list, label="batch")
    else:
        adata_merged = adata_list[0]

    return adata_merged


# read and merge h5 files
def read_and_merge_h5_files(output_dir):
    files = os.listdir(output_dir)

    # Detect .h5 and .h5 files
    exclude_keywords = ['meta', 'process', 'normalized', 'annotations', 'celltype', 'adt']
    h5_files = sorted([f for f in files if f.endswith('.h5')
                       and not any(exclude in f.lower() for exclude in exclude_keywords)])
    
    if not h5_files:
        return None
    
    adata_list = []
    
    # read h5 file one by one
    for h5_file in h5_files:
        file_path = os.path.join(output_dir, h5_file)
        print(f"Reading {h5_file}...")
        adata = sc.read_10x_h5(file_path)
        adata = check_and_transpose_adata(adata)
        adata = convert_symbol_to_ensembl(adata)
        adata_list.append(adata)
    
    # merge adata
    adata_list = filter_adata_by_species(adata_list)
    adata_list = filter_adata_by_gene_count(adata_list)

    #  merge h5 file
    if len(adata_list) > 1:
        print("Merging AnnData objects...")
        adata_merged = ad.concat(adata_list, label="batch")
    else:
        adata_merged = adata_list[0]
    
    return adata_merged


# read and merge h5ad files
def read_and_merge_h5ad_files(output_dir):
    files = os.listdir(output_dir)

    # Detect .h5ad and .h5ad.gz files
    exclude_keywords = ['meta', 'process', 'normalized', 'annotations', 'celltype', 'adt']
    h5ad_files = sorted([f for f in files if f.endswith('.h5ad')
                         and not any(exclude in f.lower() for exclude in exclude_keywords)])
    h5ad_gz_files = sorted([f for f in files if f.endswith('.h5ad.gz')
                            and not any(exclude in f.lower() for exclude in exclude_keywords)])
    
    if not (h5ad_files or h5ad_gz_files):
        return None
    
    adata_list = []
    
    # Read .h5ad files
    for h5ad_file in h5ad_files:
        file_path = os.path.join(output_dir, h5ad_file)
        print(f"Reading {h5ad_file}...")

        try:
            adata = sc.read_h5ad(file_path)
            adata = check_and_transpose_adata(adata)
            adata = convert_symbol_to_ensembl(adata)
            adata_list.append(adata)
        except Exception as e:
            print(f"Error reading {h5ad_file}: {e}")
            continue

    # Read .h5ad.gz files
    for h5ad_gz_file in h5ad_gz_files:
        file_path = os.path.join(output_dir, h5ad_gz_file)
        print(f"Reading {h5ad_gz_file}...")

        # Create a decompressed file path (same directory, but without .gz extension)
        decompressed_file_path = file_path[:-3]
        
        try:
            # Decompress the .h5ad.gz file to .h5ad
            with gzip.open(file_path, 'rb') as f_in:
                with open(decompressed_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Read the decompressed .h5ad file
            adata = sc.read_h5ad(decompressed_file_path)
            adata = check_and_transpose_adata(adata)
            adata = convert_symbol_to_ensembl(adata)
            adata_list.append(adata)

            os.remove(decompressed_file_path)

        except Exception as e:
            print(f"Error reading {h5ad_gz_file}: {e}")
            continue

    # merge adata
    adata_list = filter_adata_by_species(adata_list)
    adata_list = filter_adata_by_gene_count(adata_list)

    # Merge all AnnData objects
    if len(adata_list) > 1:
        print("Merging AnnData objects...")
        adata_merged = ad.concat(adata_list, label="batch")
    else:
        adata_merged = adata_list[0]
    
    return adata_merged


# read and merge loom files
def read_and_merge_loom_files(output_dir):
    files = os.listdir(output_dir)

    # Detect .h5ad and .h5ad.gz files
    exclude_keywords = ['meta', 'process', 'normalized', 'annotations', 'celltype', 'adt']
    loom_files = sorted([f for f in files if f.endswith('.loom')
                         and not any(exclude in f.lower() for exclude in exclude_keywords)])
    loom_gz_files = sorted([f for f in files if f.endswith('.loom.gz')
                            and not any(exclude in f.lower() for exclude in exclude_keywords)])
    
    if not (loom_files or loom_gz_files):
        return None
    
    adata_list = []
    
    for loom_file in loom_files:
        file_path = os.path.join(output_dir, loom_file)
        print(f"Reading {loom_file}...")

        try:
            adata = sc.read_loom(file_path, var_names='var_names')
            adata = check_and_transpose_adata(adata)
            adata = convert_symbol_to_ensembl(adata)
            adata_list.append(adata)
        except Exception as e:
            print(f"Error reading {loom_file}: {e}")
            continue

    for loom_gz_file in loom_gz_files:
        file_path = os.path.join(output_dir, loom_gz_file)
        print(f"Reading {loom_gz_file}...")

        decompressed_file_path = file_path[:-3]
        
        try:
            with gzip.open(file_path, 'rb') as f_in:
                with open(decompressed_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            adata = sc.read_loom(decompressed_file_path, var_names='var_names')
            adata = check_and_transpose_adata(adata)
            adata = convert_symbol_to_ensembl(adata)
            adata_list.append(adata)
            os.remove(decompressed_file_path)

        except Exception as e:
            print(f"Error reading {loom_gz_file}: {e}")
            continue

    if len(adata_list) > 1:
        print("Merging AnnData objects...")
        adata_merged = ad.concat(adata_list, label="batch")
    else:
        adata_merged = adata_list[0]
    
    return adata_merged


# read and merge txt/csv/tsv gz files
def read_and_merge_txt_gz_files(output_dir):

    files = os.listdir(output_dir)
    keywords = ['count', 'raw', 'data', 'umi', 'smart', 'express']
    exclude_keywords = ['meta', 'process', 'normalized', 'annotations', 'celltype', 'adt', 'cells', 'tissue']
    txt_gz_files = [f for f in files if any(keyword in f.lower() for keyword in keywords) 
                and not any(exclude in f.lower() for exclude in exclude_keywords) 
                and (f.endswith('.txt.gz') or f.endswith('.csv.gz') or f.endswith('.tsv.gz')
                     or f.endswith('.txt') or f.endswith('.csv') or f.endswith('.tsv'))]
    
    if not txt_gz_files:
        return None
    
    adata_list = []
    
    for txt_gz_file in txt_gz_files:
        file_path = os.path.join(output_dir, txt_gz_file)
        print(f"Reading {txt_gz_file}...")

        if txt_gz_file.endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                delimiter = '\t' if txt_gz_file.endswith(('.txt.gz', '.tsv.gz')) else ','
                adata = sc.read_text(f, delimiter=delimiter)
        
        else:
            with open(file_path, 'r') as f:
                delimiter = '\t' if txt_gz_file.endswith(('.txt', '.tsv')) else ','
                adata = sc.read_text(f, delimiter=delimiter)
        
        adata = check_and_transpose_adata(adata)
        adata = convert_symbol_to_ensembl(adata)
        adata_list.append(adata)
    
    # merge adata
    adata_list = filter_adata_by_species(adata_list)
    adata_list = filter_adata_by_gene_count(adata_list)

    if len(adata_list) > 1:
        print("Merging AnnData objects...")
        adata_merged = ad.concat(adata_list, label="batch")
    else:
        adata_merged = adata_list[0]
    
    return adata_merged


# auto detect and process files
def auto_detect_and_process_files(output_dir):
    
    # file list
    files = os.listdir(output_dir)
    dirs = [d for d in files if os.path.isdir(os.path.join(output_dir, d))]
    
    # define keywords
    keywords = ['count', 'raw', 'data', 'umi', 'smart', 'express']
    exclude_keywords = ['meta', 'process', 'normalized', 'annotations', 'celltype', 'ADT', 'cells', 'tissue']
    
    # detect file types
    mtx_files = [f for f in files if f.endswith('.mtx') or f.endswith('.mtx.gz')]
    h5_files = sorted([f for f in files if f.endswith('.h5')
                       and not any(exclude in f.lower() for exclude in exclude_keywords)])
    h5ad_files = sorted([f for f in files if f.endswith('.h5ad')
                         and not any(exclude in f.lower() for exclude in exclude_keywords)])
    h5ad_gz_files = sorted([f for f in files if f.endswith('.h5ad.gz')
                            and not any(exclude in f.lower() for exclude in exclude_keywords)])
    loom_files = sorted([f for f in files if f.endswith('.loom')
                         and not any(exclude in f.lower() for exclude in exclude_keywords)])
    loom_gz_files = sorted([f for f in files if f.endswith('.loom.gz')
                            and not any(exclude in f.lower() for exclude in exclude_keywords)])
    txt_gz_files = [f for f in files if any(keyword in f.lower() for keyword in keywords) 
                and not any(exclude in f.lower() for exclude in exclude_keywords) 
                and (f.endswith('.txt.gz') or f.endswith('.csv.gz') or f.endswith('.tsv.gz')
                     or f.endswith('.txt') or f.endswith('.csv') or f.endswith('.tsv'))]
    
    adata_list_combine = []

    # process files according to file types
    if mtx_files or dirs:
        print("MTX files detected, processing with read_and_merge_mtx_files...")
        adata_mtx = read_and_merge_mtx_files(output_dir)
        if adata_mtx is not None:
            adata_list_combine.append(adata_mtx)

    if h5_files:
        print("H5 files detected, processing with read_and_merge_h5_files...")
        adata_h5 = read_and_merge_h5_files(output_dir)
        if adata_h5 is not None:
            adata_list_combine.append(adata_h5)

    if h5ad_files or h5ad_gz_files:
        print("H5AD files detected, processing with read_and_merge_h5ad_files...")
        adata_h5ad = read_and_merge_h5ad_files(output_dir)
        if adata_h5ad is not None:
            adata_list_combine.append(adata_h5ad)

    if loom_files or loom_gz_files:
        print("Loom files detected, processing with read_and_merge_loom_files...")
        adata_loom = read_and_merge_loom_files(output_dir)
        if adata_loom is not None:
            adata_list_combine.append(adata_loom)

    if txt_gz_files:
        print("TXT.GZ files detected, processing with read_and_merge_txt_gz_files...")
        adata_txt_gz = read_and_merge_txt_gz_files(output_dir)
        if adata_txt_gz is not None:
            adata_list_combine.append(adata_txt_gz)

    # merge adata
    if len(adata_list_combine) > 1:
        print("Merging all detected AnnData objects...")
        adata_merged_combine = ad.concat(adata_list_combine, label="batch")
    elif len(adata_list_combine) == 1:
        adata_merged_combine = adata_list_combine[0]

    return adata_merged_combine


# save adata
def save_adata(adata):
    
    save_dir = './process'
    save_path = os.path.join(save_dir, 'merged_adata.h5ad')

    # create folders
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # save adata
    adata.write(save_path)


if __name__ == "__main__":
    
    # get arguments
    args = parse_args()
    
    # define source_id and output_dir
    geo_id = args.source_id
    output_dir = args.output_dir

    # download and extract supplementary files
    download_geo_supp_files(geo_id, output_dir)
    check_and_extract_raw_tar(output_dir)
    check_and_extract_tar_gz(output_dir)

    # process adata
    adata = auto_detect_and_process_files(output_dir)
    save_adata(adata)

    # delete folders
    shutil.rmtree(output_dir)