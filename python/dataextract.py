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
import mygene
import statistics
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
            print(f"Downloaded: {file}")

    # shut down ftp
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
        else:
            print(f"{tar_file} is not a valid tar file.")




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
        else:
            print(f"{tar_gz_file} is not a valid tar.gz file.")




# remove header
def check_and_remove_header(file_path, file_type):
    lines = []
    header = None
    barcode_column = 0  # Default to first column
    column_lengths = []
    
    # 根据文件扩展名选择打开文件的方式
    open_func = gzip.open if file_path.endswith('.gz') else open
    
    with gzip.open(file_path, 'rt') as f:
        for i, line in enumerate(f):
            if i == 0:
                first_line = line.strip().split('\t')
                if file_type == 'features':
                    header_keywords = ['gene', 'id', 'x']
                elif file_type == 'barcodes':
                    header_keywords = ['barcode', 'id', 'x']
                
                if any(keyword in first_line[0].lower() for keyword in header_keywords):
                    print(f"Header found in {file_type} file:", first_line)
                    header = first_line
                    continue
                else:
                    print(f"No header found in {file_type} file.")
            
            columns = line.strip().split('\t')
            lines.append(columns)
            
            if len(column_lengths) < len(columns):
                column_lengths.extend([[] for _ in range(len(columns) - len(column_lengths))])
            
            for j, col in enumerate(columns):
                column_lengths[j].append(len(col))

    # Analyze column lengths to identify the likely barcode column
    if file_type == 'barcodes' and len(column_lengths) > 1:
        avg_lengths = [statistics.mean(lengths) for lengths in column_lengths]
        std_lengths = [statistics.stdev(lengths) if len(lengths) > 1 else 0 for lengths in column_lengths]
        
        # Choose the column with the highest average length and lowest standard deviation
        barcode_column = max(range(len(avg_lengths)), 
                             key=lambda i: (avg_lengths[i], -std_lengths[i]))

    # Print file structure information
    print(f"\nFile structure for {file_type}:")
    print(f"Total lines: {len(lines) + (1 if header else 0)}")
    print(f"Header: {header if header else 'None'}")
    print(f"Number of columns: {len(column_lengths)}")
    print(f"Average column lengths: {[round(statistics.mean(lengths), 2) for lengths in column_lengths]}")
    if file_type == 'barcodes' and len(column_lengths) > 1:
        print(f"Identified barcode column: {barcode_column}")
    print(f"First few lines of content:")
    for line in lines[:5]:
        print('\t'.join(line))

    # Write back the file without header, using identified barcode column for barcodes
    with gzip.open(file_path, 'wt') as f_out:
        if file_type == 'barcodes':
            for line in lines:
                f_out.write(f"{line[barcode_column]}\n")
        else:
            for line in lines:
                f_out.write('\t'.join(line) + '\n')

    print(f"File saved {'with barcodes from column ' + str(barcode_column) if file_type == 'barcodes' else 'without header'}: {file_path}")
    
    return header, barcode_column, column_lengths





# check and update features
def check_and_update_features(file_path):
    
    with gzip.open(file_path, 'rt') as f:
        features = pd.read_csv(f, sep='\t', header=None)

    # Determine if there are only one columns
    if features.shape[1] == 1:
        print("Features file has only one column, adding row numbers as first column and 'Gene Expression' as third column...")

        # Add a first and third column
        features[1] = range(1, len(features) + 1)
        features[2] = 'Gene Expression'

        new_order = [1, 0, 2]
        features = features.iloc[:, new_order]


        with gzip.open(file_path, 'wt') as f_out:
            features.to_csv(f_out, sep='\t', index=False, header=False)
        print(f"Updated features file saved: {file_path}")
    
    # Determine if there are only two columns
    if features.shape[1] == 2:
        print("Features file has only two columns, adding third column 'Gene Expression'...")

        # Add a third column with the value 'Gene Expression'
        features[2] = 'Gene Expression'

        # Overwrite source file
        with gzip.open(file_path, 'wt') as f_out:
            features.to_csv(f_out, sep='\t', index=False, header=False)
        print(f"Updated features file saved: {file_path}")
    else:
        print("Features file already has three columns, no changes made.")





# identify species
def identify_species(adata):
    gene_names = adata.var_names
    print(f"Gene names: {gene_names}") 
    # Check if gene_names is empty
    if len(gene_names) == 0:
        raise ValueError("The gene names list is empty. Please check the input data.")

    try:
       is_symbol = not gene_names[0].startswith("ENSG") and not gene_names[0].startswith("ENSMUS")
    except IndexError:
       raise ValueError("Gene names are not accessible. Please check the input data.")

    # is_symbol = not gene_names[0].startswith("ENSG") and not gene_names[0].startswith("ENSMUS")

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
    
    # print AnnData count
    print(f"Human AnnData count: {human_count}")
    print(f"Mouse AnnData count: {mouse_count}")
    
    # Retaining adata according to counts
    if human_count > mouse_count:
        print("Retaining human samples only.")
        return human_adata
    else:
        print("Retaining mouse samples only.")
        return mouse_adata




# convert gene symbols to ensembl
def convert_symbol_to_ensembl(adata):
    # Ensure we're working with a copy, not a view
    adata = adata.copy()
    
    mg = mygene.MyGeneInfo()

    # Function to check if gene names are symbols or Ensembl IDs
    def is_symbol(gene_names):
        gene_prefixes = ["ENSG", "ENSMUS"]
        return not any(gene_names[0].startswith(prefix) for prefix in gene_prefixes)

    # Determine if the gene names are symbols or already Ensembl IDs
    if is_symbol(adata.var.index):
        print("Detected gene names as symbols.")

        # Identify species
        species = identify_species(adata)  # Use your species detection function here
        print(f"Detected species: {species}")

        # Determine query parameters based on species
        species_map = {
            "human": ("symbol", "ensembl.gene", "human", "grch37"),
            "mouse": ("symbol", "ensembl.gene", "mouse", "grcm38")
        }
        if species not in species_map:
            raise ValueError("Species could not be identified. Please provide valid human or mouse data.")

        # Unpack query parameters
        scopes, fields, species_query, build_version = species_map[species]

        # Perform batch query
        query_result = mg.querymany(adata.var.index.tolist(), scopes=scopes, fields=fields, as_dataframe=True, species=species_query, build=build_version)

    else:
        print("Gene names are already Ensembl IDs. No conversion necessary.")
        return adata

    # Reset index, remove duplicates, and fill missing values
    query_result = query_result.reset_index().drop_duplicates(subset='query')
    query_result['ensembl.gene'] = query_result['ensembl.gene'].fillna('unmapped')

    # Create a mapping from gene symbols to Ensembl IDs
    gene_map = query_result.set_index('query')['ensembl.gene'].to_dict()

    # Update adata.var with the mapped Ensembl IDs
    adata.var['ensembl_id'] = adata.var.index.map(gene_map).fillna('unmapped')

    # Count and remove unmapped genes
    unmapped_count = (adata.var['ensembl_id'] == 'unmapped').sum()
    adata = adata[:, adata.var['ensembl_id'] != 'unmapped']

    # Count and remove duplicate Ensembl IDs, keeping the first occurrence
    duplicate_count = adata.var['ensembl_id'].duplicated().sum()
    adata = adata[:, ~adata.var['ensembl_id'].duplicated(keep='first')]

    print(f"Number of unmapped genes removed: {unmapped_count}")
    print(f"Number of duplicate Ensembl IDs removed: {duplicate_count}")

    # Set the Ensembl IDs as the new index
    adata.var_names = adata.var['ensembl_id']
    adata.var = adata.var.set_index('ensembl_id')

    print("Gene symbols converted to Ensembl IDs. Unmapped and duplicate genes removed.")
    print(f"Original number of genes: {len(gene_map)}")
    print(f"Final number of genes: {adata.n_vars}")
    
    return adata




# detect adata_list by gene number
def filter_adata_by_gene_count(adata_list, gene_threshold=3000):

    filtered_adata_list = []
    
    for adata in adata_list:
        gene_count = adata.var.shape[0] 
        
        print(f"AnnData object has {gene_count} genes.")
        
        # save object more than 3000
        if gene_count > gene_threshold:
            filtered_adata_list.append(adata)
        else:
            print(f"AnnData object with {gene_count} genes discarded (less than {gene_threshold} genes).")
    
    # print saved adata numbers
    print(f"Retained {len(filtered_adata_list)} AnnData objects with gene count greater than {gene_threshold}.")
    
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
    var_has_genes = contains_gene_names(adata.var_names)
    obs_has_genes = contains_gene_names(adata.obs_names)

    # If var_names already contain gene names, do nothing
    if var_has_genes:
        print("var_names contains gene names. No action required.")
    # If obs_names contain gene names, transpose the data
    elif obs_has_genes:
        print("obs_names contains gene names. Transposing the data...")
        adata = adata.transpose()
    else:
        print("No gene names detected in var_names or obs_names.")

    return adata




# filter single adata by gene count
def filter_single_adata_by_gene_count(adata, gene_threshold=3000):
    gene_count = adata.var.shape[0]
    if gene_count > gene_threshold:
        return adata
    else:
        return None 




# detect and process metadata
def detect_and_process_metadata(output_dir, adata):
    files = os.listdir(output_dir)
    print(f"Files in directory: {files}")
    
    meta_keywords = ['meta', 'metadata', 'annotation', 'celltype', 'cluster', 'cell', 'readme']
    count_keywords = ['count', 'raw', 'umi', 'smart', 'express', 'matrix']
    meta_files = [f for f in files if any(keyword in f.lower() for keyword in meta_keywords) 
                  and not any(keyword in f.lower() for keyword in count_keywords)
                  and f.lower().endswith(('.txt.gz', '.tsv.gz', '.csv.gz', '.txt', '.tsv', '.csv', '.xlsx'))]
    
    print(f"Detected metadata files: {meta_files}")
    
    processed_files = []
    
    def process_metadata_file(file, file_path):
        try:
            if file.endswith('.xlsx'):
                meta_df = pd.read_excel(file_path)
            elif file.endswith('.gz'):
                meta_df = pd.read_csv(file_path, compression='gzip', sep=None, engine='python')
            else:
                meta_df = pd.read_csv(file_path, sep=None, engine='python')
            
            # Store the original metadata as a DataFrame in uns
            adata.uns[f'original_metadata_{file}'] = meta_df.copy()
            
            potential_id_columns = [col for col in meta_df.columns if col.lower() in ['cell', 'barcode', 'cell_id', 'cell.id', 'cell.name', 'cell_name']]
            
            if potential_id_columns:
                id_column = potential_id_columns[0]
                print(f"Using column '{id_column}' as cell barcodes for {file}")
                meta_df.set_index(id_column, inplace=True)
            elif meta_df.index.name is None and meta_df.index.dtype == 'object':
                print(f"Using row names as cell barcodes for {file}")
                meta_df.index.name = 'cell_barcode'
            else:
                print(f"Warning: No cell ID column found in {file}. Using the index as is.")
            
            meta_df.index = meta_df.index.astype(str)
            
            matching_indices = adata.obs.index.astype(str).isin(meta_df.index)
            matched_cells = matching_indices.sum()
            
            if matched_cells == 0:
                print(f"Warning: No cell barcodes from {file} matched with the AnnData object. No changes made to adata.obs.")
                return
            
            for col in meta_df.columns:
                if col not in adata.obs.columns:
                    # Determine the appropriate dtype for the new column
                    if pd.api.types.is_numeric_dtype(meta_df[col]):
                        adata.obs[col] = pd.Series(dtype=meta_df[col].dtype)
                    else:
                        # For non-numeric types, use object dtype to accommodate various data types
                        adata.obs[col] = pd.Series(dtype='object')
                    
                    # Add the data, preserving the original dtype
                    adata.obs.loc[matching_indices, col] = meta_df.loc[adata.obs.index[matching_indices].astype(str), col].values
            
            print(f"Metadata from {file} has been added to the AnnData object and stored in uns.")
            print(f"Matched {matched_cells} cells out of {adata.n_obs} total cells.")
            processed_files.append(file)
        except Exception as e:
            print(f"Error processing {file}: {str(e)}")
            import traceback
            traceback.print_exc()
    
    if meta_files:
        print("Metadata files detected:")
        for i, file in enumerate(meta_files):
            print(f"{i+1}. {file}")
        
        # Process metadata files
        for file in meta_files:
            file_path = os.path.join(output_dir, file)
            process_metadata_file(file, file_path)
        
        print(f"Metadata processing completed. {len(processed_files)} files were processed.")
    else:
        print("No metadata files detected automatically. Checking for other potential metadata files...")
        other_files = [f for f in files if not any(keyword in f.lower() for keyword in count_keywords) 
                       and f.lower().endswith(('.txt.gz', '.tsv.gz', '.csv.gz', '.txt', '.tsv', '.csv', '.xlsx'))]
        
        if other_files:
            print("Other potential metadata files found:")
            for i, file in enumerate(other_files):
                print(f"{i+1}. {file}")
            
            while True:
                choice = input("Enter the number of the file you want to inspect as potential metadata (or 'done' if finished): ")
                if choice.lower() == 'done':
                    break
                try:
                    file_to_inspect = other_files[int(choice) - 1]
                    file_path = os.path.join(output_dir, file_to_inspect)
                    
                    # Read and display the first few rows
                    if file_to_inspect.endswith('.xlsx'):
                        df = pd.read_excel(file_path, nrows=5)
                    elif file_to_inspect.endswith('.gz'):
                        df = pd.read_csv(file_path, compression='gzip', sep=None, engine='python', nrows=5)
                    else:
                        df = pd.read_csv(file_path, sep=None, engine='python', nrows=5)
                    
                    print(f"\nFirst few rows of {file_to_inspect}:")
                    print(df)
                    
                    use_as_metadata = input("Do you want to use this file as metadata? (yes/no): ")
                    if use_as_metadata.lower() == 'yes':
                        process_metadata_file(file_to_inspect, file_path)
                except Exception as e:
                    print(f"Error processing file: {str(e)}")
        else:
            print("No potential metadata files found.")
    
    # Ensure adata.uns['metadata'] exists even if empty
    if 'metadata' not in adata.uns:
        adata.uns['metadata'] = []
    
    # Log processed files
    adata.uns['processed_metadata_files'] = processed_files
    
    return adata


# get prefixes
def extract_prefix(file_list):
    return set([f.rsplit('_', 1)[0] for f in file_list])

# # read and merge 10x files
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

    # # 调试输出文件列表
    # print("Files found:", files)
    # print("Features files:", features_files)
    # print("Barcodes files:", barcodes_files)
    # print("Matrix files:", matrix_files)
    # get prefixes
    features_prefixes = extract_prefix(features_files)
    barcodes_prefixes = extract_prefix(barcodes_files)
    matrix_prefixes = extract_prefix(matrix_files)
    
  
    prefixes = features_prefixes & barcodes_prefixes & matrix_prefixes
    
    adata_list = []
    
    for prefix in prefixes:
        feature_file = [f for f in features_files if f.startswith(prefix)][0]
        barcode_file = [f for f in barcodes_files if f.startswith(prefix)][0]
        matrix_file = [f for f in matrix_files if f.startswith(prefix)][0]

        # create folders
        folder_name = os.path.join(output_dir, prefix.strip('_'))
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        # move files and rename
        shutil.move(os.path.join(output_dir, feature_file), os.path.join(folder_name, 'features.tsv.gz'))
        shutil.move(os.path.join(output_dir, barcode_file), os.path.join(folder_name, 'barcodes.tsv.gz'))
        shutil.move(os.path.join(output_dir, matrix_file), os.path.join(folder_name, 'matrix.mtx.gz'))
        
        # check features and barcodes
        _, _, _ = check_and_remove_header(os.path.join(folder_name, 'features.tsv.gz'), file_type='features')
        barcodes_header, barcode_column, column_lengths = check_and_remove_header(os.path.join(folder_name, 'barcodes.tsv.gz'), file_type='barcodes')
        check_and_update_features(os.path.join(folder_name, 'features.tsv.gz'))

        # read 10x files
        print(f"Reading {folder_name}...")
        adata = sc.read_10x_mtx(folder_name, var_names='gene_symbols', cache=True)
        
        # If there were additional columns in the barcodes file, add them to adata.obs
        if len(column_lengths) > 1:
            with gzip.open(os.path.join(folder_name, 'barcodes.tsv.gz'), 'rt') as f:
                original_barcodes_data = [line.strip().split('\t') for line in f]
            
            if barcodes_header:
                column_names = barcodes_header
            else:
                column_names = [f'barcode_col_{i}' for i in range(len(column_lengths))]
            
            for i, col_name in enumerate(column_names):
                if i != barcode_column:
                    adata.obs[col_name] = [row[i] if i < len(row) else '' for row in original_barcodes_data]

        adata = filter_single_adata_by_gene_count(adata)
        if adata is None:
            continue
        adata = convert_symbol_to_ensembl(adata)
        adata_list.append(adata)

    # directly read folders
    for folder in dirs:
        folder_path = os.path.join(output_dir, folder)
        print(f"Directly reading from folder {folder_path}...")
        adata = sc.read_10x_mtx(folder_path, var_names='gene_symbols', cache=True)
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
        print("No .h5 files found in the directory.")
        return None
    
    adata_list = []
    
    # read h5 file one by one
    for h5_file in h5_files:
        file_path = os.path.join(output_dir, h5_file)
        print(f"Reading {h5_file}...")

        try:
            adata = sc.read_10x_h5(file_path)
            adata = check_and_transpose_adata(adata)
            adata = convert_symbol_to_ensembl(adata)
            adata_list.append(adata)
        except Exception as e:
            continue   

        #adata = sc.read_10x_h5(file_path)
        #adata = check_and_transpose_adata(adata)
        #adata = convert_symbol_to_ensembl(adata)
        #adata_list.append(adata)*/
    
    # merge adata
    adata_list = filter_adata_by_species(adata_list)
    adata_list = filter_adata_by_gene_count(adata_list)
    #  merge h5 file
    if not adata_list:
        return None
    elif len(adata_list) == 1:
        return adata_list[0]
    else:
        print("Merging AnnData objects...")
        adata_merged = ad.concat(adata_list, label="batch")     
        return adata_merged
    
    #  merge h5 file
    #if len(adata_list) > 1:
        #print("Merging AnnData objects...")
        #adata_merged = ad.concat(adata_list, label="batch")
    #else:
        #adata_merged = adata_list[0]
    
    #return adata_merged




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
        print("No .h5ad or .h5ad.gz files found in the directory.")
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
        print("No .loom or .loom.gz files found in the directory.")
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
        print("No count-related .txt, .csv, .tsv, .txt.gz, .csv.gz, or .tsv.gz files found in the directory.")
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

# read and merge excel files
def read_and_merge_excel_files(output_dir):
    files = os.listdir(output_dir)
    keywords = ['count', 'raw', 'data', 'umi', 'smart', 'express']
    exclude_keywords = ['meta', 'process', 'normalized', 'annotations', 'celltype', 'adt', 'cells', 'tissue']
    excel_files = [f for f in files if any(keyword in f.lower() for keyword in keywords) 
                   and not any(exclude in f.lower() for exclude in exclude_keywords) 
                   and f.endswith(('.xlsx', '.xls'))]
    
    if not excel_files:
        print("No count-related Excel files found in the directory.")
        return None
    
    adata_list = []
    
    for excel_file in excel_files:
        file_path = os.path.join(output_dir, excel_file)
        print(f"Reading {excel_file}...")

        try:
            # Read Excel file
            df = pd.read_excel(file_path, engine='openpyxl' if excel_file.endswith('.xlsx') else 'xlrd')
            
            # Create AnnData object
            adata = sc.AnnData(X=df.values, obs=pd.DataFrame(index=df.index), var=pd.DataFrame(index=df.columns))
            
            adata = check_and_transpose_adata(adata)
            adata = convert_symbol_to_ensembl(adata)
            adata_list.append(adata)
        except Exception as e:
            print(f"Error reading {excel_file}: {str(e)}")
            continue
    
    # merge adata
    adata_list = filter_adata_by_species(adata_list)
    adata_list = filter_adata_by_gene_count(adata_list)

    if len(adata_list) > 1:
        print("Merging AnnData objects...")
        adata_merged = ad.concat(adata_list, label="batch")
    elif len(adata_list) == 1:
        adata_merged = adata_list[0]
    else:
        print("No valid AnnData objects created from Excel files.")
        return None
    
    return adata_merged




# auto detect and process files
def auto_detect_and_process_files(output_dir):
    files = os.listdir(output_dir)
    
    # Detect different file types
    mtx_files = [f for f in files if f.endswith('.mtx') or f.endswith('.mtx.gz')]
    h5_files = [f for f in files if f.endswith('.h5')]
    h5ad_files = [f for f in files if f.endswith('.h5ad')]
    h5ad_gz_files = [f for f in files if f.endswith('.h5ad.gz')]
    loom_files = [f for f in files if f.endswith('.loom')]
    loom_gz_files = [f for f in files if f.endswith('.loom.gz')]
    excel_files = [f for f in files if f.endswith(('.xlsx', '.xls'))]
    
    keywords = ['count', 'raw', 'data', 'umi', 'smart', 'express']
    exclude_keywords = ['meta', 'process', 'normalized', 'annotations', 'celltype', 'adt', 'cells', 'tissue']
    txt_gz_files = [f for f in files if any(keyword in f.lower() for keyword in keywords) 
                and not any(exclude in f.lower() for exclude in exclude_keywords) 
                and (f.endswith('.txt.gz') or f.endswith('.csv.gz') or f.endswith('.tsv.gz')
                     or f.endswith('.txt') or f.endswith('.csv') or f.endswith('.tsv'))]
    
    dirs = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]
    
    adata_list_combine = []

    # Process files according to file types
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
        print("TXT/CSV/TSV files detected, processing with read_and_merge_txt_gz_files...")
        adata_txt_gz = read_and_merge_txt_gz_files(output_dir)
        if adata_txt_gz is not None:
            adata_list_combine.append(adata_txt_gz)

    if excel_files:
        print("Excel files detected, processing with read_and_merge_excel_files...")
        adata_excel = read_and_merge_excel_files(output_dir)
        if adata_excel is not None:
            adata_list_combine.append(adata_excel)

    # Merge adata
    if len(adata_list_combine) > 1:
        print("Merging all detected AnnData objects...")
        adata_merged_combine = ad.concat(adata_list_combine, label="batch")
    elif len(adata_list_combine) == 1:
        adata_merged_combine = adata_list_combine[0]
    else:
        print("No valid data found.")
        return None

    # Process metadata
    adata_merged_combine = detect_and_process_metadata(output_dir, adata_merged_combine)

    return adata_merged_combine



# save adata
def save_adata(adata):
    save_dir = './process'
    save_path = os.path.join(save_dir, 'merged_adata.h5ad')

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        print(f"Directory '{save_dir}' created.")

    try:
        # Convert non-string objects to strings in all DataFrames in adata.uns
        for key, value in adata.uns.items():
            if isinstance(value, pd.DataFrame):
                adata.uns[key] = value.astype(str)
        
        # Ensure all columns in adata.obs are compatible with HDF5 storage
        for col in adata.obs.columns:
            if adata.obs[col].dtype == 'object':
                adata.obs[col] = adata.obs[col].astype(str)
        
        adata.write(save_path)
        print(f"AnnData object saved to '{save_path}'.")
    except Exception as e:
        print(f"Error saving AnnData object: {str(e)}")
        print("Detailed error information:")
        import traceback
        traceback.print_exc()

    # Print the structure of uns
    print("\nStructure of uns:")
    for key, value in adata.uns.items():
        print(f"{key}: {type(value)}")
        if isinstance(value, dict):
            print(f"  First few keys: {list(value.keys())[:5]}")
        elif isinstance(value, pd.DataFrame):
            print(f"  DataFrame shape: {value.shape}")
            print(f"  DataFrame columns: {value.columns.tolist()}")
            print(f"  DataFrame dtypes:\n{value.dtypes}")




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




adata




adata.var.head()





adata.obs.head()




# Print the keys in uns
print("Keys in adata.uns:")
print(adata.uns.keys())

# Print the contents of each key in uns
for key in adata.uns.keys():
    print(f"\nContents of adata.uns['{key}']:")
    print(adata.uns[key])

# If you want to see the structure and types of the uns contents:
def print_uns_structure(uns, indent=0):
    for key, value in uns.items():
        print('  ' * indent + f"{key}: ", end='')
        if isinstance(value, dict):
            print()
            print_uns_structure(value, indent + 1)
        elif isinstance(value, list):
            print(f"List of {len(value)} items")
            if value:
                print('  ' * (indent + 1) + f"First item type: {type(value[0])}")
                if isinstance(value[0], dict):
                    print_uns_structure(value[0], indent + 2)
        else:
            print(f"{type(value)} - {str(value)[:50]}{'...' if len(str(value)) > 50 else ''}")

print("\nStructure of adata.uns:")
print_uns_structure(adata.uns)

