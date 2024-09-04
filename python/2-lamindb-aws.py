#!/usr/bin/env python
# coding: utf-8

# import packages
import lamindb as ln
import anndata as ad
import bionty as bt
import scanpy as sc
import pandas as pd
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import subprocess
import argparse
from GPT import *   



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


subprocess.run(['lamin', 'init', '--storage', 's3://cartabio/ai/data/test', '--schema', 'bionty'])
ln.setup.settings.instance._keep_artifacts_local = True
# ln.settings.storage_local = "./lamindb/"


# load adata
# adata = sc.read_h5ad('../test/dataforload/kang_processing.h5ad')
adata = sc.read_h5ad('./dataforload/kang_processing.h5ad')

adata
obs_columns = adata.obs.columns.tolist()


# Defines the column name to add
required_columns = [
    "dataset_id", 
    "assay", 
    "cell_type_original", # "cell_type_ontology", "cell_type_ontology_id",
    "development_stage_original", #"development_stage_ontology", "development_stage_ontology_id", 
    "disease_original",# "disease_ontology", "disease_ontology_id", 
    "tissue_original",# "tissue_ontology", "tissue_ontology_id",
    "donor_id",
    "sex", 
    "is_primary"
]


# Initializes an empty dictionary to store match results
obs_df = adata.obs.copy()

''' fuzzy match
matched_columns = {col: None for col in required_columns}
# Iterate over existing columns and perform fuzzy matching
for col in obs_df.columns:
    for req_col in required_columns:
        if req_col in col:
            matched_columns[req_col] = col
'''

matched_columns = GPT_for_column(obs_columns)


# Print matching results
print("Match:")
for req_col, matched_col in matched_columns.items():
    print(f"{req_col}: {matched_col}")


# Update the obs data frame
for req_col, matched_col in matched_columns.items():
    if matched_col is not None:
        obs_df[req_col] = obs_df[matched_col]
    else:
        obs_df[req_col] = None


# Update the obs data frame to adata
adata.obs = obs_df

# mapping assays
bionty = bt.ExperimentalFactor.public()  # access the public ontology through bionty
name_mapper = {}
for name in adata.obs['assay'].unique():
    if name is not None and isinstance(name, str) and name.strip():  # Check whether name is empty or invalid
        # search the public ontology and use the ontology id of the top match
        search_result = bionty.search(name)
        if not search_result.empty:
            ontology_id = bionty.search(name).iloc[0].ontology_id
            # create a record by loading the top match from bionty
            record = bt.ExperimentalFactor.from_source(ontology_id=ontology_id)
            name_mapper[name] = record.name  # map the original name to standardized name
            record.save()
            record.add_synonym(name)
        else:
            name_mapper[name] = "unknown"
    else:
        name_mapper[name] = "unknown"

# Apply the mapping results to adata.obs
adata.obs['assay'] = adata.obs['assay'].map(name_mapper)


# mapping celltype
bionty = bt.CellType.public()  # access the public ontology through bionty

name_mapper = {}
ontology_id_mapper = {}


for name in adata.obs['cell_type_original'].unique():
    if name is not None and isinstance(name, str) and name.strip():  # Check whether name is empty or invalid
        # search the public ontology and use the ontology id of the top match
        search_result = bionty.search(name)
        if not search_result.empty:
            ontology_id = search_result.iloc[0].ontology_id
            # create a record by loading the top match from bionty
            record = bt.CellType.from_source(ontology_id=ontology_id)
            name_mapper[name] = record.name  # map the original name to standardized name
            ontology_id_mapper[name] = ontology_id 
            record.save()
            record.add_synonym(name)
        else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"
    else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"
    

# Apply the mapping results to adata.obs
adata.obs['cell_type_ontology'] = adata.obs['cell_type_original'].map(name_mapper)
adata.obs['cell_type_ontology_id'] = adata.obs['cell_type_original'].map(ontology_id_mapper)


# mapping DevelopmentalStage
bionty = bt.DevelopmentalStage.public()  # access the public ontology through bionty

name_mapper = {}
ontology_id_mapper = {}

for name in adata.obs['development_stage_original'].unique():
    if name is not None and isinstance(name, str) and name.strip():  # Check whether name is empty or invalid
        # search the public ontology and use the ontology id of the top match
        search_result = bionty.search(name)
        if not search_result.empty:
            ontology_id = search_result.iloc[0].ontology_id
            # create a record by loading the top match from bionty
            record = bt.DevelopmentalStage.from_source(ontology_id=ontology_id)
            name_mapper[name] = record.name  # map the original name to standardized name
            ontology_id_mapper[name] = ontology_id 
            record.save()
            record.add_synonym(name)
        else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"
    else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"

# Apply the mapping results to adata.obs
adata.obs['development_stage_ontology'] = adata.obs['development_stage_original'].map(name_mapper)
adata.obs['development_stage_ontology_id'] = adata.obs['development_stage_original'].map(ontology_id_mapper)


# mapping Disease
bionty = bt.Disease.public()  # access the public ontology through bionty

name_mapper = {}
ontology_id_mapper = {}

for name in adata.obs['disease_original'].unique():
    if name is not None and isinstance(name, str) and name.strip():  # Check whether name is empty or invalid
        # search the public ontology and use the ontology id of the top match
        search_result = bionty.search(name)
        if not search_result.empty:
            ontology_id = search_result.iloc[0].ontology_id
            # create a record by loading the top match from bionty
            record = bt.Disease.from_source(ontology_id=ontology_id)
            name_mapper[name] = record.name  # map the original name to standardized name
            ontology_id_mapper[name] = ontology_id 
            record.save()
            record.add_synonym(name)
        else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"
    else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"

# Apply the mapping results to adata.obs
adata.obs['disease_ontology'] = adata.obs['disease_original'].map(name_mapper)
adata.obs['disease_ontology_id'] = adata.obs['disease_original'].map(ontology_id_mapper)


# mapping Tissue
bionty = bt.Tissue.public()  # access the public ontology through bionty

name_mapper = {}
ontology_id_mapper = {}

for name in adata.obs['tissue_original'].unique():
    if name is not None and isinstance(name, str) and name.strip():  # Check whether name is empty or invalid
        # search the public ontology and use the ontology id of the top match
        search_result = bionty.search(name)
        if not search_result.empty:
            ontology_id = search_result.iloc[0].ontology_id
            # create a record by loading the top match from bionty
            record = bt.Tissue.from_source(ontology_id=ontology_id)
            name_mapper[name] = record.name  # map the original name to standardized name
            ontology_id_mapper[name] = ontology_id 
            record.save()
            record.add_synonym(name)
        else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"
    else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"

# Apply the mapping results to adata.obs
adata.obs['tissue_ontology'] = adata.obs['tissue_original'].map(name_mapper)
adata.obs['tissue_ontology_id'] = adata.obs['tissue_original'].map(ontology_id_mapper)

# All elements of the data box are coverted to string type
adata.obs = adata.obs.applymap(str)

# Define categorical variables and their mappings
categoricals = {
      adata.obs.assay.name: bt.ExperimentalFactor.name,
      adata.obs.cell_type_ontology.name: bt.CellType.name,
      adata.obs.cell_type_ontology_id.name: bt.CellType.ontology_id,
      adata.obs.development_stage_ontology.name: bt.DevelopmentalStage.name,
      adata.obs.development_stage_ontology_id.name: bt.DevelopmentalStage.ontology_id,
      adata.obs.disease_ontology.name: bt.Disease.name,
      adata.obs.disease_ontology_id.name: bt.Disease.ontology_id,
      adata.obs.tissue_ontology.name: bt.Tissue.name,
      adata.obs.tissue_ontology_id.name: bt.Tissue.ontology_id,
}


# Identify species
gene_ids = adata.var_names.unique()

# Check the prefixes of Ensembl gene ids
def detect_species(gene_ids):
    human_count = sum(1 for gene in gene_ids if gene.startswith('ENSG'))
    mouse_count = sum(1 for gene in gene_ids if gene.startswith('ENSMUSG'))

    if human_count > mouse_count:
        return 'human'
    elif mouse_count > human_count:
        return 'mouse'
    else:
        return 'unknown'

species = detect_species(gene_ids)


# create an Curate object to guide validation and annotation
curate = ln.Curate.from_anndata(
    adata, 
    var_index=bt.Gene.ensembl_gene_id,
    categoricals=categoricals, 
    organism=species,
)


# Add new records to the database from the variable index
curate.add_new_from_var_index()


# checks data against the defined criteria
curate.validate()


# save artifact
artifact = curate.save_artifact(description=source_id)


# list the artifact in the S3 bucket
ln.setup.settings.storage.root.view_tree()


# upload artifact to S3 server
artifact.path
artifact.save(upload=True)


# View catalog
ln.Artifact.df()

# Exit lamindb
subprocess.run(['lamin', 'close'])

