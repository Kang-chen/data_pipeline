#!/usr/bin/env python
# coding: utf-8

# In[1]:
# import packages
import lamindb as ln
import anndata as ad
import bionty as bt
import scanpy as sc
import pandas as pd
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import subprocess

subprocess.run(['lamin', 'init', '--storage', 's3://cartabio/ai/data/fujing_test', '--schema', 'bionty'])
ln.setup.settings.instance._keep_artifacts_local = True
ln.settings.storage_local = "./lamindb/"


# In[2]:
# load adata
# adata = sc.read_h5ad('./testdata/lamindb.pbmc3k.h5ad')
adata = sc.read_h5ad('./dataforload/kang_processing.h5ad')
adata


# In[3]:
# 定义要添加的列名
required_columns = [
    "dataset_id", "assay", 
    "cell_type_original",  "cell_type_ontology", "cell_type_ontology_id",
    "development_stage_original", "development_stage_ontology", "development_stage_ontology_id", 
    "disease_original", "disease_ontology", "disease_ontology_id", 
    "tissue_original", "tissue_ontology", "tissue_ontology_id",
    "donor_id", "sex", "is_primary"
]


# In[5]:
# 初始化一个空的字典来存储匹配结果
obs_df = adata.obs.copy()
matched_columns = {col: None for col in required_columns}

# 遍历现有的列并进行模糊匹配
for col in obs_df.columns:
    for req_col in required_columns:
        if req_col in col:
            matched_columns[req_col] = col

# 打印匹配结果
print("Match:")
for req_col, matched_col in matched_columns.items():
    print(f"{req_col}: {matched_col}")

# 更新obs数据框
for req_col, matched_col in matched_columns.items():
    if matched_col is not None:
        obs_df[req_col] = obs_df[matched_col]
    else:
        obs_df[req_col] = 'unknown'

# 将更新后的obs数据框更新到adata中
adata.obs = obs_df


# In[6]:
# adata.obs['cell_type_original'] = adata.obs['cell_type']


# In[7]:
# 检查更新后的obs数据框
adata.obs.head()


# In[8]:
# mapping celltype
bionty = bt.CellType.public()  # access the public ontology through bionty

name_mapper = {}
ontology_id_mapper = {}

for name in adata.obs['cell_type_original'].unique():
    # search the public ontology and use the ontology id of the top match
    search_result = bionty.search(name)
    if not search_result.empty:
        ontology_id = search_result.iloc[0].ontology_id
        # create a record by loading the top match from bionty
        record = bt.CellType.from_public(ontology_id=ontology_id)
        name_mapper[name] = record.name  # map the original name to standardized name
        ontology_id_mapper[name] = ontology_id 
        record.save()
        record.add_synonym(name)
    else:
        name_mapper[name] = "Unknown"
        ontology_id_mapper[name] = "Unknown"

# 将映射结果应用到adata.obs中
adata.obs['cell_type_ontology'] = adata.obs['cell_type_original'].map(name_mapper)
adata.obs['cell_type_ontology_id'] = adata.obs['cell_type_original'].map(ontology_id_mapper)

# 检查更新后的obs数据框
adata.obs.head()


# In[9]:
# mapping DevelopmentalStage
bionty = bt.DevelopmentalStage.public()  # access the public ontology through bionty

name_mapper = {}
ontology_id_mapper = {}

for name in adata.obs['development_stage_original'].unique():
    # 确保name不是None或空字符串
    if pd.notna(name):
        # search the public ontology and use the ontology id of the top match
        search_result = bionty.search(name)
        if not search_result.empty:
            ontology_id = search_result.iloc[0].ontology_id
            # create a record by loading the top match from bionty
            record = bt.DevelopmentalStage.from_public(ontology_id=ontology_id)
            name_mapper[name] = record.name  # map the original name to standardized name
            ontology_id_mapper[name] = ontology_id 
            record.save()
            record.add_synonym(name)
        else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"

# 将映射结果应用到adata.obs中
adata.obs['development_stage_ontology'] = adata.obs['development_stage_original'].map(name_mapper)
adata.obs['development_stage_ontology_id'] = adata.obs['development_stage_original'].map(ontology_id_mapper)

# 检查更新后的obs数据框
adata.obs.head()


# In[10]:
# mapping Disease
bionty = bt.Disease.public()  # access the public ontology through bionty

name_mapper = {}
ontology_id_mapper = {}

for name in adata.obs['disease_original'].unique():
    # 确保name不是None或空字符串
    if pd.notna(name):
        # search the public ontology and use the ontology id of the top match
        search_result = bionty.search(name)
        if not search_result.empty:
            ontology_id = search_result.iloc[0].ontology_id
            # create a record by loading the top match from bionty
            record = bt.Disease.from_public(ontology_id=ontology_id)
            name_mapper[name] = record.name  # map the original name to standardized name
            ontology_id_mapper[name] = ontology_id 
            record.save()
            record.add_synonym(name)
        else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"

# 将映射结果应用到adata.obs中
adata.obs['disease_ontology'] = adata.obs['disease_original'].map(name_mapper)
adata.obs['disease_ontology_id'] = adata.obs['disease_original'].map(ontology_id_mapper)

# 检查更新后的obs数据框
adata.obs.head()


# In[11]:
# mapping Tissue
bionty = bt.Tissue.public()  # access the public ontology through bionty

name_mapper = {}
ontology_id_mapper = {}

for name in adata.obs['tissue_original'].unique():
    # 确保name不是None或空字符串
    if pd.notna(name):
        # search the public ontology and use the ontology id of the top match
        search_result = bionty.search(name)
        if not search_result.empty:
            ontology_id = search_result.iloc[0].ontology_id
            # create a record by loading the top match from bionty
            record = bt.Tissue.from_public(ontology_id=ontology_id)
            name_mapper[name] = record.name  # map the original name to standardized name
            ontology_id_mapper[name] = ontology_id 
            record.save()
            record.add_synonym(name)
        else:
            name_mapper[name] = "unknown"
            ontology_id_mapper[name] = "unknown"

# 将映射结果应用到adata.obs中
adata.obs['tissue_ontology'] = adata.obs['tissue_original'].map(name_mapper)
adata.obs['tissue_ontology_id'] = adata.obs['tissue_original'].map(ontology_id_mapper)

# 检查更新后的obs数据框
adata.obs.head()


# In[38]:
# 定义分类变量及其映射
categoricals = {
    # 'dataset_id': 'Dataset.name',
    # 'assay': 'Assay.name',
      adata.obs.cell_type_ontology.name: bt.CellType.name,
      adata.obs.cell_type_ontology_id.name: bt.CellType.ontology_id,
      adata.obs.development_stage_ontology.name: bt.DevelopmentalStage.name,
      adata.obs.development_stage_ontology_id.name: bt.DevelopmentalStage.ontology_id,
      adata.obs.disease_ontology.name: bt.Disease.name,
      adata.obs.disease_ontology_id.name: bt.Disease.ontology_id,
      adata.obs.tissue_ontology.name: bt.Tissue.name,
      adata.obs.tissue_ontology_id.name: bt.Tissue.ontology_id,
    # 'donor_id': 'Donor.name',
    # 'sex': 'Sex.name',
    # adata.obs.is_primary.name: 'boolean'  # is_primary 是一个布尔值，不需要映射到注册表
}


# In[39]:
# create an Curate object to guide validation and annotation
curate = ln.Curate.from_anndata(
    adata, 
   # var_index=bt.Gene.symbol,
    var_index=bt.Gene.ensembl_gene_id,
    categoricals=categoricals, 
   # organism="human",
    organism="mouse",
)


# In[40]:


# Add new records to the database from the variable index
curate.add_new_from_var_index()


# In[41]:


# checks data against the defined criteria
curate.validate()


# In[42]:


# save artifact
artifact = curate.save_artifact(description="kang2")


# In[43]:


# Seed a collection
#collection = ln.Collection(
#    artifact, name="scRNA-seq collection", version="1"
#)
#collection.save()
#collection.describe()
#collection.view_lineage()


# In[44]:


# list the artifact in the S3 bucket
ln.setup.settings.storage.root.view_tree()


# In[45]:


# upload artifact to S3 server
artifact.path
artifact.save(upload=True)


# In[46]:


# 查看目录
ln.Artifact.df()


# In[47]:


# 退出lamindb
subprocess.run(['lamin', 'close'])

