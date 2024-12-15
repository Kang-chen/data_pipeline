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
import shutil
from GPT import *
import warnings

# Ignore the warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Observation names are not unique")
warnings.filterwarnings("ignore", category=UserWarning, message="Variable names are not unique")
warnings.filterwarnings("ignore", category=UserWarning, message="Trying to modify attribute `.obs` of view, initializing view as actual")

class DataProcessor:
    def __init__(self, adata_path,source_id):
        self.adata = sc.read_h5ad(adata_path)
        self.obs_df = self.adata.obs.copy()
        self.description = source_id
    
    def fuzzy_match_columns(self, obs_columns, required_columns):
        matched_columns = GPT_for_column(obs_columns)
        print("Match:")
        for req_col, matched_col in matched_columns.items():
            print(f"{req_col}: {matched_col}")
            if matched_col is not None:
                self.obs_df[req_col] = self.obs_df[matched_col]
            else:
                self.obs_df[req_col] = None
        self.adata.obs = self.obs_df

    def map_ontology(self, column_name, ontology_class, original_col, mapped_col, ontology_id_col):
        bionty = ontology_class
        name_mapper = {}
        ontology_id_mapper = {}

        self.adata.obs[original_col] = self.adata.obs[original_col].fillna("unknown")

        for name in self.adata.obs[original_col].unique():
            if name is not None and isinstance(name, str) and name.strip():
                record = bionty.from_source(name=name)
                if record is not None:
                   name_mapper[name] = record.name
                   ontology_id_mapper[name] = str(record.ontology_id) 
            else:
                for search_func in [
                    lambda: bionty.filter(name=name).df(),
                    lambda: bionty.search(name, field="synonyms", limit=3).df() if name else None,
                    lambda: bionty.public().search(name, limit=3) if name else None
                ]:
                    search_result = search_func()
                    if search_result is not None and hasattr(search_result, 'empty') and not search_result.empty:
                        ontology_id = search_result.iloc[0].ontology_id
                        record = ontology_class.from_source(ontology_id=ontology_id)
                        name_mapper[name] = record.name
                        ontology_id_mapper[name] = str(ontology_id)
                        record.save()
                        try:
                            record.add_synonym(name)
                        except ValueError as e:
                            if "synonym that is already associated with a record" in str(e):
                                print(e)
                            else:
                                raise
                        break  # Break out of the `for search_func` loop if successful
                    # if search_result is not None and hasattr(search_result, 'empty') and not search_result.empty:
                    #    ontology_id = search_result.iloc[0].ontology_id
                    #    record = ontology_class.from_source(ontology_id=ontology_id)
                    #    name_mapper[name] = record.name
                    #    ontology_id_mapper[name] = ontology_id
                    #    record.save()
                    #     try:
                    #         record.add_synonym(name)
                    #     except ValueError as e:
                    #         if "synonym that is already associated with a record" in str(e):
                    #             print(e)
                    #         else:
                    #             raise
                    #     break
                    # # if search_result and not search_result.empty:
                    # #     ontology_id = search_result.iloc[0].ontology_id
                    # #     record = ontology_class.from_source(ontology_id=ontology_id)
                    # #     name_mapper[name] = record.name
                    # #     ontology_id_mapper[name] = ontology_id
                    # #     record.save()
                        
                else:
                    name_mapper[name] = "unknown"
                    ontology_id_mapper[name] = "unknown"
        #                 try:
        #                     record.add_synonym(name)
        #                 except ValueError as e:
        #                     if "synonym that is already associated with a record" not in str(e):
        #                         raise
        #                 break
        # for name in self.adata.obs[original_col].unique():
        #     if name is not None and isinstance(name, str) and name.strip():
        #         record = bionty.from_source(name=name)
        #     if record is not None:
        #         name_mapper[name] = record.name
        #         ontology_id_mapper[name] = record.ontology_id
        #     else:
        #         for search_func in [
        #             lambda: bionty.filter(name=name).df(),
        #             lambda: bionty.search(name, field="synonyms", limit=3).df(),
        #             lambda: bionty.public().search(name, limit=3)
        #         ]:
        #             search_result = search_func()
        #             if not search_result.empty:
        #                 ontology_id = search_result.iloc[0].ontology_id
        #                 record = ontology_class.from_source(ontology_id=ontology_id)
        #                 name_mapper[name] = record.name
        #                 ontology_id_mapper[name] = ontology_id
        #                 record.save()
        #                 try:
        #                     record.add_synonym(name)
        #                 except ValueError as e:
        #                     if "synonym that is already associated with a record" in str(e):
        #                         # Skip this step if the specific ValueError occurs
        #                         print(e)
        #                     else:
        #                         # Re-raise the exception if it's not the specific one we're checking for
        #                         raise
        #                 break
        #         else:
        #             name_mapper[name] = "unknown"
        #             ontology_id_mapper[name] = "unknown"

            #         if not search_result.empty:
            #             ontology_id = search_result.iloc[0].ontology_id
            #             record = ontology_class.from_source(ontology_id=ontology_id)
            #             name_mapper[name] = record.name
            #             ontology_id_mapper[name] = ontology_id
            #             record.save()
            #             try:
            #                 record.add_synonym(name)
            #             except ValueError as e:
            #                 if "synonym that is already associated with a record" not in str(e):
            #                     raise
            #             break
            #     else:
            #         name_mapper[name] = "unknown"
            #         ontology_id_mapper[name] = "unknown"
            #        else:
            # name_mapper[name] = "unknown"
            # ontology_id_mapper[name] = "unknown"

        self.adata.obs[mapped_col] = self.adata.obs[original_col].map(name_mapper)
        self.adata.obs[ontology_id_col] = self.adata.obs[original_col].map(ontology_id_mapper)

    # def map_ontology(self, column_name, ontology_class, original_col, mapped_col, ontology_id_col):
    #     # bionty = ontology_class.public()
    #     bionty = ontology_class
    #     name_mapper = {}
    #     ontology_id_mapper = {}

    #     for name in self.adata.obs[original_col].unique():
    #         record = bionty.from_source(name = name)
    #         if record is not None:
    #                 name_mapper[name] = record.name
    #                 ontology_id_mapper[name] = record.ontology_id
    #                 continue
    #         if name is not None and isinstance(name, str) and name.strip():
    #             for search_func in [
    #             lambda: bionty.filter(name=name).df(),
    #             lambda: bionty.search(name, field="synonyms", limit=3).df(),
    #             lambda: bionty.public().search(name, limit=3)
    #             ]:
    #             search_result = search_func()  # 执行当前的匿名函数
    #             if not search_result.empty:  # 如果返回结果不为空
    #                 ontology_id = search_result.iloc[0].ontology_id
    #                 record = ontology_class.from_source(ontology_id=ontology_id)
    #                 name_mapper[name] = record.name
    #                 ontology_id_mapper[name] = ontology_id
    #                 record.save()
    #                 try:
    #                     record.add_synonym(name)
    #                 except ValueError as e:
    #                     print(f"Skipping add_synonym for {name}. Error: {e}")
    #                 break
    #         else:
    #             # 如果没有找到匹配的记录，将值设置为'unknown'
    #             name_mapper[name] = "unknown"
    #             ontology_id_mapper[name] = "unknown"    break
    #             if not search_result.empty:
    #                 ontology_id = search_result.iloc[0].ontology_id
    #                 # print(name)   
    #                 record = ontology_class.from_source(ontology_id=ontology_id)
    #                 name_mapper[name] = record.name
    #                 ontology_id_mapper[name] = ontology_id
    #                 record.save()
    #                 try:
    #                     record.add_synonym(name)
    #                 except ValueError as e:
    #                     if "synonym that is already associated with a record" in str(e):
    #                         # Skip this step if the specific ValueError occurs
    #                         print(e)
    #                     else:
    #                         # Re-raise the exception if it's not the specific one we're checking for
    #                         raise
    #             else:
    #                 name_mapper[name] = "unknown"
    #                 ontology_id_mapper[name] = "unknown"
    #     else:
    #         name_mapper[name] = "unknown"
    #         ontology_id_mapper[name] = "unknown"

    # else:    
    #     self.adata.obs[mapped_col] = self.adata.obs[original_col].map(name_mapper)
    #     self.adata.obs[ontology_id_col] = self.adata.obs[original_col].map(ontology_id_mapper)
    
    def detect_species(self):
        gene_ids = self.adata.var_names.unique()
        human_count = sum(1 for gene in gene_ids if gene.startswith('ENSG'))
        mouse_count = sum(1 for gene in gene_ids if gene.startswith('ENSMUSG'))

        if human_count > mouse_count:
            return 'human'
        elif mouse_count > human_count:
            return 'mouse'
        else:
            return 'unknown'

    def process_data(self):
        required_columns = [
            "dataset_id", 
            "organisms",
            "assay", 
            "cell_type_original", 
            "development_stage_original", 
            "disease_original",
            "tissue_original",
            "donor_id",
            "sex", 
            "is_primary"
        ]
        self.fuzzy_match_columns(self.adata.obs.columns.tolist(), required_columns)

        # Map various ontologies
        self.map_ontology("assay", bt.ExperimentalFactor, "assay", "assay_ontology", "assay_ontology_id")
        self.map_ontology("cell_type_original", bt.CellType, "cell_type_original", "cell_type_ontology", "cell_type_ontology_id")
        self.map_ontology("development_stage_original", bt.DevelopmentalStage, "development_stage_original", "development_stage_ontology", "development_stage_ontology_id")
        self.map_ontology("disease_original", bt.Disease, "disease_original", "disease_ontology", "disease_ontology_id")
        self.map_ontology("tissue_original", bt.Tissue, "tissue_original", "tissue_ontology", "tissue_ontology_id")

        # Convert all elements of the data box to string type
        self.adata.obs = self.adata.obs.applymap(str)
        diseases_ontology_value, diseases_ontology_id, tissue_ontology_value, tissue_ontology_id, assay_ontology_value, assay_ontology_id = get_standardized_ontologies(self.description)

        # 更新相关列，如果值为 'unknown' 则使用标准化值
        self.adata.obs.loc[self.adata.obs['disease_ontology'] == 'unknown', 'disease_ontology'] = diseases_ontology_value
        self.adata.obs.loc[self.adata.obs['tissue_ontology'] == 'unknown', 'tissue_ontology'] = tissue_ontology_value
        self.adata.obs.loc[self.adata.obs['assay_ontology'] == 'unknown', 'assay_ontology'] = assay_ontology_value

        # 更新 ontology_id 列
        self.adata.obs.loc[self.adata.obs['disease_ontology_id'] == 'unknown', 'disease_ontology_id'] = diseases_ontology_id
        self.adata.obs.loc[self.adata.obs['tissue_ontology_id'] == 'unknown', 'tissue_ontology_id'] = tissue_ontology_id
        self.adata.obs.loc[self.adata.obs['assay_ontology_id'] == 'unknown', 'assay_ontology_id'] = assay_ontology_id
        
        
        # Define categorical variables and their mappings
        categoricals = {
            self.adata.obs.assay_ontology.name: bt.ExperimentalFactor.name,
            self.adata.obs.cell_type_ontology.name: bt.CellType.name,
            # self.adata.obs.cell_type_ontology_id.name: bt.CellType.ontology_id,
            self.adata.obs.development_stage_ontology.name: bt.DevelopmentalStage.name,
            # self.adata.obs.development_stage_ontology_id.name: bt.DevelopmentalStage.ontology_id,
            self.adata.obs.disease_ontology.name: bt.Disease.name,
            # self.adata.obs.disease_ontology_id.name: bt.Disease.ontology_id,
            self.adata.obs.tissue_ontology.name: bt.Tissue.name,
            # self.adata.obs.tissue_ontology_id.name: bt.Tissue.ontology_id,
        }

        # Detect species
        species = self.detect_species()

        # Create a Curate object to guide validation and annotation
        curate = ln.Curate.from_anndata(
            self.adata, 
            var_index=bt.Gene.ensembl_gene_id,
            categoricals=categoricals, 
            organism=species,
        )

        # Add new records to the database from the variable index
        curate.add_new_from_var_index()

        # Add new records to the database from columns
        curate.add_new_from('assay_ontology')
        curate.add_new_from('cell_type_ontology')
        curate.add_new_from('development_stage_ontology')
        curate.add_new_from('disease_ontology')
        curate.add_new_from('tissue_ontology')

        # Validate the data
        curate.validate()

        # Save artifact
        artifact = curate.save_artifact(description=self.description)

        return artifact

def map_ontology_string(ontology_class, original_string):
    if not isinstance(original_string, str):
        raise ValueError(f"Expected 'original_string' to be a string, got {type(original_string)}")

    # Initialize mapping dictionaries
    name_mapper = {}
    ontology_id_mapper = {}

    name = original_string.strip()
    if not name:
        name_mapper[name] = "unknown"
        ontology_id_mapper[name] = "unknown"
        return ontology_class.from_source(name="unknown")
    
    
    record = ontology_class.from_source(name = name)    
    if record is not None:
        name_mapper[name] = record.name
        ontology_id_mapper[name] = record.ontology_id
        return record
    

    # Search functions to try
    search_funcs = [
        lambda: ontology_class.filter(name=name).df(),
        lambda: ontology_class.search(name, field="synonyms", limit=5).df(),
        lambda: ontology_class.public().search(name, limit=5)
    ]

    results = []
    for search_func in search_funcs:
        search_result = search_func()

        # Ensure 'name' column is present
        if 'name' not in search_result.columns:
            if 'name' in search_result.index.name:
                search_result['name'] = search_result.index  # Assign the index to the 'name' column
            else:
                search_result['name'] = pd.NA  # If index is also missing, fill with NaN

        # Standardize columns
        search_result = search_result.rename(columns={'definition': 'description'}).fillna(pd.NA)
        search_result = search_result[['name', 'ontology_id', 'description', 'synonyms']]

        results.append(search_result)

        # If results are found, stop further searching
        if not search_result.empty:
            break

    # Concatenate all results and process
    final_result = pd.concat(results, ignore_index=True)

    if final_result.empty:
        name_mapper[name] = "unknown"
        ontology_id_mapper[name] = "unknown"
        return ontology_class.from_source(name="unknown")

    # Process ontology ID from results
    ontology_id = GPT_for_ontology(final_result, name)
    if ontology_id is None:
        name_mapper[name] = "unknown"
        ontology_id_mapper[name] = "unknown"
        return ontology_class.from_source(name="unknown")

    # Create ontology record
    record = ontology_class.from_source(ontology_id=ontology_id)
    record.save()
    name_mapper[name] = record.name
    ontology_id_mapper[name] = ontology_id

    try:
        record.add_synonym(name)
    except ValueError as e:
        # Print error message instead of raising
        print(f"Skipping add_synonym for {name}. Error: {e}")

    return record


def split_excel_field(field):
    return [elem for sublist in [item.split('/') for item in field] for elem in sublist]


def standardize_ontologies(diseases, tissue_type, library_protocol):
    diseases_ontology = [
        map_ontology_string(bt.Disease, term)
        for term in split_excel_field([diseases])
    ]
    tissue_ontology = [
        map_ontology_string(bt.Tissue, term)
        for term in split_excel_field([tissue_type])
    ]
    assay_ontology = [
        map_ontology_string(bt.ExperimentalFactor, term)
        for term in split_excel_field([library_protocol])
    ]
    return diseases_ontology, tissue_ontology, assay_ontology


def get_standardized_ontologies(source_id):
    # Load the Excel file
    file_path = './Data_collection.xlsx'
    df = pd.read_excel(file_path)
    
    # Find the corresponding row based on the source_id
    row = df[df['source_id'] == source_id].squeeze()
    if row.empty:
        raise ValueError(f"No data found for source_id: {source_id}")

    # Extract the required column values
    diseases = row['diseases']
    tissue_type = row['tissue_type']
    library_protocol = row['library_protocol']

    # Standardize the ontologies (returns both names and IDs)
    diseases_ontology_list, tissue_ontology_list, assay_ontology_list = standardize_ontologies(
        diseases, tissue_type, library_protocol
    )

    # Get the standardized names and IDs (assuming single values for simplicity)
    diseases_ontology_value = diseases_ontology_list[0].name if diseases_ontology_list else 'unknown'
    diseases_ontology_id = diseases_ontology_list[0].ontology_id if diseases_ontology_list else 'unknown'

    tissue_ontology_value = tissue_ontology_list[0].name if tissue_ontology_list else 'unknown'
    tissue_ontology_id = tissue_ontology_list[0].ontology_id if tissue_ontology_list else 'unknown'

    assay_ontology_value = assay_ontology_list[0].name if assay_ontology_list else 'unknown'
    assay_ontology_id = assay_ontology_list[0].ontology_id if assay_ontology_list else 'unknown'

    return diseases_ontology_value, diseases_ontology_id, tissue_ontology_value, tissue_ontology_id, assay_ontology_value, assay_ontology_id



def load_meta_from_init_excel(artifact, source_id):
    
    # Loading Excel file
    file_path = './Data_collection.xlsx'
    df = pd.read_excel(file_path)
            
    # Find the corresponding row based on the source_id
    row = df[df['source_id'] == source_id].squeeze()
        
    if row.empty:
        raise ValueError(f"No data found for source_id: {source_id}")


    # Add information to the Feature
    ln.Feature(name='data_type', dtype='cat[ULabel]').save()
    ln.Feature(name='species', dtype='cat[ULabel]').save()
    ln.Feature(name='diseases', dtype='cat[ULabel]').save()
    ln.Feature(name='tissue_type', dtype='cat[ULabel]').save()
    ln.Feature(name='has_disease_status', dtype='bool').save()
    ln.Feature(name='has_cell_type', dtype='bool').save()
    ln.Feature(name='has_gender', dtype='bool').save()
    ln.Feature(name='has_age', dtype='bool').save()
    ln.Feature(name='has_raw_data', dtype='cat[ULabel]').save()
    ln.Feature(name='library_protocol', dtype='cat[ULabel]').save()
    ln.Feature(name='pubmed_id', dtype='cat[ULabel]').save()
    ln.Feature(name='publication_title', dtype='cat[ULabel]').save()
    ln.Feature(name='disease_ontology', dtype='cat[bionty.Disease]').save()
    ln.Feature(name='tissue_ontology', dtype='cat[bionty.Tissue]').save()
    ln.Feature(name='assay_ontology', dtype='cat[bionty.ExperimentalFactor]').save()

    # Extract the required column values
    data_type = row['data_type']
    species = row['species']
    diseases = row['diseases']
    tissue_type = row['tissue_type']
    library_protocol = row['library_protocol']
    has_raw_data = str(row['has_raw_data'])
    pubmed_id = str(row['pubmed_id']) 
    publication_title = row['publication_title']

    # 转换为ULabel对象并保存
    data_type_label = ln.ULabel.from_values([data_type], create=True)
    species_label = ln.ULabel.from_values([species], create=True)
    diseases_label = ln.ULabel.from_values([diseases], create=True)
    tissue_type_label = ln.ULabel.from_values([tissue_type], create=True)
    library_protocol_label = ln.ULabel.from_values([library_protocol], create=True)
    has_raw_data_label = ln.ULabel.from_values([has_raw_data], create=True)
    pubmed_id_label = ln.ULabel.from_values([pubmed_id], create=True)
    publication_title_label = ln.ULabel.from_values([publication_title], create=True)

    diseases_ontology, tissue_ontology, assay_ontology = standardize_ontologies(
        diseases, tissue_type, library_protocol
    )

    # 保存这些值
    ln.save(data_type_label)
    ln.save(species_label)
    ln.save(diseases_label)
    ln.save(tissue_type_label)
    ln.save(library_protocol_label)
    ln.save(has_raw_data_label)
    ln.save(pubmed_id_label)
    ln.save(publication_title_label)

    # # 将这些值添加到artifact中
    artifact.features.add_values({
        "data_type": data_type_label,
        "species": species_label,
        "diseases": diseases_label,
        "disease_ontology": diseases_ontology,
        
        "tissue_type": tissue_type_label,
        "tissue_ontology": tissue_ontology,
        
        "has_disease_status": bool(row['has_disease_status']),
        "has_cell_type": bool(row['has_cell_type']),
        "has_gender": bool(row['has_gender']),
        "has_age": bool(row['has_age']),
        "has_raw_data": has_raw_data_label,
        "library_protocol": library_protocol_label,
        "assay_ontology": assay_ontology,
        
        "pubmed_id": pubmed_id_label,
        "publication_title": publication_title_label
    })
    artifact.save()


    return artifact

class LaminDBManager:
    def __init__(self, storage_path):
        self.storage_path = storage_path

    def initialize(self):
        subprocess.run(['lamin', 'init', '--storage', self.storage_path, '--schema', 'bionty'])
        ln.setup.settings.instance._keep_artifacts_local = True

    def close(self):
        subprocess.run(['lamin', 'close'])

    def list_artifacts(self):
        return ln.Artifact.df()

    def upload_artifact(self, artifact):
        artifact.save(upload=True)

    def view_tree(self):
        return ln.setup.settings.storage.root.view_tree()


# Main workflow
def main():

    # Create a parser object
    parser = argparse.ArgumentParser(description='Process some integers.')

    # Add parameter
    parser.add_argument('--source_id', type=str, required=True, help='The source ID to process')

    # Analytic parameter
    args = parser.parse_args()

    # Get parameter values
    source_id = args.source_id

    # Print parameter value
    lamin_db_manager = LaminDBManager('s3://cartabio/ai/data/fujing_test2')
    lamin_db_manager.initialize()
    ln.settings.storage_local = "/home/fujing/fujingge/fujing_test2/"

    processor = DataProcessor('./dataforload/qc_adata.h5ad', source_id)
    artifact = processor.process_data()

    # add descrition
    artifact = load_meta_from_init_excel(artifact, source_id)
    
    # print(lamin_db_manager.list_artifacts())
    lamin_db_manager.upload_artifact(artifact)
    lamin_db_manager.view_tree()

    lamin_db_manager.close()

    # delete folders
    shutil.rmtree('./dataforload')


if __name__ == "__main__":
    main()
