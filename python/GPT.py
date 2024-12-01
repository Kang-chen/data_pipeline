# -*- coding: utf-8 -*-

import os
from openai import OpenAI
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from config import OPEN_API_KEY
from pathlib import Path
import re
import logging


"""
Function: GPT3_turbo
This function uses the GPT-3.5-turbo model from OpenAI.

Input:
    content: string
    A string that represents the user's message to the model. 
    
Output:
    ChatCompletion
    A ChatCompletion object that contains the model's response to the user's message. 

"""
def GPT3_turbo(content):
    client = OpenAI(
        api_key=OPEN_API_KEY,
    )
    
    logging.basicConfig(filename='query_GPT.log', level=logging.INFO)
    logging.info("Query Content:\n" + content)  # Log the content to the file

    chat_completion = client.chat.completions.create(
        model="gpt-4o", #gpt-4o  gpt-3.5-turbo
        messages=[
            {
                "role": "system",
                "content": "You are a single cell bioinformatics specialist."
            },
            {
                "role": "user",
                "content": content
            },
        ],
    )
    return chat_completion


"""
Function: convert_to_dict
This function takes a list of strings in a specific format and converts it into a dictionary.

Input:
output: list
A list of strings where each string is in the format '- key -> value'.

Output:
dict
A dictionary where each key-value pair corresponds to the 'key -> value' in the input list.

"""
def convert_to_dict(output):
    output_dict = {}
    lines = output.split('\n')
    for line in lines:
        if '->' in line:
            key, value = line.split(' -> ')
            output_dict[key.strip()] = None if value.strip() == 'None' else value.strip()
    return output_dict


"""
Function: check_exists
This function Checks if the matched results in input_column, sometimes GPT output the wrong matching results
"""
def check_exists(output, input_column):
    flag = 0
    for value in output.values():
        if value not in input_column and value != None:
            flag = 1
    return flag

"""
Function: blur_search
This function performs a fuzzy search on a list of target columns against a list of observed columns. It uses the fuzzy matching algorithm to find the best match for each target column in the observed columns. If the match score is above a certain threshold (80 in this case), the match is considered valid and added to the result. If no valid match is found, None is added to the result.

Input:
obs_columns: list
A list of strings where each string is an observed column name.

Output:
dict
A dictionary where the keys are the target column names and the values are the best match from the observed columns. If no valid match is found, the value is set to None.
"""
def blur_search(obs_columns):
    target_columns = [
    "dataset_id", "assay", "cell_type_original", 
    "development_stage_original", "disease_original", 
    "donor_id", "sex", "tissue_original", "is_primary"]
    
    matched_columns = {}
    for target in target_columns:
        match, score = process.extractOne(target, obs_columns, scorer=fuzz.partial_ratio)
        if score > 80:
            matched_columns[target] = match
        else:
            matched_columns[target] = None
    return matched_columns

"""
Function: GPT_for_column
This function uses the GPT-3 model to generate text based on an input column.

Input:
    input_column: list
    A list of strings where each string is a column name. 
    Example: ['development_stage', 'disease_original', 'disease_ontology', 'donor_id', 'sex', 'tissue_ontology', 'tissue_original', 'is_primary']

Output:
    dict
    A dictionary where the keys are the input column names and the values are the corresponding outputs from the GPT-3 model. 
    If the model does not generate an output for a column name, the value is set to None.
    Example: {'development_stage': None, 'disease_original': None, 'disease_ontology': None, 'donor_id': None, 'sex': None, 'tissue_ontology': None, 'tissue_original': None, 'is_primary': 'is_primary1'}

"""
def GPT_for_column(input_column):
    # 1. Read the prompt for column match
    script_dir = Path(__file__).parent
    file_name = 'prompt_column_match_2.txt'
    file_path = script_dir / file_name
    with file_path.open('r') as f:
        prompt = f.readline().strip('\n')
    
    # 2. Obtain the full input info by add the input_column information
    content = prompt + "input_column: " + str(input_column)
    
    # 3. Ask the GPT for the answer
    chat_completion = GPT3_turbo(content)
    
    # 4. Obtain the GPT output string
    chat_completion_dict = chat_completion.to_dict()
    output = chat_completion_dict['choices'][0]['message']['content']
    
    # 5. Convert the info to dict
    try:
        output = convert_to_dict(output)
        print('GPT column match finish')
    except:
        output = blur_search(input_column)
        print('Fuzzy column match finish')
    
    # 6. Check the results in GPT3.5
    # Check if the matched results in input_column, sometimes GPT output the wrong matching results
    flag = check_exists(output, input_column)
    if flag == 1:
        output = blur_search(input_column)
        print('Wrong GPT results, Fuzzy column match finish')
    
    return output



def GPT_for_ontology(final_result, name):
    # 1. Define the prompt for ontology match
    prompt = f"Please determine which ontology is most likely related to {name} based on the descriptions provided below. Only return the corresponding ontology_id, without any other information. If nothing match return none :"

    # 2. Convert the entire final_result DataFrame into a string (without the index)
    result_str = final_result.to_string(index=False)  # Convert DataFrame to string without index
    
    # 3. Combine the prompt with the final_result data
    content = f"{prompt} {result_str}"
    
    # 4. Call GPT API to get the response
    chat_completion = GPT3_turbo(content)
    
    # 5. Get the GPT output
    chat_completion_dict = chat_completion.to_dict()
    output = chat_completion_dict['choices'][0]['message']['content'].strip()  # Remove any surrounding whitespace
    
    # 6. Define the regex pattern to find ontology IDs (e.g., MONDO:0005155)
    ontology_id_pattern = r'[A-Za-z]+:\d{7}'  # Matches formats like MONDO:0005155
    
    # 7. Search for a matching ontology ID in the output
    matches = re.findall(ontology_id_pattern, output)
    
    if matches:
        # If matches are found, return the first valid ontology ID
        print(f"Ontology ID found: {matches[0]}")
        return matches[0]
    else:
        # If no valid ontology ID is found, print an error message and return None
        print(f"Invalid output, no valid ontology ID found. Output: {output}")
        return None

''' For test
if __name__ == '__main__':
    input_column = ['orig.ident',
      'nCount_RNA',
      'nFeature_RNA',
      'is_primary1',
      'n_genes',
      'n_genes_by_counts',
      'total_counts',
      'total_counts_mt',
      'pct_counts_mt',
      'doublet_scores',
      'predicted_doublets',
      'n_features',
      'sample_cell_count',
      'is_primary2',
      'is_primary']
    output = GPT_for_column(input_column)
    print(output)
'''















