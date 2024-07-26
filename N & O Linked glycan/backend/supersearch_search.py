import requests
import pandas as pd
import sys
import re
from functools import reduce
import time
import random
import os
import logging
from io import StringIO
from csv import DictReader, DictWriter

# Define URLs and payloads
SEARCH_URL = "https://beta-api.glygen.org/supersearch/search"
LIST_URL = "https://beta-api.glygen.org/supersearch/list/"
GLYCAN_URL = "https://beta-api.glygen.org/glycan/detail/"
OUTPUT_PATH = "./output/supersearch_results.tsv"
REDENDMONO_URL = "https://github.com/glygen-glycan-data/PyGly/raw/master/smw/glycandata/export/redendmono.tsv"

PAYLOAD = {
    "concept_query_list": [
        {
            "concept": "site",
            "query": {
                "aggregator": "$and",
                "aggregated_list": [],
                "unaggregated_list": [
                    {
                        "path": "glycosylation_flag",
                        "order": 0,
                        "operator": "$eq",
                        "string_value": "true"
                    }
                ]
            }
        }
    ]
}

def supersearch_search():
    response = requests.post(url=SEARCH_URL, json=PAYLOAD)
    search_data = response.json()  # converting response to a json file
    list_id = search_data["results_summary"]["site"]["list_id"]  # obtaining the id
    list_res = requests.post(LIST_URL, json={"id":list_id})
    list_data = list_res.json()
    total_length = list_data["pagination"]["total_length"]
    return total_length, list_id

# Function to make API calls with exponential backoff
def fetch_with_backoff(url, payload, max_retries=5):
    retry_delay = 1
    for attempt in range(max_retries):
        try:
            response = requests.post(url, json=payload)
            response.raise_for_status()
            return response.json()
        except requests.RequestException:
            time.sleep(retry_delay)
            retry_delay *= 2
            retry_delay += random.uniform(0, 1)
    raise Exception("Maximum retry attempts reached")



# Function to extract data from an item
def extract_data_from_item(item):
    glycosylation_type = item.get("glycosylation_type", "").title() # Correct casing for the result 
    glycosylation_type = glycosylation_type.replace("O-Glcnacylation", "O-GlcNAcylation")    
    glycosylation_subtype = ""  

    if "|" in glycosylation_type:        
        split_values = glycosylation_type.split("|")        
        glycosylation_type = split_values[0]
        glycosylation_subtype = split_values[1].split(';')[0].strip()
              

    return {
            "uniprotkb_canonical_ac": item.get("uniprot_canonical_ac", ""),
            "protein_name": item.get("protein_name", ""),
            "gene_name": "",  # Placeholder for later population
            "site": item.get("start_pos", ""),
            "amino_acid": item.get("residue", ""),
            "glycosylation_type": glycosylation_type,
            #"glycosylation_subtype": glycosylation_subtype.replace("-", ""),  # Remove any remaining '-' characters       
            "glycosylation_subtype": glycosylation_subtype,
            "peptide_seq_five_before": item.get("up_seq", "")[-5:],
            "peptide_seq_five_after": item.get("down_seq", "")[:5] if len(item.get("down_seq", "")) >= 5 else "",
            "peptide_seq_ten_before": item.get("up_seq", "")[-10:],
            "peptide_seq_ten_after": item.get("down_seq", "")[:10] if len(item.get("down_seq", "")) >= 10 else "",
            "organism": item.get("organism", ""),
            "taxonomy_id": item.get("tax_id", ""),
            "molecular_function": "",  # Placeholder for later population
            "biological_process": "",  # Placeholder for later population
            "cellular_component": "",  # Placeholder for later population
            "domain":  "",  # Placeholder for later population
            "range":  "",  # Placeholder for later population
            "status": "known_glycosite",  # Add the status field

    }

# Initial API call to get the list_id
search_data = fetch_with_backoff(SEARCH_URL, PAYLOAD)
list_id = search_data["results_summary"]["site"]["list_id"]

def all_supersearch_list_data(use_cache, logger):
    if os.path.isfile("./output/supersearch_results.tsv") and use_cache:
        logger.info("found cache file")
        return
    
    logger.info("overwriting cache")
    limit = 100
    offset = 1
    total_length, list_id = supersearch_search()

# Write initial headers to the output file
    columns = [ "uniprotkb_canonical_ac", "protein_name", "gene_name", "site", "amino_acid", "glycosylation_type",
                "glycosylation_subtype", "glytoucan_ac", "reducing_end_monosaccharide",
                "core_fucosylated", "source_tissue", "peptide_seq_five_before", "peptide_seq_five_after", 
                "peptide_seq_ten_before", "peptide_seq_ten_after","organism", "taxonomy_id","molecular_function",
                "biological_process","cellular_component", "domain", "range", "status"]
   
    with open(OUTPUT_PATH, 'w', newline='') as outfile:
        writer = DictWriter(outfile, fieldnames=columns, delimiter='\t')
        writer.writeheader()

# Fetch data and write rows
    while offset < 100:
    #while offset < total_length:
        list_payload = {"id": list_id, "offset": offset, "sort": "hit_score", "limit": limit, "order": "desc", "filters": []}
        full_list_data = fetch_with_backoff(LIST_URL, list_payload, max_retries=5)
        filtered_list = [item for item in full_list_data['results'] if "S-linked glycosylation" not in item.get("S-linked", "") and "C-linked glycosylation" not in item.get("C-linked", "")]

        with open(OUTPUT_PATH, 'a', newline='') as outfile:
            writer = DictWriter(outfile, fieldnames=columns, delimiter='\t')
            for item in filtered_list:
                extracted_data = extract_data_from_item(item)
                writer.writerow(extracted_data)

        offset += limit
        print(f"Processed up to offset: {offset}")


