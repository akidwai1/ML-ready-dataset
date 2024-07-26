import requests
import pandas as pd
import re
from functools import reduce
import time
import random 
import os
from io import StringIO
import logging
import json
from csv import DictReader, DictWriter


# Define URLs and payloads
SEARCH_URL = "https://beta-api.glygen.org/supersearch/search"
LIST_URL = "https://beta-api.glygen.org/supersearch/list/"
GLYCAN_URL = "https://beta-api.glygen.org/glycan/detail/"
OUTPUT_PATH = "./output/supersearch_results.tsv" 
REDENDMONO_URL = "https://github.com/glygen-glycan-data/PyGly/raw/master/smw/glycandata/export/redendmono.tsv"


def fetch_glycan_with_backoff(ac, max_retries=5):
    retry_delay = 1
    for attempt in range(max_retries):
        try:
            response = requests.post(url=f"{GLYCAN_URL}{ac}", json={"glytoucan_ac": ac}) # to get GO:
            response.raise_for_status()
            return response.json()
        except requests.RequestException:
            time.sleep(retry_delay)
            retry_delay *= 2
            retry_delay += random.uniform(0, 1)
            raise Exception("Maximum retry attempts reached")

def get_glycan_data(logger):
    # Load the existing data
    df = pd.read_csv(REDENDMONO_URL, sep='\t')
    df1 = pd.read_csv(OUTPUT_PATH, delimiter='\t') if os.path.exists(OUTPUT_PATH) else pd.DataFrame()
    #df1 = pd.read_csv(OUTPUT_PATH, delimiter='\t')
    merged_df = pd.concat([df1, df], axis=1)
    #merged_df = df1.merge(df, left_on='glytoucan_ac', right_on='accession', how='left')
    df1['reducing_end_monosaccharide'] = merged_df['monosaccharide']
    # Get unique uniprotkb_canonical_ac values
    unique_ac = df1['glytoucan_ac'].unique()
    unique_ac_list = unique_ac.tolist()
    print(unique_ac_list)
    unique_ac_list = [x for x in unique_ac_list if str(x).lower() != "nan"]

    print(unique_ac)
    print(len(unique_ac_list))

    columns = [ "uniprotkb_canonical_ac", "protein_name","gene_name",	"site",	"amino_acid","glycosylation_type", 
                "glycosylation_subtype", "glytoucan_ac", "reducing_end_monosaccharide",
                "core_fucosylated","source_tissue", "peptide_seq_five_before", "peptide_seq_five_after", 
                "peptide_seq_ten_before", "peptide_seq_ten_after","organism", "taxonomy_id","molecular_function",
                "biological_process","cellular_component", "domain", "range", "status"]
    

 
# Core fucosylation and source tissue info
    subtypes = []
    core_fucosylated = None
    for ac in unique_ac_list:
        print(ac)
        glycan_data = fetch_glycan_with_backoff(ac, max_retries=5) # Get glycan details
        
        for classification in glycan_data.get("classification", []):
            subtype_name = classification.get("subtype", {}).get("name")
            if subtype_name:
                subtypes.append(subtype_name)
            if subtype_name == "Core-fucosylated":
                core_fucosylated = "Y"
            else:
                core_fucosylated = "N"

        
        if glycan_data.get("expression"):
            source_tissue_info = ""
            source_tissue = []
            for expression_info in glycan_data["expression"]:
                if expression_info.get("category") == "tissue":
                    tissue = expression_info.get("tissue", {})
                    name = tissue.get("name", "")
                    namespace = tissue.get("namespace", "")
                    id = tissue.get("id", "")                     
                    source_tissue.append(f"{name} ({namespace} : {id})")
            source_tissue_info = ";".join(source_tissue)
                        # Add the core_fucosylated and source_tissue fields to the DataFrame df1
        
            df1.loc[df1.glytoucan_ac == ac, ['core_fucosylated']] = core_fucosylated
            df1.loc[(df1.glytoucan_ac == ac) & (df1.uniprotkb_canonical_ac), ['source_tissue']] = source_tissue_info  #both "glytoucan_ac" and "uniprot_canonical_ac" columns match the source_tissue_info values
       
        
        # Save the final DataFrame to a file
    df1.to_csv(OUTPUT_PATH, sep='\t', index=False)
