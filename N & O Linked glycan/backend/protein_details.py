import requests
import pandas as pd
import re
from functools import reduce
import time
import random
import os
import logging
from io import StringIO
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from csv import DictReader, DictWriter


DETAIL_URL = "https://api.glygen.org/protein/detail/"

# Function to get protein details and additional data
#def get_protein_details(protein_accession, max_retries=5):
def get_protein_details(protein_accession, max_retries=5):
    retry_delay = 1
    for attempt in range(max_retries):
        try:
            protein_response = requests.post(url=f"{DETAIL_URL}{protein_accession}", json={"uniprotkb_canonical_ac": protein_accession})
            protein_response.raise_for_status()
            return protein_response.json()
        except requests.RequestException:
            time.sleep(retry_delay)
            retry_delay *= 2
            retry_delay += random.uniform(0, 1)
    raise Exception("Maximum retry attempts reached")

# Function to get protein details and additional data
def get_domain_details(protein_accession, max_retries=5):
#def get_domain_details(protein_accession, max_retries=5):
    retry_delay = 1  # Initial delay in seconds
    for attempt in range(max_retries):
        try:
            protein_response = requests.get(url=f"https://www.ebi.ac.uk/proteins/api/proteins/{protein_accession}", headers={'Accept': 'application/json'})
            protein_response.raise_for_status()
            return protein_response.json()
        except requests.RequestException:
            print(protein_accession)
            time.sleep(retry_delay)
            retry_delay *= 2  # Double the delay for the next attempt
            retry_delay += random.uniform(0, 1)  # Add jitter

    raise Exception("Maximum retry attempts reached")

def find_gene_names_and_organism(protein_data):
    gene_name = ""
    glytoucan_ac = ""
    organism = ""

    if protein_data.get('gene_names'):
        for g in protein_data['gene_names']:
            if (g['resource'] == 'UniProtKB') and (g['type'] == 'recommended'):
                gene_name = g['name']
                break
            else:
                if (g['resource'] == 'UniProtKB'):
                    gene_name = g['name']
                    break

    if protein_data.get('glycosylation'):
        for gly in protein_data["glycosylation"]:
            glytoucan_ac = gly.get('glytoucan_ac')
    
    if protein_data.get('species'):
        for o in protein_data['species']:
            organism = o['name']
    return gene_name, glytoucan_ac, organism

def find_domain_and_range(domain_data):
    domain_dict = {}
    
    if domain_data.get('features'):
        for d in domain_data['features']:
            if d['type'] == "DOMAIN":
                begin = d.get('begin')
                end = d.get('end')
                end = end.lstrip(">")
                desc = d.get('description')
                range = f'{begin}-{end}' 
                domain_dict[range] = desc
    return domain_dict

def find_range_and_value(data, value):
    try:
        for key, val in data.items():
            start, end = map(int, key.split('-'))
            if start <= value <= end:
                return key, val
        return "", ""
    except Exception as e:
        return '', ''


# Function to find GO annotations

def find_go_annotations(protein_accession, organism):
    linked_csv_files = {
        'Homo sapiens': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/human_protein_go_annotation.csv',
        'Mus musculus': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/mouse_protein_go_annotation.csv',
        'Gallus gallus': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/chicken_protein_go_annotation.csv',
        'Dictyostelium discoideum': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/dicty_protein_go_annotation.csv',
        'Drosophila melanogaster': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/fruitfly_protein_go_annotation.csv',
        'Hepacivirus C genotype 1a': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/hcv1a_protein_go_annotation.csv',
        'Hepacivirus C genotype 1b': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/hcv1b_protein_go_annotation.csv',
        'Sus scrofa': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/pig_protein_go_annotation.csv',
        'Rattus norvegicus': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/rat_protein_go_annotation.csv',
        'Severe acute respiratory syndrome coronavirus 1': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/sarscov1_protein_go_annotation.csv',
        'Severe acute respiratory syndrome coronavirus 2': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/sarscov2_protein_go_annotation.csv',
        'Saccharomyces cerevisiae S288C': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/yeast_protein_go_annotation.csv'
    }

    logger = logging.getLogger(__name__)
    
    link = linked_csv_files.get(organism)
    if not link:
        logger.error(f"No CSV file link found for organism: {organism}")
        return '', '', ''

    try:
        response = requests.get(link, verify=True)  # Ensure SSL verification
        response.raise_for_status()  # Raise an exception for bad responses

        # Use StringIO to convert the response content to a file-like object
        csv_data = response.content.decode('utf-8')
        linked_df = pd.read_csv(StringIO(csv_data))  # Read CSV from string
    except Exception as e:
        logger.error(f"Error reading CSV file for organism: {organism} from link: {link}. Error: {e}")
        return '', '', ''
    biological_list, molecular_list, cellular_list = [], [], []
    categories = ["molecular_function", "biological_process", "cellular_component"]

    try:
        matched_df = linked_df[linked_df['uniprotkb_canonical_ac'] == protein_accession]
    except Exception as e:
        logger.error(f"Error reading CSV file for organism: {organism} from link: {link}. Error: {e}")
        return '', '', ''

    matched_df.loc[:, 'go_term_id'] = matched_df['go_term_id'].astype(str).str.replace('_', ':')
    for category in categories:  # extracting info based on categories into the corresponding list
        go_terms = matched_df.loc[matched_df['go_term_category'] == category, ['go_term_label', 'go_term_id']]
        if category == "molecular_function":
            for index, row in go_terms.iterrows():
                molecular_list.append(f"{row['go_term_label']} ({row['go_term_id']})")
        elif category == "biological_process":
            for index, row in go_terms.iterrows():
                biological_list.append(f"{row['go_term_label']} ({row['go_term_id']})")
        elif category == "cellular_component":
            for index, row in go_terms.iterrows():
                cellular_list.append(f"{row['go_term_label']} ({row['go_term_id']})")

    return '; '.join(map(str, molecular_list)), '; '.join(map(str, biological_list)), '; '.join(map(str, cellular_list))

def process_protein(row, cache, logger):
    start_time = time.time()
    uniprot = row["uniprotkb_canonical_ac"]
    protein_ac = uniprot.replace("-1", "")
    protein_ac = re.sub(r'-\d+', '', protein_ac)
    protein_data = get_protein_details(uniprot, max_retries=5) # Your function to fetch protein details
    gene_name, glytoucan_ac, organism = find_gene_names_and_organism(protein_data)  # Your function to find gene names and kinases
    molecular_function, biological_process, cellular_component = find_go_annotations(uniprot, organism)  # Your function to find GO annotations
    if protein_ac:
        try:
            domain_data = get_domain_details(protein_ac, max_retries=5)
            domain_dict = find_domain_and_range(domain_data)
            range, domain = find_range_and_value(domain_dict, int(row['site']))
        except Exception as e:
            domain = ''
            range = ''
            domain_dict = {}
    gene_name, glytoucan_ac, organism = find_gene_names_and_organism(protein_data)
    molecular_function, biological_process, cellular_component = find_go_annotations(uniprot, organism)
    fasta_seq = protein_data["sequence"]['sequence']

    target_data = {
        "gene_name": gene_name,
        "organism": organism,
        "glytoucan_ac":glytoucan_ac,
        "molecular_function": molecular_function,
        "biological_process": biological_process,
        "cellular_component": cellular_component,
        "fasta_seq": fasta_seq,
        "domain_dict" : domain_dict
    }

    row["gene_name"] = gene_name
    row["glytoucan_ac"] = target_data["glytoucan_ac"]
    row['organism'] = organism
    row["molecular_function"] = molecular_function
    row["biological_process"] = biological_process
    row["cellular_component"] = cellular_component
    row["range"] = range
    row["domain"] = domain

    logger.info(f"Processed {uniprot} in {time.time() - start_time} seconds")

    with open("./output/cache.json", "w") as cache_file:
        json.dump(cache, cache_file)

    return row, target_data

def get_all_protein_data(use_cache, logger):
    
    if use_cache and os.path.isfile("./output/cache.json"):
         with open("./output/cache.json", "r") as cache_file:
            cache = json.load(cache_file)
    else:
        cache = {}

    cache_checkpoint = 100
# Write initial headers to the output file
    columns = [ "uniprotkb_canonical_ac", "protein_name","gene_name", "site", "amino_acid","glycosylation_type", 
                 "glycosylation_subtype", "glytoucan_ac", "reducing_end_monosaccharide",
                "core_fucosylated", "source_tissue", "peptide_seq_five_before", "peptide_seq_five_after", 
                "peptide_seq_ten_before","peptide_seq_ten_after","organism", "taxonomy_id","molecular_function",
                "biological_process","cellular_component","domain", "range", "status"]

    outfile_path = "./output/supersearch_results.tsv"

    with open(outfile_path, "r") as infile, open(outfile_path + ".tmp", "w", newline='') as outfile:
        reader = DictReader(infile, delimiter='\t')
        writer = DictWriter(outfile, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        api_hit = 0
#######
        # Use ThreadPoolExecutor for concurrent API requests
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = []

            for i, row in enumerate(reader):

                uniprot = row["uniprotkb_canonical_ac"]
                if uniprot in cache:
                    target_data = cache[uniprot]
                    row["gene_name"] = target_data["gene_name"]
                    row["glytoucan_ac"] = target_data["glytoucan_ac"]
                    row['organism'] = target_data['organism']
                    row["molecular_function"] = target_data["molecular_function"]
                    row["biological_process"] = target_data["biological_process"]
                    row["cellular_component"] = target_data["cellular_component"]
                    range, domain = find_range_and_value(target_data['domain_dict'], int(row['site']))
                    row["range"] = range
                    row["domain"] = domain
                    writer.writerow(row)
                else:
                    # Submit the task to the executor
                    futures.append(executor.submit(process_protein, row, cache, logger))
                
                if i % cache_checkpoint == 0:
                    results = [future.result() for future in as_completed(futures)]
                    for row, target_data in results:
                        writer.writerow(row)
                        cache[row["uniprotkb_canonical_ac"]] = target_data
                        api_hit += 1

                    # Clear the futures list after processing
                    futures.clear()

                    # Save the cache to file
                    with open("./output/cache.json", "w") as cache_file:
                        json.dump(cache, cache_file)
                    
                    logger.info("Hit checkpoint at index " + str(i))
                    logger.info("API calls " + str(api_hit))
                    logger.info("cache hits " + str(cache_checkpoint - api_hit))
                    logger.info("============================================")
                
                    api_hit = 0

            if futures:
                results = [future.result() for future in as_completed(futures)]
                for row, target_data in results:
                    writer.writerow(row)
                    cache[row["uniprotkb_canonical_ac"]] = target_data
                    api_hit += 1

            # Clear futures list after final processing
            futures.clear()

    # Replace the original file with the updated file
    os.replace(outfile_path + ".tmp", outfile_path) # Replace the original file with the updated fil

    # Save cache
    with open("./output/cache.json", "w") as cache_file:
        json.dump(cache, cache_file)


