import requests
import pandas as pd
import re
from functools import reduce
import time
import random
import os
import logging
import json
from csv import DictReader, DictWriter

def get_positions_excluding_sites(fasta_seq, df, protein_accession):
    # Get the specific sites to exclude
    sites_to_exclude = df.loc[df.uniprotkb_canonical_ac == protein_accession, 'site'].astype(int).values
    
    # Find all motifs "NXSX" or "NXTX" where X is not P
    motif_pattern = re.compile(r'N[^P](S|T)[^P]')
    motifs = [(m.start(), m.group()) for m in motif_pattern.finditer(fasta_seq)]
    
    # Find positions of 'S' excluding specific sites and those in motifs
    s_res = [m.start() for m in re.finditer('S', fasta_seq) if (m.start() + 1) not in sites_to_exclude and not any(m.start() in range(pos, pos + 4) for pos, _ in motifs)]
    
    # Find positions of 'T' excluding specific sites and those in motifs
    t_res = [m.start() for m in re.finditer('T', fasta_seq) if (m.start() + 1) not in sites_to_exclude and not any(m.start() in range(pos, pos + 4) for pos, _ in motifs)]
    
    # Find positions of 'Y' excluding specific sites and those in motifs
    y_res = [m.start() for m in re.finditer('Y', fasta_seq) if (m.start() + 1) not in sites_to_exclude and not any(m.start() in range(pos, pos + 4) for pos, _ in motifs)]
    
    return s_res, t_res, y_res, motifs

def sites_data(logger, df):
    # Get unique uniprotkb_canonical_ac values
    unique_ac = df['uniprotkb_canonical_ac'].unique()
    unique_ac_list = unique_ac.tolist()
    print(len(unique_ac_list))

    # Load the cache from the JSON file
    with open("./output/cache.json", "r") as cache_file:
        cache = json.load(cache_file)

    count = 1
    start_time = time.time()

    for protein_accession in unique_ac_list:
        if count % 1000 == 0:
            logger.info("elapse time " + str(time.time() - start_time))
            logger.info("Hit checkpoint at index " + str(count))
            logger.info("============================================")
            start_time = time.time()
            
        target_data = cache[protein_accession]
        fasta_seq = target_data['fasta_seq']
        gene_name = target_data['gene_name']
        organism = target_data['organism']
        biological_process = target_data['biological_process']
        molecular_function = target_data['molecular_function']
        cellular_component = target_data['cellular_component']

        s_res, t_res, y_res, motifs = get_positions_excluding_sites(fasta_seq, df, protein_accession)

        protein_name = df.loc[df.uniprotkb_canonical_ac == protein_accession, 'protein_name'].values[0]
        tax_id = str(df.loc[df.uniprotkb_canonical_ac == protein_accession, 'taxonomy_id'].values[0])
        
        for s in s_res:
            s_list = [
                protein_accession,
                protein_name,
                gene_name,  # Placeholder for later population
                str(s + 1),
                "Ser",
                "",
                "",
                "",
                "",
                "",
                "",
                fasta_seq[max(0, s - 5):s],
                fasta_seq[s + 1:s + 6],
                fasta_seq[max(0, s - 10):s],
                fasta_seq[s + 1:s + 11], # Placeholder for later population
                organism,
                tax_id,
                molecular_function,  # Placeholder for later population
                biological_process,  # Placeholder for later population
                cellular_component,
                "",
                "",
                "unknown_glycosite"  # Placeholder for later population
            ]
        
            s_row = "\t".join(s_list) + "\n"
            with open("./output/supersearch_results.tsv", "a") as file:
                file.write(s_row)

        for t in t_res:
            t_list = [
                protein_accession,
                protein_name,
                gene_name,
                str(t + 1),
                "Thr",
                "",
                "",
                "",
                "",
                "",
                "",
                fasta_seq[max(0, t - 5):t],
                fasta_seq[t + 1:t + 6],
                fasta_seq[max(0, t - 10):t],
                fasta_seq[t + 1:t + 11],
                organism,
                tax_id,
                molecular_function,
                biological_process,
                cellular_component,
                "",
                "",
                "unknown_glycosite"
            ]
            t_row = "\t".join(t_list) + "\n"
            with open("./output/supersearch_results.tsv", "a") as file:
                file.write(t_row)

        for y in y_res:
            y_list = [
                protein_accession,
                protein_name,
                gene_name,
                str(y + 1),
                "Tyr",
                "",
                "",
                "",
                "",
                "",
                "",
                fasta_seq[max(0, y - 5):y],
                fasta_seq[y + 1:y + 6],
                fasta_seq[max(0, y - 10):y],
                fasta_seq[y + 1:y + 11],
                organism,
                tax_id,
                molecular_function,
                biological_process,
                cellular_component,
                "",
                "",
                "unknown_glycosite"
            ]
            y_row = "\t".join(y_list) + "\n"
            with open("./output/supersearch_results.tsv", "a") as file:
                file.write(y_row)

        for motif_pos, motif in motifs:
            if motif[2] == 'S':
                res_type = "Ser ARG"
            else:
                res_type = "Thr ARG"
            motif_list = [
                protein_accession,
                protein_name,
                gene_name,
                str(motif_pos + 1),
                res_type,
                "",
                "",
                "",
                "",
                "",
                "",
                fasta_seq[max(0, motif_pos - 5):motif_pos],
                fasta_seq[motif_pos + 1:motif_pos + 6],
                fasta_seq[max(0, motif_pos - 10):motif_pos],
                fasta_seq[motif_pos + 1:motif_pos + 11],
                organism,
                tax_id,
                molecular_function,
                biological_process,
                cellular_component,
                "",
                "",
                "unknown_glycosite"
            ]
            motif_row = "\t".join(motif_list) + "\n"
            with open("./output/supersearch_results.tsv", "a") as file:
                file.write(motif_row)

        count += 1
