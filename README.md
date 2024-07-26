**DOCUMENTATION**

| **Field** | **Description** | **Example** | **Value** |
| --- | --- | --- | --- |
| uniprotkb_ac | Canonical Accession of the protein from UniProtKB | P00533-1 | Single |
| protein_name | Name of the protein from UniProtKB | Epidermal growth factor receptor | Single or no value |
| gene_name | HGNC Gene name for the gene coding protein | EGFR | Single or no value |
| site | Position on the protein for the O-lined site that is known to be modified or is unknown to be modified. | 430 | Single (remove entries with no values) |
| amino_acid | Amino acid (aa) that is modified in three-letter code | Thr | Single (remove entries with no values) |
| glycosylation_type | Glycosylation modification present on the protein | O-Linked | Single |
| glycosylation_subtype | Specifies the subtype or specific form of glycosylation on the protein | O-GlcNAcylation | Single or other |
| glytoucan_ac | Unique identifier for the specific glycan structure attached to the protein | G49108TO | Single or no value |
| reducing_end_monosacchride | Indicates the monosaccharide present at the reducing end of the glycan structure | GlcNAc | Single |
| core_fucosylated | Indicates whether the protein contains a core fucose modification | Y or N | Single |
| source_tissue | This denotes the tissues from which the protein sample was obtained | gizzard (UBERON : 0005052) |     |
| o-glcnacsylation | O-Glcnacsylation observed on the residue | Y or N | Single |
| peptide_seq_five_before | Peptide sequence 5 aa before the modified residue (Upstream towards N terminal) | IIRGR | Single or no value |
| peptide_seqfive_after | Peptide sequence 5 aa after the modified residue (Downstream towards C terminal) | KQHGQ | Single or no value |
| peptide_seq_ten_before | Peptide sequence 10 aa before the modified residue (towards N terminal) | FENLEIIRGR | Single or no value |
| peptide_seqten_after | Peptide sequence 10 aa after the modified residue (towards C terminal) | KQHGQFSLAV | Single or no value |
| organism | Scientific name of the organism to which the protein belongs | Homo sapiens | Single |
| taxonomy_id | NCBI Taxonomy ID for the scientific name of the organism | 9606 | Single |
| molecular_function | Gene Ontology Molecular Function name and ID | ATPase binding (GO:0051117)\|ATP binding (GO:0005524) | Single or no value or multiple value |
| biological_process | Gene Ontology Biological Function name and ID | activation of phospholipase C activity (GO:0007202)\|astrocyte activation (GO:0048143) | Single or no value or multiple value |
| cellular_component | Gene Ontology Molecular Function name and ID | apical plasma membrane (GO:0016324)\|basal plasma membrane (GO:0009925) | Single or no value or multiple value |
| domain | Description of the domain | Protein kinase | Single or no value |
| range | Range of Amino Acids within the Domain | 150-408 | Single or no value |
| status | Whether the site is known or unknown to be glycosylated | Known_glycoste<br><br>potential_glycosite | Single or no value |

**Installation & Running Guide**

**Overview**

I created the code for Project 3 to generate an ML-ready dataset that helps predict N-linked and O-linked glycosylation sites on proteins. The motivation behind this project is to facilitate the identification of glycosylation sites, which play critical roles in various biological processes. By curating a comprehensive dataset with precise information such as protein accession numbers, glycosylation types, amino acid modifications, and relevant Gene Ontology terms, the code provides a structured format suitable for machine learning applications.

The command-line tool facilitates easy execution and generates an output file in a spreadsheet format, ensuring the dataset is user-friendly and accessible for further analysis in glycosylation prediction studies.

**Prerequisites**

**\- Python:** Ensure Python 3.x is installed on your system. Download from \[python.org\](<https://www.python.org/downloads/>).

**\- Python Packages:** Install the necessary packages using pip:

pip install requests pandas

**Installation Steps**

**1\. Clone the Repository:**

\- Clone the repository containing the Python scripts you need.

**2\. Install Dependencies:**

\- Open a terminal and run:

pip install requests pandas

\- Use \`pip3\` if you have Python 3 or above installed:

pip3 install requests pandas

**Running Instructions**

\- To run the program, execute the following command:

python main.py

\- For Python 3, use:

python3 main.py

**\- Options: You can add command-line options as needed:**

- No Cache: To restart with no cache, add \`-n\` or \`--nocache\`:

python3 main.py --nocache

- Output Format: By default, the output is in TSV format. To save as CSV, use:

python3 main.py -o csv

- API Version: The program uses the production API by default. To switch to the beta API, add:

python3 main.py -a beta

**\- Combining Options:** Multiple options can be combined. For example:

python3 main.py --nocache -o csv -a beta

**Understanding the Scripts**

- The scripts are designed to fetch and process protein and glycan information, handling API requests and outputting processed data.

**Output Files**

- Processed data is saved in TSV format by default. You can specify CSV format if needed.

**Customization**

- Adjust the scripts based on your requirements, such as modifying API endpoints or data processing logic.

**Troubleshooting**

- Ensure a stable internet connection and check for errors in script outputs or logs if issues arise.

**References**

- Ensure proper API authentication and access permissions for external APIs used in the scripts.

Overall, By following these steps, you can effectively install and run the provided Python scripts for processing protein and glycan information.
