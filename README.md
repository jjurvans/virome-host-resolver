# Wastewater Virome Host-Resolution Pipeline






<img width="1579" height="360" alt="Githubi 001" src="https://github.com/user-attachments/assets/9fd7a3ad-9b25-48ed-9961-743f1a9d53a8" />






This repository contains script for processing wastewater virome detections, inferring most probable viral hosts using NCBI Entrez metadata, extracting unique host names, and resolving them to canonical taxonomy using GNverifier and GBIF.

## Overview

The pipeline performs the following steps:

1. Combine viral detection tables and group records by accession or label.
2. Query NCBI Entrez (GenBank) to extract host metadata.
3. Infer the most probable host per virus group.
4. Export:
   - `most_probable_hosts_pooled.tsv` — pooled groups with inferred hosts  
   - `pool_members.tsv` — mapping of rows to pooled groups
5. Extract unique host names.
6. Resolve each host name to a canonical scientific name and full taxonomic classification using:
   - GNverifier API  
   - GBIF Species Match API  

## Requirements
Python ≥ 3.9  
Required libraries:
- pandas
- biopython
- requests
- tqdm

## Requirements
These scripts use:
- NCBI Entrez API (Biopython)
- GNverifier API
- GBIF Species Match API

You must set:
- Entrez.email = "your_email@example.com"
- Entrez.api_key = "YOUR_NCBI_API_KEY"

## Input file
EsViritu pipeline output with an accession column or species name column.

## Scripts
1. Host inference Generates the two main output files:
most_probable_hosts_pooled.tsv
pool_members.tsv
Includes:
Accession-based host extraction
Label-based fallback search
Majority-vote host assignment
Confidence scoring
Evidence tiering (T1 accession > T2 label > T0 unknown)
2. Unique-host extractor reads most_probable_hosts_pooled.tsv and produces:
unique_most_probable_hosts.txt
Alphabetically sorted.
3. Taxonomic resolver uses GNverifier + GBIF to generate:
resolved_hosts_full_taxonomy_final.csv
Adds scientific name, taxonomic rank, kingdom → species classification, and flags: ok, no_family, or cannot_interpret

## Citation
If you use this pipeline, please cite the repository and the APIs used:
- Biopython: Cock PJA, Antao T, Chang JT, et al. Biopython: Freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009.
- NCBI Entrez: Sayers EW, et al. Database resources of the National Center for Biotechnology Information. Nucleic Acids Research 48(D1):D9–D16 (2020). https://doi.org/10.1093/nar/gkz899
- GNverifier: Pyle RL. Global Names Architecture: A unified approach to managing the diversity of life. Bionomina 2016, 10: 1–11. https://doi.org/10.11646/bionomina.10.1.1
- GBIF Species API: GBIF Secretariat. GBIF Species API. https://api.gbif.org/v1/species

### Contact

This repository is maintained by Jaana Jurvansuu, PhD, Adjunct Professor, AIvoa — Applied AI and Modelling for Biological Sciences.  

For collaboration or implementation enquiries, please connect via [LinkedIn](https://www.linkedin.com/in/jaanajurvansuu/).

### **Note**

The original scripts were generated with assistance from ChatGPT (OpenAI) and has been subsequently reviewed, modified, and tested by the repository author.

