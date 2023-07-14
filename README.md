### Introduction
This reposiotory contains the custom code used to perform genetic context analysis for the manuscript COMMSBIO-23-0244A. Alongside the main script used to run the analysis ("context_analysis.smk"), files containing user instructions, auxiliary scripts and example data are also provided.

### Dependencies
To run the scripts in this repository, the following software is required:
- Python >= 3.11.4
    - pandas >= 2.0.3
    - Biopython >= 1.81
    - gdown >= 4.7.1
- Snakemake >= 5.32.0
- R >= 4.2.0
    - Taxonomizr >= 0.10.2

Additionally, the following software should be located in your $PATH:
- Diamond >= 0.9.24
- Prodigal >= 2.6.3
- CD-hit >= 4.7
- mafft >= 7.3.10
- FastTree >= 2.1.9
- HMMER >= 3.3.2
- EMBOSS >= 6.5.7
- BLAST >= 2.10.1

### Tutorial
Below is a step-by-step guide on how to run genetic context analysis on a small subset of predicted APH(3') sequences. Note that some databases (e.g. NCBI Taxonomy) have been updated since the study described in COMMSBIO-23-0244A was performed, therefore some inconsistencies may appear. The pipeline is designed to run on Unix based servers.

1. Clone the repository using

```
git clone https://github.com/davidgllund/ARG_analysis_scripts.git
```

2. Download the files necessary to run the pipeline using

```
python scripts/download_auxiliary_files.py
```
    
3. Cluster the protein sequences into families (<70\% amino acid identity) using

```
snakemake -s scripts/preprocess_fargene_output.smk --cores [number of cores] all
```

4. For each cluster, generate files containing nuelotide and proteins sequence(s), as well as taxonomy and NCBI accession ids of host genome(s).

```
snakemake -s scripts/compile_cluster_data.smk --cores [number of cores] all
```

5. Retrieve genetic regions of up to 10kb up- and downstream of the genes in each cluster from their host genome(s) using 

```
source scripts/retrieve_contexts.sh
```

6. Run genetic context analysis using

```
snakemake -s scripts/context_analysis.smk --cores [number of cores] --use-conda --conda-frontend conda all
```

The results of the analysis can then be found in the file called "context_analysis_results.txt". This is a large table describing the closest known homolog (with \% amino acid identity), host taxonomy (and if any pathogenic species are present among the hosts), and any identified mobile genetic elements and/or co-localized mobile resistance genes for the genes in each cluster.
