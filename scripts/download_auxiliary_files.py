#!/usr/bin/env/python
import subprocess
import os

if not os.path.isdir("databases"):
    subprocess.run("mkdir databases", shell=True)

# Download resistance gene sequencs from ResFinder
if not os.path.isfile("databases/resfinder_translated/resfinder.fasta"):
    subprocess.run("git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git", shell=True)
    subprocess.run("mv resfinder_db/ databases", shell=True)
    subprocess.run("cat databases/resfinder_db/all.fsa | tr -d '(' | tr -d ')' > databases/resfinder.fsa", shell=True)
    subprocess.run("mkdir databases/resfinder_translated", shell=True)
    subprocess.run("transeq databases/resfinder.fsa databases/resfinder_translated/resfinder.fasta -frame=1", shell=True)
    subprocess.run("makeblastdb -in databases/resfinder_translated/resfinder.fasta -dbtype prot", shell=True)

# Download sequences from ISfinder
if not os.path.isfile("databases/ISFinder-sequences/IS.fna"):
    subprocess.run("git clone https://github.com/thanhleviet/ISfinder-sequences.git", shell=True)
    subprocess.run("mkdir databases/ISfinder-sequences", shell=True)
    subprocess.run("mv ISfinder-sequences databases", shell=True)
    subprocess.run("makeblastdb -in databases/ISfinder-sequences/IS.fna -dbtype nucl", shell=True)
    
if not os.path.isdir("index_files"):
    subprocess.run("mkdir index_files", shell=True)

# Download list of pathogenic species from PATRIC
if not os.path.isfile("index_files/pathogen_index.txt"):
    subprocess.run("wget ftp://ftp.bvbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt", shell=True)
    subprocess.run("cat PATRIC_genomes_AMR.txt | tail -n +2 | grep '.' | cut -f 2 | grep -v 'sp.' | cut -d ' ' -f 1,2 | sort | uniq > index_files/pathogen_list.txt", shell=True)
    subprocess.run("mv PATRIC_genomes_AMR.txt index_files", shell=True)

# Retrieve taxonomy of bacterial hosts
if not os.path.isfile("index_files/taxonomy.txt"):
    subprocess.run("Rscript scripts/get_taxonomy.R", shell=True)

# Download hidden Markov models from Conjscan
if not os.path.isdir("models"):
    subprocess.run("mkdir models", shell=True)
    subprocess.run("git clone https://github.com/gem-pasteur/Macsyfinder_models.git", shell=True)
    subprocess.run("mv Macsyfinder_models/models/Conjugation/ models/conjscan_models", shell=True)
    subprocess.run("rm -rf Macsyfinder_models", shell=True)

# Download GEnView v0.1
if not os.path.isdir("genview_dir"):
    subprocess.run("mkdir genview_dir", shell=True)
    subprocess.run("mkdir genview_dir/db", shell=True)
    subprocess.run("git clone https://github.com/EbmeyerSt/GEnView.git", shell=True)
    os.chdir("GEnView")
    subprocess.run("git checkout d428943", shell=True)
    os.chdir("..")
    subprocess.run("mv GEnView/ genview_dir", shell=True)

