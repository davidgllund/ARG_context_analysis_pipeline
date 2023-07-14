from Bio import SeqIO
import os

DIR, = glob_wildcards("clusters/{dir}")

def make_labels(list, name):
    labels = []
    for row in list:
        labels.append(name + "_" + row + "\n")
    
    return labels

def change_fasta_headers(new_headers, input_fasta, output_fasta):
    with open(input_fasta) as original_fasta, open(output_fasta, 'w') as corrected_fasta:
        records = SeqIO.parse(original_fasta, "fasta")
        i = 0
        for record in records:
            record.id = new_headers[i]
            record.description = ""
            i += 1
            SeqIO.write(record, corrected_fasta, "fasta")

rule all:
    input:
        "clusters"

# Change headerlines to appropriate format for GEnView
rule change_headers:
    input:
        "predicted-orfs-amino.fasta"
    output:
        acc = temp("acc.txt"),
        fasta = temp("predicted-orfs-amino-new-headers.fna")
    run:
        shell("grep '>' {input} | cut -d '_' -f 2,3 > {output.acc}")

        file1 = open(output[0], 'r')
        ids = file1.readlines()
        file1.close()

        headerlines_updated = make_labels(ids, "arg")

        change_fasta_headers(headerlines_updated, input[0], output[1])

# Cluster protein sequences at 70% AA identity
rule cluster_seqs:
    input:
        "predicted-orfs-amino-new-headers.fna"
    output:
        dir = directory("clusters"),
        lines = temp("cluster_names.txt")
    shell:
        """
        mkdir {output.dir}
        usearch -cluster_fast {input} -id 0.7 -clusters {output.dir}/cluster_
        ls {output.dir} > {output.lines}
        while read line; do mv clusters/$line clusters/$line.fna; mkdir clusters/$line; mv clusters/$line.fna clusters/$line/proteins.fna; done<{output.lines}
        """
