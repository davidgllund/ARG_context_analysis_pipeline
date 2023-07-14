from Bio import SeqIO
import os

DIR, = glob_wildcards("clusters/{dir}/proteins.fna")

def extract_fasta_subset(headers_subset, headers_complete, input_fasta, output_fasta):
    with open(input_fasta) as original_fasta, open(output_fasta, "w") as corrected_fasta:
        records = SeqIO.parse(original_fasta, "fasta")
        i = 0
        for record in records:
            if headers_complete[i] not in headers_subset:
                i += 1
                continue
            else:
                SeqIO.write(record, corrected_fasta, "fasta")
                i += 1
                continue

rule all:
    input:
        expand("clusters/{dir}/nucleotides.fna", dir=DIR),
        expand("clusters/{dir}/assembly_accessions.txt", dir=DIR),
        expand("clusters/{dir}/species.txt", dir=DIR),
        expand("clusters/{dir}/unique_proteins.fna", dir=DIR),
        expand("clusters/{dir}/pathogens.txt", dir=DIR)

rule extract_nucl_headers:
    input:
        "predicted-orfs.fasta"
    output:
        temp("header_nucl.txt")
    shell:
        """
        grep '>' {input} > {output}
        """

# Extract the nucleotide sequences corresponding to the proteins in each cluster
rule get_nucl_fasta:
    input:
        prot = "clusters/{dir}/proteins.fna",
        header1 = "header_nucl.txt",
        nucl = "predicted-orfs.fasta"
    output:
        acc = temp("clusters/{dir}/acc.txt"),
        header2 = temp("clusters/{dir}/header2.txt"),
        fasta = "clusters/{dir}/nucleotides.fna"
    run:
        shell("grep '>' {input.prot} | cut -d '_' -f 2,3 > {output.acc}")
        shell("while read line; do grep $line {input.header1} >> {output.header2}; done<{output.acc}")

        file1 = open(input[1])
        headers_complete = file1.readlines()
        file1.close()

        file2 = open(output[1])
        headers_subset = file2.readlines()
        file2.close()

        extract_fasta_subset(headers_subset, headers_complete, input[2], output[2])

rule extract_accno:
    input:
        "clusters/{dir}/nucleotides.fna"
    output:
        "clusters/{dir}/assembly_accessions.txt"
    shell:
        """
        set +o pipefail
        grep '>' {input} | cut -d '_' -f 1,2 | tr -d '>' > {output}
        """

# Get taxonomic information of hosts for the proteins in each cluster
rule get_taxonomy:
    input:
        fasta = "clusters/{dir}/proteins.fna",
        tax = "index_files/taxonomy.txt"
    output:
        ids = temp("clusters/{dir}/acc.txt"),
        tax = temp("clusters/{dir}/taxonomy.txt"),
        sp = "clusters/{dir}/species.txt",
        ph = "clusters/{dir}/phylum.txt",
	    cl = "clusters/{dir}/class.txt"
    shell:
        """
        grep '>' {input.fasta} | cut -d '_' -f 2 > {output.ids}

        set +o pipefail
        while read line; do 
            grep $line {input.tax} | head -1 >> {output.tax} 
        done<{output.ids}

        cat {output.tax} | cut -f 8 > {output.sp}
        cat {output.tax} | cut -f 4 > {output.cl}
        cat {output.tax} | cut -f 3 > {output.ph}
        """

rule cluster_proteins:
    input:
        "clusters/{dir}/proteins.fna"
    output:
        "clusters/{dir}/unique_proteins.fna"
    shell:
        """
        usearch -cluster_fast {input} -id 1 -centroids {output}
        """

rule find_pathogens:
    input:
        db = "index_files/pathogen_list.txt",
        sp = "clusters/{dir}/species.txt"
    output:
        "clusters/{dir}/pathogens.txt"
    shell:
        """
        grep -f {input.db} {input.sp} >> {output}
        """