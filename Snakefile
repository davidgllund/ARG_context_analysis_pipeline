#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORT PACKAGES AND DEFINE FUNCTIONS
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import re
import sys
import csv
from Bio import SeqIO

# Define list of wildcards to make rules apply to all subdirectories
DIR, = glob_wildcards("{dir}/contexts.fna")

# Function that, given a list, counts the number of unique entries and returns the top 3 most frequently occurring entries from the list.
def get_top_items(list):
    unique_names = []

    for name in list:
        if name not in unique_names:
            unique_names.append(name)

    number_of_occurrences = [0]*len(unique_names)

    for i in range(len(list)):
        for j in range(len(unique_names)):
            if list[i] == unique_names[j]:
                number_of_occurrences[j] += 1

    index = np.argsort(-np.array(number_of_occurrences))
    
    if len(unique_names) > 3:
        top_index = index[0:3]
    else:
        top_index = index

    names = ""

    for i in range(len(top_index)):
        if i == (len(top_index)-1) and len(unique_names) > 3:
            names += unique_names[top_index[i]].rstrip() + "..."
        elif i == (len(top_index)-1) and len(unique_names) <= 3:
            names += unique_names[top_index[i]].rstrip()
        else:
            names += unique_names[top_index[i]].rstrip() + "\n"

    if names == "":
        names += "NA"

    return names
    
def summarize_blast_results(results,threshold):
    elements=[]
    summary=[]

    # Define lists to hold input data
    query=[]
    subject=[]
    identity=[]
    aln_length=[]
    s_start=[]
    s_end=[]

    with open(results) as tsv:
        for line in csv.reader(tsv, delimiter="\t"):
            # Remove extraneous lines from results
            hit = re.search("#", str(line))
            if hit is None:
                # Only keep hits that are above threshold in terms of AA identity
                if float(line[2]) >= float(threshold):
                    query.append(line[0])
                    subject.append(line[1])
                    identity.append(float(line[2]))
                    aln_length.append(float(line[3]))
                    s_start.append(int(line[8]))
                    s_end.append(int(line[9]))

    # Divide hits according to context in which they were found
    contexts=[]

    # Remove duplicate context names
    for name in subject:
        if name not in contexts:
            contexts.append(name)

    # For each genetic context, find the best hit among overlapping hits, report these together with their corresponding AA identity
    for name in contexts:
        tmp_query=[]
        tmp_id=[]
        tmp_start=[]
        tmp_end=[]
        tmp_range=[]

        for i in range(len(subject)):
            hit = re.search(subject[i], name)
            if hit is not None:
                tmp_query.append(query[i])
                tmp_id.append(identity[i])
                tmp_start.append(min(s_start[i],s_end[i]))
                tmp_end.append(max(s_start[i],s_end[i]))
                tmp_range.append((min(s_start[i],s_end[i]),max(s_start[i],s_end[i])))

        df = np.transpose(pd.DataFrame([tmp_query,tmp_id,tmp_start,tmp_end]))

        # Sort the intervals by lower bound
        sorted_intervals = sorted(tmp_range, key=lambda tup: tup[0])
        overlapping = []

        for higher in sorted_intervals:
            if not overlapping:
                # Add first intervall to the list of overlapping intervals
                overlapping.append(higher)
            else:
                # Get current interval with the highest values
                lower = overlapping[-1]
                # Check if intervals overlap (lower bound of new interval <= higher bound of current interval)
                if higher[0] <= lower[1]:
                    # Define new upper bound of overlapping interval
                    upper_limit = max(lower[1], higher[1])
                    # Replace current highest interval with combination of overlapping intervals
                    overlapping[-1] = (lower[0], upper_limit)
                else:
                    # If no overlap found, add new interval to list of overlapping intervals
                    overlapping.append(higher)

        for interval in overlapping:
            overlap_query=[]
            overlap_id=[]
            for i in range(len(df.iloc[:,0])):
                if df.iloc[i,2] in range(interval[0],interval[1]) or df.iloc[i,3] in range(interval[0],interval[1]):
                    overlap_query.append(df.iloc[i,0])
                    overlap_id.append(df.iloc[i,1])

            hits=np.transpose(pd.DataFrame([overlap_query, overlap_id]))
            best_hits=hits.sort_values(by=[1], ascending=False)

            elements.append((best_hits.iloc[0,0], best_hits.iloc[0,1]))

    # Compile final results
    unique = []

    for item in elements:
        if item[0] not in unique:
            unique.append(item[0])

    for name in unique:
        j = 0
        percentages = []

        for item in elements:
            hit = re.search(name,item[0])

            if hit is not None:
                j += 1
                percentages.append(item[1])

        if min(percentages) == max(percentages):
            summary.append(name + " [" + str(max(percentages)) + "%] (" + str(j) + ")" + "\n")
        else:
            summary.append(name + " [" + str(max(percentages)) + "-" + str(min(percentages)) + "%] (" + str(j) + ")" + "\n")

    return summary


def summarize_blast_results_distance_considered(results,position,threshold,distance):
    elements=[]
    summary=[]
    
    df_pos = pd.read_csv(position, delimiter = '\t', names=("Context","Start","End"))

    # Define lists to hold input data
    query=[]
    subject=[]
    identity=[]
    aln_length=[]
    s_start=[]
    s_end=[]

    with open(results) as tsv:
        for line in csv.reader(tsv, delimiter="\t"):
            # Remove extraneous lines from results
            hit = re.search("#", str(line))
            if hit is None:
                # Only keep hits that are above threshold in terms of AA identity
                if float(line[2]) >= float(threshold):
                    query.append(line[0])
                    subject.append(line[1])
                    identity.append(float(line[2]))
                    aln_length.append(float(line[3]))
                    s_start.append(int(line[8]))
                    s_end.append(int(line[9]))
                
    # Divide hits according to context in which they were found
    contexts=[]

    # Remove duplicate context names
    for name in subject:
        if name not in contexts:
            contexts.append(name)
        
    # For each genetic context, find the best hit among overlapping hits, report these together with their corresponding AA identity
    for name in contexts:
        bool=df_pos.Context.str.contains(name)
        arg_start = df_pos[bool].min(axis=1)
        arg_end = df_pos[bool].max(axis=1)
    
        if sum(bool) == 0:
            continue
        else:
            tmp_query=[]
            tmp_id=[]
            tmp_start=[]
            tmp_end=[]
            tmp_range=[]

            for i in range(len(subject)):
                if arg_end.iloc[0] > max(s_start[i],s_end[i]) and abs(arg_start.iloc[0]-max(s_start[i],s_end[i]))<int(distance) or arg_start.iloc[0] < min(s_start[i],s_end[i]) and abs(min(s_start[i],s_end[i])-arg_end.iloc[0])<int(distance):
                    hit = re.search(subject[i], name)
                    if hit is not None:
                        tmp_query.append(query[i])
                        tmp_id.append(identity[i])
                        tmp_start.append(min(s_start[i],s_end[i]))
                        tmp_end.append(max(s_start[i],s_end[i]))
                        tmp_range.append((min(s_start[i],s_end[i]),max(s_start[i],s_end[i])))

        df = np.transpose(pd.DataFrame([tmp_query,tmp_id,tmp_start,tmp_end]))

        # Sort the intervals by lower bound
        sorted_intervals = sorted(tmp_range, key=lambda tup: tup[0])
        overlapping = []

        for higher in sorted_intervals:
            if not overlapping:
                # Add first intervall to the list of overlapping intervals
                overlapping.append(higher)
            else:
                # Get current interval with the highest values
                lower = overlapping[-1]
                # Check if intervals overlap (lower bound of new interval <= higher bound of current interval)
                if higher[0] <= lower[1]:
                    # Define new upper bound of overlapping interval
                    upper_limit = max(lower[1], higher[1])
                    # Replace current highest interval with combination of overlapping intervals
                    overlapping[-1] = (lower[0], upper_limit)
                else:
                    # If no overlap found, add new interval to list of overlapping intervals
                    overlapping.append(higher)

        for interval in overlapping:
            overlap_query=[]
            overlap_id=[]
            for i in range(len(df.iloc[:,0])):
                if df.iloc[i,2] in range(interval[0],interval[1]) or df.iloc[i,3] in range(interval[0],interval[1]):
                    overlap_query.append(df.iloc[i,0])
                    overlap_id.append(df.iloc[i,1])

            hits=np.transpose(pd.DataFrame([overlap_query, overlap_id]))
            best_hits=hits.sort_values(by=[1], ascending=False)

            elements.append((best_hits.iloc[0,0], best_hits.iloc[0,1]))
        
    # Compile final results
    unique = []

    for item in elements:
        if item[0] not in unique:
            unique.append(item[0])

    for name in unique:
        j = 0
        percentages = []

        for item in elements:
            hit = re.search(name,item[0])

            if hit is not None:
                j += 1
                percentages.append(item[1])

        if min(percentages) == max(percentages):
            summary.append(name + " [" + str(max(percentages)) + "%] (" + str(j) + ")" + "\n")
        else:
            summary.append(name + " [" + str(max(percentages)) + "-" + str(min(percentages)) + "%] (" + str(j) + ")" + "\n")

        return summary

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
# RULES FOR COMPILING METADATA ABOUT GENE FAMILIES
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
# "ALL"-rule, run to get final results
#----------------------------------------------------------------------------------------------------------------------------------------------------------
rule all:
    input:
        expand("{dir}/metadata.txt", dir=DIR),
        "metadata_complete.txt"

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate auxillary data
#----------------------------------------------------------------------------------------------------------------------------------------------------------
rule translate_contexts:
    input:
        "{dir}/contexts.fna"
    output:
        prot = temp("{dir}/contexts_protein.fna")
    shell:
        """
        set +o pipefail
        transeq {input} {output.prot} -frame=6
        """

rule create_blast_databases:
    input:
        nucl = "{dir}/contexts.fna",
        prot = "{dir}/contexts_protein.fna"
    output:
        nucl = temp("{dir}/db_nucl/contexts.fna"),
	prot = temp("{dir}/db_prot/contexts.fna")
    shell:
        """
        set +o pipefail

        if [[ ! -d {wildcards.dir}/db_nucl ]]; then
            mkdir {wildcards.dir}/db_nucl
        fi
        
        cp {input.nucl} {output.nucl}
        makeblastdb -in {output.nucl} -dbtype nucl

        if [[ ! -d {wildcards.dir}/db_prot ]]; then
            mkdir {wildcards.dir}/db_prot
        fi
        
        cp {input.prot} {output.prot}
        makeblastdb -in {output.prot} -dbtype prot
        """

rule find_arg_positions:
    input:
        genes = "{dir}/nucleotides.fna",
        contexts = "{dir}/db_nucl/contexts.fna"
    output:
        out = "{dir}/blastout_position.txt",
        tmp = temp("{dir}/pos_tmp.txt"),
        positions = "{dir}/ARG_positions.txt"
    shell:
        """
        set +o pipefail

        blastn -query {input.genes} -db {input.contexts} -out {output.out} -evalue 0.00000001 -qcov_hsp_perc 100 -outfmt 7
        grep -v '^#' {output.out} | cut -f 2,9,10 > {output.tmp}
        python /home/dlund/scripts/miscellaneous_scripts/remove_duplicate_lines.py {output.tmp} {output.positions}
        """

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Search for conjugative elements
#----------------------------------------------------------------------------------------------------------------------------------------------------------
rule apply_conjscan:
    input:
        "{dir}/contexts_protein.fna"
    output:
        mge = temp("{dir}/conj_elements.txt"),
        hits = temp("{dir}/conj_hits.txt"),
        dir = directory("{dir}/conjscan")
    shell:
        """
        set +o pipefail
        mkdir {output.dir}

        for model in /home/dlund/models/conjscan_models/*
        do 
            name=$(echo $model | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1)
            hmmsearch --domtblout {output.dir}/$name -E 0.0000000001 $model {input}
        done

        for f in {output.dir}/*
        do 
            hit=$(grep -v '^#' $f | wc -l)
            if [[ $hit == 0 ]]; then 
                rm $f 
            else
                echo $hit >> {output.hits}
            fi 
        done

        if [[ ! -f {output.hits} ]]; then
            touch {output.hits}
        fi

        ls -l {output.dir} | tail -n +2 | rev | cut -d ' ' -f 1 | rev  > {output.mge}
        """

rule adjust_conjscan:
    input:
        mge = "{dir}/conj_elements.txt",
        hits = "{dir}/conj_hits.txt",
        list = "{dir}/assembly_accessions.txt"
    output:
        "{dir}/conj_sum.txt"
    run:
        file1 = open(input[0], 'r')
        mge = file1.readlines()
        file1.close()

        file2 = open(input[1], 'r')
        hits = file2.readlines()
        file2.close()

        file3 = open(input[2], 'r')
        list = file3.readlines()
        file3.close()

        sum = ""

        for i in range(len(mge)):
            if i == (len(mge)-1):
                sum += str(mge[i].rstrip()) + " (" + str(hits[i]).rstrip() + ")"
            else:
                sum += str(mge[i].rstrip()) + " (" + str(hits[i]).rstrip() + ")" + "\n"

        if len(mge) == 0:
            sum += "NA"

        outfile = open(output[0], 'w')
        outfile.write(sum)
        outfile.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Search for IS elements
#----------------------------------------------------------------------------------------------------------------------------------------------------------
rule is_blast:
    input:
        "{dir}/db_nucl/contexts.fna"
    output:
        out = "{dir}/blastout_is.txt"
    shell:
        """
        set +o pipefail

        blastn -query /home/dlund/index_files/ISfinder-sequences/IS.fna -out {output.out} -db {input} -evalue 0.001 -qcov_hsp_perc 50 -outfmt 7
        """

rule compile_is:
    input:
        results = "{dir}/blastout_is.txt",
        positions = "{dir}/ARG_positions.txt"
    output:
        temp("{dir}/is_elements.txt")
    run:      
        is_information = summarize_blast_results_distance_considered(input[0],input[1],90.00,1000)        

        is_output = ""

        if is_information is None:
            is_information = []

        for row in is_information:
            is_output += str(row)

        if is_output == "":
            is_output += "NA"

        outfile = open(output[0], 'w')
        outfile.write(is_output)
        outfile.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Search for integrons
#----------------------------------------------------------------------------------------------------------------------------------------------------------
rule apply_integron_finder:
    input:
        "{dir}/contexts.fna"
    output:
        sum = "{dir}/Results_Integron_Finder_contexts/contexts.summary",
        id = "{dir}/integron_id.txt",
        complete = "{dir}/complete.txt"
    conda:
        "/home/dlund/envs/IFv2.yaml"
    shell:
        """
        set +o pipefail

        size=$(grep '>' {input} | wc -l)

        if [[ $size -gt 100 ]]; then
            grep '>' {input} > {wildcards.dir}/lines.txt

            python /home/dlund/scripts/miscellaneous_scripts/get_random_entries.py {wildcards.dir}/lines.txt 100 {wildcards.dir}/sample.txt

            python /home/dlund/scripts/miscellaneous_scripts/extract_subset.py {wildcards.dir}/sample.txt {wildcards.dir}/lines.txt {input} {wildcards.dir}/context_subset.fna

            if [[ ! -f {output.sum} ]]; then
                integron_finder --outdir {wildcards.dir} --mute {wildcards.dir}/context_subset.fna

#                mv {wildcards.dir}/Results_Integron_Finder_context_subset/context_subset.summary {output.sum}
            fi

            if [[ ! -f {output.sum} ]]; then
                touch {output.sum}
            fi

            rm {wildcards.dir}/lines.txt {wildcards.dir}/sample.txt {wildcards.dir}/context_subset.fna

        else
	    if [[ ! -f {output.sum} ]]; then
                integron_finder --outdir {wildcards.dir} --mute {input}
            fi

            if [[ ! -f {output.sum} ]]; then
                touch {output.sum} 
            fi
        fi

        cat {wildcards.dir}/Results_Integron_Finder*/*.summary | cut -f 2 | tail -n +2 > {output.id}
        cat {wildcards.dir}/Results_Integron_Finder*/*.summary | cut -f 3 | tail -n +2 > {output.complete}
        """

rule compile_integrons:
    input:
        id = "{dir}/integron_id.txt",
        complete = "{dir}/complete.txt"
    output:
        "{dir}/integrons.txt"
    run:
        file1 = open(input[0], 'r')
        ids = file1.readlines()
        file1.close()

        file2 = open(input[1], 'r')
        complete = file2.readlines()
        file2.close()

        unique_names = []

        for name in ids:
            if name not in unique_names:
                unique_names.append(name)

        sum_complete = np.zeros((len(ids),len(unique_names)))

        for i in range(len(unique_names)):
            for j in range(len(ids)):
                if ids[j] == unique_names[i]:
                    sum_complete[j,i] = int(complete[j])

        sum_total = []

        for i in range(len(unique_names)):
            sum_total.append(unique_names[i].rstrip() + ": " + "Complete (" + str(sum(sum_complete[:,i])) + ")" + "\n")

        integrons = ""

        for row in sum_total:
            integrons += str(row)

        if integrons == "":
            integrons += "NA"

        outfile = open(output[0], 'w')
        outfile.write(integrons)
        outfile.close()
                    
#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Identify co-localized ARGs
#----------------------------------------------------------------------------------------------------------------------------------------------------------
rule blast_co_loc:
    input:
        "{dir}/db_nucl/contexts.fna"
    output:
        "{dir}/blastout_co_loc.txt"
    shell:
        """
        blastn -query /home/dlund/index_files/resfinder.fsa -db {input} -out {output} -evalue 0.001 -qcov_hsp_perc 75 -outfmt 7
        """

rule compile_co_loc:
    input:
        results = "{dir}/blastout_co_loc.txt",
        positions = "{dir}/ARG_positions.txt" 
    output:
        temp("{dir}/co_loc_ARGs.txt")
    run:
        info_co_loc = summarize_blast_results(input[0],90.00)

        output_args = ""

        if info_co_loc is None:
            info_co_loc = []

        for row in info_co_loc:
            output_args += str(row) 

        if output_args == "":
            output_args += "NA"

        outfile = open(output[0], 'w')
        outfile.write(output_args)
        outfile.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Compile information about the identifed hosts of each gene family
#----------------------------------------------------------------------------------------------------------------------------------------------------------
rule compile_taxonomy:
    input:
        species = "{dir}/species.txt",
        pathogens = "{dir}/pathogens.txt",
        phylum = "{dir}/phylum.txt"
    output:
        species = temp("{dir}/top_species.txt"),
        pathogens = temp("{dir}/top_pathogens.txt"),
        phylum = temp("{dir}/top_phylum.txt")
    run:
        infile1 = open(input[0], 'r')
        species = infile1.readlines()
        infile1.close()

        infile2 = open(input[1], 'r')
        pathogens = infile2.readlines()
        infile2.close()

        infile3 = open(input[2], 'r')
        phylum = infile3.readlines()
        infile3.close()

        top_species = get_top_items(species)
        top_pathogens = get_top_items(pathogens)
        top_phylum = get_top_items(phylum)

        outfile1 = open(output[0], 'w')
        outfile1.write(str(top_species))
        outfile1.close()

        outfile2 = open(output[1], 'w')
        outfile2.write(str(top_pathogens))
        outfile2.close()

        outfile3 = open(output[2], 'w')
        outfile3.write(str(top_phylum))
        outfile3.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Find similarity to known genes
#----------------------------------------------------------------------------------------------------------------------------------------------------------
rule blast:
    input:
        "{dir}/proteins.fna"
    output:
        "{dir}/blastout.txt"
    shell:
        """
        set +o pipefail
        blastp -query {input} -out {output} -db /home/dlund/index_files/blast_db/reference_ARGs.fsa -qcov_hsp_perc 50 -outfmt 7 -max_target_seqs 1
        """

rule adjust_blast:
    input:
        "{dir}/blastout.txt"
    output:
        temp("{dir}/blast_info.txt")
    run:
        subject = []
        identity = []
        summary = []

        with open(input[0]) as tsv:
            for line in csv.reader(tsv, delimiter="\t"):
                hit = re.search("#", str(line))
                if hit is None:
                    subject.append(line[1].replace("'","").replace("(","").replace(")",""))
                    identity.append(line[2])

        unique_genes = []
        
        for name in subject:
            if name not in unique_genes:
                unique_genes.append(name)

        for name in unique_genes:
            percentages = []

            for i in range(len(subject)):
                hit = re.search(name, subject[i])
                if hit is not None:
                    percentages.append(float(identity[i]))

            if min(percentages) == max(percentages):
                summary.append(name + " [" + str(max(percentages)) + "%]" + "\n")
            else:
                summary.append(name + " [" + str(min(percentages)) + "-" + str(max(percentages)) + "%]" + "\n")


        blast_info = ""
        
        for row in summary:
            blast_info += str(row)

        if blast_info == "":
            blast_info += "NA"

        outfile = open(output[0], 'w')
        outfile.write(blast_info)
        outfile.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Compile metadata into one file describing each gene family, and combine these into a table describing all families
#----------------------------------------------------------------------------------------------------------------------------------------------------------
rule compile_metadata:
    input:
        conj = "{dir}/conj_sum.txt",
        iss = "{dir}/is_elements.txt",
        integrons = "{dir}/integrons.txt",
        co_loc = "{dir}/co_loc_ARGs.txt",
        blast ="{dir}/blast_info.txt",
        list = "{dir}/assembly_accessions.txt",
        species = "{dir}/top_species.txt",
        pathogens = "{dir}/top_pathogens.txt",
        phylum = "{dir}/top_phylum.txt"
    output:
        data = "{dir}/metadata.txt",
        name = temp("{dir}/name.txt"),
        number = temp("{dir}/number.txt")
    shell:
        """
        echo {wildcards.dir} > {output.name}
        cat {input.list} | wc -l > {output.number}
        paste -d '\t' {output.name} {output.number} {input.blast} {input.phylum} {input.species} {input.pathogens} {input.conj} {input.iss} {input.integrons} {input.co_loc} > {output.data}
        echo '\n' >> {output.data}
        """

rule make_table_head:
    output:
        temp("head.txt")
    run:
        head = ""
        head += 'Family' + '\t' + 'Genes' + '\t' + 'Closest known' + '\t' + 'Phylum' + '\t' + 'Species' + '\t' + 'Pathogens' + '\t' + 'Conjugative elements' + '\t' + 'IS elements' + '\t' + 'Integrons' + '\t' + 'Co-localized ARGs' +'\n' + '\n'

        file = open(output[0], 'w')
        file.write(head)
        file.close()

rule make_table:
    input:
        data=expand("{dir}/metadata.txt", dir=DIR),
        head="head.txt"
    output:
        "metadata_complete.txt"
    shell:
        "cat {input.head} */metadata.txt > {output}"

#----------------------------------------------------------------------------------------------------------------------------------------------------------