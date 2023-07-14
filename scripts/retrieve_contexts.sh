#!/bin/bash
for dir in clusters/cluster* 
do
	# Retrieve genetic contexts of args within each cluster
    hit=$(cat $dir/assembly_accessions.txt | wc -l)
    if [[ ! -f $dir/unique_proteins.dmnd && $hit != 0 ]]; then

		python genview_dir/GEnView/genview_create_db.py -d genview_dir/db -db $dir/unique_proteins.fna -p 10 -id 100 -scov 100 --acc_list $dir/assembly_accessions.txt --split 5 --assemblies
		
		python genview_dir/GEnView/genview_extract.py -genes 'arg' -o context_dir -db genview_dir/db/context_db_flank.db -id 100

	    mv context_dir/arg_100_analysis/arg_contexts.fna $dir/contexts.fna

		rm genview_dir/db/GCA_* genview_dir/db/orfs_clustered* genview_dir/db/all* genview_dir/db/context_db_flank.db
		rm -r genview_dir/db/integrons_tmp context_dir
    fi
done
