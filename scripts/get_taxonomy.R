#!/usr/bin/env/Rscript
library(taxonomizr)

# Prepare database
if (!(file.exists("databases/accessionTaxa.sql"))) {
    prepareDatabase("databases/accessionTaxa.sql") 
}

# Convert NCBI accession IDs to full lineage
accno <- system(paste("grep '>' predicted-orfs-amino.fasta | cut -d '_' -f 2"), intern=TRUE)
accno <- unique(accno)

id <- accessionToTaxa(accno,"databases/accessionTaxa.sql") 
taxonomy <- data.frame(getTaxonomy(id,"databases/accessionTaxa.sql"))
rownames(taxonomy) <- accno

write.table(taxonomy, "index_files/taxonomy.txt",sep="\t", quote=FALSE, col.names=FALSE) 
