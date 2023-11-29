#!/bin/bash

# Making a reduced reference genome for black spruce (Pmarianav1), based on BLAST'ing the fasta sequences of the baits used in Capblancq et al. (2020) based on the Picea glauca genome

# Make a BLAST db for the black spruce reference genome (only need to do this once, and then comment out)
# makeblastdb -in ~/mydata/datashare/ReferenceGenomes/Picea/Pmarianav1/assembly/GCA_032191745.1_Pmariana_40-10-1_v1_genomic.fna -dbtype nucl -parse_seqids -blastdb_version 5

cd ~/mydata/datashare/Spruce

# Blasting the baits sequences against the norway spruce genome 
/popgen/blast/ncbi-blast-2.7.1+/bin/blastn \
-query ./exome_capture/WES_mapping/Reference_BaitsWES/UVM_131901_RG_0806_probes_ref_WS77111_Pglauca_Scaffolds.fasta \
-db ../ReferenceGenomes/Picea/Pmarianav1/assembly/GCA_032191745.1_Pmariana_40-10-1_v1_genomic.fna \
-task blastn -num_threads 15 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
-out ~/blast_baits_pmariana.out -max_target_seqs 1
#-out ./exome_capture/WES_mapping/Reference_BaitsWES/blast_baits_pmariana.out -max_target_seqs 1


# Listing the scaffolds of P. mariana genome where the baits mapped
cat ~/blast_baits_pmariana.out | awk '{print $2}' | sort | uniq > ~/list_scaffolds_pmariana.txt
#cat ./exome_capture/WES_mapping/Reference_BaitsWES/blast_baits_pmariana.out | awk '{print $2}' | sort | uniq > ./exome_capture/WES_mapping/Reference_BaitsWES/list_scaffolds_pmariana.txt

# Reducing the P. mariana reference with only the matching scaffolds
samtools faidx ../ReferenceGenomes/Picea/Pmarianav1/assembly/GCA_032191745.1_Pmariana_40-10-1_v1_genomic.fna `cat ~/list_scaffolds_pmariana.txt` > ~/Pmariana1.0-genome_reduced.fa
#samtools faidx ../ReferenceGenomes/Picea/Pmarianav1/assembly/GCA_032191745.1_Pmariana_40-10-1_v1_genomic.fna `cat ./exome_capture/WES_mapping/Reference_BaitsWES/list_scaffolds_pmariana.txt` > ./exome_capture/WES_mapping/ReferenceGenomes/Pmariana1.0-genome_reduced.fa

# Number of scaffolds retained and total length 
cat ~/Pmariana1.0-genome_reduced.fa | grep '>' - | wc -l # how many scaffolds? 28,953 
grep -v ">" ~/Pmariana1.0-genome_reduced.fa | wc | awk '{print $3-$1}' # how many bp? 2,927,599,912


#cat ./exome_capture/WES_mapping/ReferenceGenomes/Pmariana1.0-genome_reduced.fa | grep '>' - | wc -l # how many scaffolds? 
#grep -v ">" ./exome_capture/WES_mapping/ReferenceGenomes/Pmariana1.0-genome_reduced.fa | wc | awk '{print $3-$1}' # how many bp?


