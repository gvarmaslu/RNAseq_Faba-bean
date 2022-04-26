#!/bin/bash -l


# TF Annotation
#######################

#####################################################################
#######
#Sequence homology searches

#BLAST
#https://biohpc.cornell.edu/lab/doc/Trinity_exercise2.pdf

NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin/"

WORKDIR="/home/gala0002/proj/proj_dir/"

# PlantTFDBv5.0
######################

#create DB
$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/Glycine_Max_TF/PlantTFDBv5.0/Gma_pep_pars.fas -out /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/Glycine_Max_TF/PlantTFDBv5.0/Gma_pep -input_type fasta -dbtype prot
#create DB
$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/Medicago_Truncatula_TF/PlantTFDBv5.0/Mtr_pep_pars.fas -out /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/Medicago_Truncatula_TF/PlantTFDBv5.0/Mtr_pep -input_type fasta -dbtype prot
#create DB
$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Arabidopsis-thaliana/PlantTFDBv5.0/Ath_pep_pars.fas -out /home/gala0002/proj/proj_dir/REF/Arabidopsis-thaliana/PlantTFDBv5.0/Ath_pep -input_type fasta -dbtype prot

######################
###### --Glycine-max_soybean
######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/Glycine_Max_TF/PlantTFDBv5.0/Gma_pep -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_PlantTFDBv5.0_Glycine-max.blastx.outfmt6

######################
###### --Medicago-truncatula

######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/Medicago_Truncatula_TF/PlantTFDBv5.0/Mtr_pep -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_PlantTFDBv5.0_Medicago-truncatula.blastx.outfmt6

######################
###### --Arabidopsis-thaliana

######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF/Arabidopsis-thaliana/PlantTFDBv5.0/Ath_pep -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_PlantTFDBv5.0_Arabidopsis-thaliana.blastx.outfmt6


######################
###### --Pisum_sativum



###############################################################


#iTAKv18.12

#cat Thale_cress-transcription-factor.fa | awk -F"\t" '{print ">"$1"_"$2"\n"$3}' > Thale_cress-transcription-factor.fasta

#cat 'Barrel_medic-transcription factor.txt' 'Barrel_medic-transcriptional regulator.txt' | awk -F"\t" '{print ">"$1"_"$2"\n"$3}' > Barrel_medic_TF-TR_prot.fa
#cat Soybean-transcription_factor.txt Soybean-transcriptional_regulator.txt | awk -F"\t" '{print ">"$1"_"$2"\n"$3}' > Soybean_TF-TR_prot.fa


######
#BLAST
#https://biohpc.cornell.edu/lab/doc/Trinity_exercise2.pdf

NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin/"

WORKDIR="/home/gala0002/proj/proj_dir/"

# PlantTFDBv5.0
######################
######################
###### --Glycine-max_soybean

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/Glycine_Max_TF/iTAKv18.12/Soybean_TF-TR_prot.fa -out /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/Glycine_Max_TF/iTAKv18.12/Soybean_TF-TR_prot -input_type fasta -dbtype prot
#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/Medicago_Truncatula_TF/iTAKv18.12/Barrel_medic_TF-TR_prot.fa -out /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/Medicago_Truncatula_TF/iTAKv18.12/Barrel_medic_TF-TR_prot -input_type fasta -dbtype prot
#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Arabidopsis-thaliana/iTAKv18.12/Thale_cress_TF-TR_prot.fa -out /home/gala0002/proj/proj_dir/REF/Arabidopsis-thaliana/iTAKv18.12/Thale_cress_TF-TR_prot -input_type fasta -dbtype prot


######################
###### --Glycine-max_soybean
######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/Glycine_Max_TF/iTAKv18.12/Soybean_TF-TR_prot -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_iTAKv18.12_Glycine-max.blastx.outfmt6

######################
###### --Medicago-truncatula

######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/Medicago_Truncatula_TF/iTAKv18.12/Barrel_medic_TF-TR_prot -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_iTAKv18.12_Medicago-truncatula.blastx.outfmt6

######################
###### --Arabidopsis-thaliana

######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF/Arabidopsis-thaliana/iTAKv18.12/Thale_cress_TF-TR_prot -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_iTAKv18.12_Arabidopsis-thaliana.blastx.outfmt6

######################
###### --Pisum_sativum

###############################################################



