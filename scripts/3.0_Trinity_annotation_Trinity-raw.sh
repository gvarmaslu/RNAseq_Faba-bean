#!/bin/bash -l


# De novo RNA-Seq Assembly, Annotation, and Analysis Using Trinity and Trinotate
#https://bioinformaticsdotca.github.io/RNAseq_2019_Module8_lab

#######################
########
#for Trinity
TRINITY_HOME="/bioinfo/trinityrnaseq-Trinity-v2.8.5/"
WORKDIR="/home/gala0002/proj/proj_dir/"
cd ${WORKDIR}
mkdir ${WORKDIR}3.1_align_and_estimate_abundances_Trinity
#######################
######
#Transcript expression quantitation using Salmon
######

${TRINITY_HOME}util/align_and_estimate_abundance.pl --seqType fq \
 --samples_file ${WORKDIR}sampleinfo_all_files.txt --transcripts ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta \
 --est_method salmon --trinity_mode --prep_reference --thread_count 70 --output_dir ${WORKDIR}3.1_align_and_estimate_abundances_Trinity


######
find 3.1_align_and_estimate_abundances_Trinity/* -name "quant.sf" | tee quant_files.list 

######
${TRINITY_HOME}util/abundance_estimates_to_matrix.pl --est_method salmon \
--out_prefix Trinity --name_sample_by_basedir \
--quant_files quant_files.list \
--gene_trans_map ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta.gene_trans_map


######
head -n20 Trinity.isoform.counts.matrix | column -t
head -n20 Trinity.isoform.TMM.EXPR.matrix | column -t


#######################
######
#Functional Annotation of Assembled Transcripts Using Trinotate
######

mkdir ${WORKDIR}4.1_Trinotate
cd ${WORKDIR}4.1_Trinotate

TRANSDECODER_HOME="/bioinfo/TransDecoder-TransDecoder-v5.5.0/"

#Identification of likely protein-coding regions in transcripts

$TRANSDECODER_HOME/TransDecoder.LongOrfs -t ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta
$TRANSDECODER_HOME/TransDecoder.Predict -t ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta


#######
#Identification of likely protein-coding regions in transcripts

$TRANSDECODER_HOME/TransDecoder.LongOrfs -t ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta
$TRANSDECODER_HOME/TransDecoder.Predict -t ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta


######
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4028052/
#TGI Clustering tools (TGICL): a software system for fast clustering of large EST datasets

cd ${WORKDIR}4.1_Trinotate

TGICL="/bioinfo/trans-assembly-prune/TGICL-2.1/bin/"

#perl5.16
#/localperl/bin/

${TGICL}tgicl ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta -p 99 -l 50 -v 100 > ${WORKDIR}4.1_Trinotate/TGICL_output.txt


######
#CAP3
#http://seq.cs.iastate.edu/cap3.html

CAP3="/bioinfo/trans-assembly-prune/CAP3/"

${CAP3}cap3 ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta > ${WORKDIR}4.1_Trinotate/cap3_output.txt

#######################
#######
#Sequence homology searches

#BLAST
#https://biohpc.cornell.edu/lab/doc/Trinity_exercise2.pdf

NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin/"

BLASTDB="/mnt/udda-backup/data/data_backup_2020-Jun/DATABASE/"
HMMER="/bioinfo/HMMER/hmmer-3.3/src/"

WORKDIR="/home/gala0002/proj/proj_dir/"

######################
######-- uniprot_sprot
${NCBIBLAST}blastp -query ${WORKDIR}4.1_Trinotate/Trinity.fasta.transdecoder.pep \
-db ${BLASTDB}uniprot_sprot.fasta -num_threads 65 \
-max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/Trinity.fasta.transdecoder.pep_swissprot.blastp.outfmt6

######################
###### --Glycine-max_soybean

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein.faa -out /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein -input_type fasta -dbtype prot

######-- search db
${NCBIBLAST}blastp -query ${WORKDIR}4.1_Trinotate/Trinity.fasta.transdecoder.pep \
-db /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein -num_threads 65 \
-max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/Trinity.fasta.transdecoder.pep_Glycine-max.blastp.outfmt6

######################
###### --Medicago-truncatula

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_protein.faa -out /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_protein -input_type fasta -dbtype prot

######-- search db
${NCBIBLAST}blastp -query ${WORKDIR}4.1_Trinotate/Trinity.fasta.transdecoder.pep \
-db /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_protein -num_threads 65 \
-max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/Trinity.fasta.transdecoder.pep_Medicago-truncatula.blastp.outfmt6

######################
###### --Pisum_sativum

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_prot.fa -out /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_prot -input_type fasta -dbtype prot

######-- search db
${NCBIBLAST}blastp -query ${WORKDIR}4.1_Trinotate/Trinity.fasta.transdecoder.pep \
-db /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_prot -num_threads 65 \
-max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/Trinity.fasta.transdecoder.pep_Pisum-sativum.blastp.outfmt6


######################
###### --Arabidopsis-thaliana

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.3_TAIR10_protein.faa -out /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.3_TAIR10_protein -input_type fasta -dbtype prot

######-- search db
${NCBIBLAST}blastp -query ${WORKDIR}4.1_Trinotate/Trinity.fasta.transdecoder.pep \
-db /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.3_TAIR10_protein -num_threads 65 \
-max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/Trinity.fasta.transdecoder.pep_Arabidopsis-thaliana.blastp.outfmt6


######

${HMMER}hmmscan --cpu 65 --domtblout ${WORKDIR}4.2_blastout/Trinity.fasta.transdecoder.pep_TrinotatePFAM.out \
${BLASTDB}Pfam-A.hmm \
${WORKDIR}4.1_Trinotate/Trinity.fasta.transdecoder.pep


#####################################################################

#######################
######
#Computational prediction of sequence features
#signalp -f short -n signalp.out Trinity.fasta.transdecoder.pep








