#!/bin/bash -l


# De novo RNA-Seq Assembly, Annotation, and Analysis Using Trinity and Trinotate
#https://bioinformaticsdotca.github.io/RNAseq_2019_Module8_lab

#######################
########
#for Trinity
TRINITY_HOME="/bioinfo/trinityrnaseq-Trinity-v2.8.5/"
WORKDIR="/home/gala0002/proj/proj_dir/"
cd ${WORKDIR}
mkdir ${WORKDIR}3.2_align_and_estimate_abundances_trinity_GSM
#######################
######
#Transcript expression quantitation using Salmon
######

${TRINITY_HOME}util/align_and_estimate_abundance.pl --seqType fq \
 --samples_file ${WORKDIR}sampleinfo_all_files.txt --transcripts ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
 --est_method salmon --trinity_mode --prep_reference --thread_count 70 --output_dir ${WORKDIR}3.2_align_and_estimate_abundances_trinity_GSM


######
find 3.2_align_and_estimate_abundances_trinity_GSM/* -name "quant.sf" | tee quant_files_trinity_GSM.list 

######
${TRINITY_HOME}util/abundance_estimates_to_matrix.pl --est_method salmon \
--out_prefix trinity_GSM --name_sample_by_basedir \
--quant_files quant_files_trinity_GSM.list \
--gene_trans_map none

#--gene_trans_map ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta.gene_trans_map


######
head -n20 Trinity.isoform.counts.matrix | column -t
head -n20 Trinity.isoform.TMM.EXPR.matrix | column -t


#########
# Function to keep transcripts > 2 reads (sum of transcript reads should be >2 and each sample should have at east 2 reads)

# Test case: remove 0 sum or  <1 rows:
#cat genefile.txt | awk '{f=0; for(i=2;i<=NF;i++) if($i+0!=0 && $i+0>=1)f=1} f'
#cat genefile.txt | awk '{for(i=2;i<=NF;i++) if($i >= 3){print; next}}' 

### 
#cat trinity_GSM.isoform.counts.matrix | awk '{f=0; for(i=2;i<=NF;i++) if($i+0!=0 && $i+0>=2)f=2} f' | wc -l 

#########
# Function to round decimal values  
cat trinity_GSM.isoform.counts.matrix | perl -pe 's/(\d*\.\d*)/int($1+0.5)/ge' > trinity_GSM.isoform.counts_clean_round.matrix

#cat trinity_GSM.isoform.counts.matrix | awk '{f=0; for(i=2;i<=NF;i++) if($i+0!=0 && $i+0>=2)f=2} f' | perl -pe 's/(\d*\.\d*)/int($1+0.5)/ge' > trinity_GSM.isoform.counts_filter_round.matrix

# Fitered out transcripts (unique in file1)
#comm -23 <(cat file1) <(cat file2) 
comm -23 <(cat trinity_GSM.isoform.counts.matrix | cut -f1 | sort -u ) <(cat trinity_GSM.isoform.counts_filter_round.matrix | cut -f1 | sort -u) | wc -l 
#30952
#comm -23 <(cat trinity_GSM.isoform.counts.matrix | cut -f1 | sort -u ) <(cat trinity_GSM.isoform.counts_clean_round.matrix | cut -f1 | sort -u) > trinity_GSM.isoform.counts_clean_round.matrix_uniq
#grep -f trinity_GSM.isoform.counts_filter_round.matrix_uniq -A 1 2.0_trinity_out_all_sampleinfo/Trinity.fasta  > trinity_GSM.isoform.counts_filter_round.matrix_uniq.fa
#grep -f trinity_GSM.isoform.counts_filter_round.matrix_uniq 4.2_blastout/trinity_gene_splice_modeler.fasta.transdecoder.pep_Medicago-truncatula.blastp.outfmt6_anno.tsv

#######################
######
#Identification of likely protein-coding regions in transcripts

#TransDecoder = Functional Annotation of Assembled Transcripts 
######

mkdir ${WORKDIR}4.1_Trinotate
cd ${WORKDIR}4.1_Trinotate

TRANSDECODER_HOME="/bioinfo/TransDecoder-TransDecoder-v5.5.0/"

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


######
#BORF
#https://github.com/betsig/borf

#pip install borf
WORKDIR="/home/gala0002/proj/proj_dir/"

mkdir ${WORKDIR}4.1.2_BORF
cd ${WORKDIR}4.1.2_BORF

borf ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta --output_path ${WORKDIR}4.1.2_BORF

#####################################################################
#######
#Sequence homology searches

#BLAST
#https://biohpc.cornell.edu/lab/doc/Trinity_exercise2.pdf

NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin/"

BLASTDB="/mnt/udda-backup/data/data_backup_2020-Jun/DATABASE/"
HMMER="/bioinfo/HMMER/hmmer-3.3/src/"

WORKDIR="/home/gala0002/proj/proj_dir/"

######################

######################
###### --Glycine-max_soybean

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein.faa -out /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein -input_type fasta -dbtype prot

######-- search db
nice -n 5 ${NCBIBLAST}blastp -query ${WORKDIR}4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep \
-db /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein -num_threads 65 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta.transdecoder.pep_Glycine-max.blastp.outfmt6

######################
###### --Medicago-truncatula

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_protein.faa -out /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_protein -input_type fasta -dbtype prot

######-- search db
nice -n 5 ${NCBIBLAST}blastp -query ${WORKDIR}4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep \
-db /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_protein -num_threads 65 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta.transdecoder.pep_Medicago-truncatula.blastp.outfmt6

######################
###### --Pisum_sativum

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_prot.fa -out /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_prot -input_type fasta -dbtype prot

######-- search db
nice -n 5 ${NCBIBLAST}blastp -query ${WORKDIR}4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep \
-db /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_prot -num_threads 65 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta.transdecoder.pep_Pisum-sativum.blastp.outfmt6


######################
###### --Arabidopsis-thaliana

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.3_TAIR10_protein.faa -out /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.3_TAIR10_protein -input_type fasta -dbtype prot

######-- search db
nice -n 5 ${NCBIBLAST}blastp -query ${WORKDIR}4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep \
-db /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.3_TAIR10_protein -num_threads 65 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta.transdecoder.pep_Arabidopsis-thaliana.blastp.outfmt6


#######################
######
#Computational prediction of sequence features
#signalp -f short -n signalp.out Trinity.fasta.transdecoder.pep

#####################################################################
#######
#Sequence homology searches


#BLAST
#https://biohpc.cornell.edu/lab/doc/Trinity_exercise2.pdf

NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin/"

BLASTDB="/mnt/udda-backup/data/data_backup_2020-Jun/DATABASE/"

WORKDIR="/home/gala0002/proj/proj_dir/"

######################
######################
###### --Glycine-max_soybean

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein.faa -out /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein -input_type fasta -dbtype prot

######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_Glycine-max.blastx.outfmt6

######################
###### --Medicago-truncatula

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_protein.faa -out /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_protein -input_type fasta -dbtype prot

######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_protein -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_Medicago-truncatula.blastx.outfmt6

######################
###### --Pisum_sativum

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_prot.fa -out /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_prot -input_type fasta -dbtype prot

######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_prot -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_Pisum-sativum.blastx.outfmt6


######################
###### --Arabidopsis-thaliana

#create DB
#$NCBIBLAST/makeblastdb -in /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.3_TAIR10_protein.faa -out /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.3_TAIR10_protein -input_type fasta -dbtype prot

######-- search db
nice -n 5 ${NCBIBLAST}blastx -query ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
-db /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.3_TAIR10_protein -num_threads 68 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 \
> ${WORKDIR}4.2_blastout/trinity_gene_splice_modeler.fasta_Arabidopsis-thaliana.blastx.outfmt6

###############################################################


#######################
######
#https://bioinformaticsdotca.github.io/RNAseq_2019_Module8_lab
#Functional Annotation of Assembled Transcripts Using Trinotate
######
#Preparing and Generating a Trinotate Annotation Report
#Preparing Trinotate (loading the database)

#cp ../data/trinotate_data/Trinotate.boilerplate.sqlite  Trinotate.sqlite

#cp ../data/trinotate_data/Trinotate.boilerplate.sqlite  Trinotate.sqlite
#wget -c https://data.broadinstitute.org/Trinity/__deprecated_trinotate_resources/Trinotate_v2.0_RESOURCES/Trinotate.sprot_uniref90.20150131.boilerplate.sqlite.gz

TRINITY_HOME="/bioinfo/trinityrnaseq-Trinity-v2.8.5/"
WORKDIR="/home/gala0002/proj/proj_dir/"
cd ${WORKDIR}4.1.2_Trinotate/

#Load your Trinotate.sqlite database with your Trinity transcripts and predicted protein sequences:

#/bioinfo/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite init \
#--transcript_fasta  ${WORKDIR}2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta \
#--transdecoder_pep ${WORKDIR}4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep


/bioinfo/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite init \
--gene_trans_map ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta.gene_trans_map \
--transcript_fasta ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta \
--transdecoder_pep ${WORKDIR}4.1_Trinotate/Trinity.fasta.transdecoder.pep


#Load in the various outputs generated earlier:

/bioinfo/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_pfam ${WORKDIR}4.2_blastout/Trinity.fasta.transdecoder.pep_TrinotatePFAM.out

#Generate the Trinotate Annotation Report

/bioinfo/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite report > Trinotate.xls


###############################################################

###############################################################


# annotation 

#Glycine-max_soybean

head trinity_gene_splice_modeler.fasta_Glycine-max.blastx.outfmt6 | cut -f2 > blast-target-out.txt

grep -f blast-target-out.txt /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_protein.faa -w | awk '{sub(/ /,"\t")}1' | cut -f2 | cut -d"[" -f 1


#######################
######
#https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats
#Trinity Transcriptome Contig Nx and ExN50 Statistics
TRINITY_HOME="/bioinfo/trinityrnaseq-Trinity-v2.6.5"


$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl \
	${WORKDIR}Trinity.isoform.TMM.EXPR.matrix ${WORKDIR}2.0_trinity_out_all_sampleinfo/Trinity.fasta | tee ExN50.stats



${TRINITY_HOME}/util/misc/plot_ExN50_statistic.Rscript  ExN50.stats  


#xpdf ExN50.stats.plot.pdf


#cat Trinity.isoform.TMM.EXPR.matrix.E-inputs |  egrep -v ^\# | awk '$1 <= 96' | wc -l







