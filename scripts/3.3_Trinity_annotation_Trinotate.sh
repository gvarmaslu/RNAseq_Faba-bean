#!/bin/bash -l


# De novo RNA-Seq Assembly, Annotation, and Analysis Using Trinity and Trinotate
#Transcriptome annotation: Workshop on RNA-Seq
#https://nbisweden.github.io/workshop-RNAseq/1906/lab_annotation.html
#http://arcsb2017.weebly.com/uploads/4/6/1/9/46199061/trinotate_install.pdf

#######################
########
#for Trinity
WORKDIR="/home/varma/proj/proj_faba-bean/"
cd ${WORKDIR}

####

cd /home/varma/proj/proj_faba-bean/
mkdir RNAseq_annotation
cd RNAseq_annotation

####
mkdir Trinotate
ln -s /home/varma/proj/proj_faba-bean/2.0_trinity_out_all_sampleinfo/Trinity_200.fa

#

cd /home/gala0002/proj/proj_dir/4.1.2_Trinotate
qfna="../2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta"
qfaa="../4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep"

cat $qfna \
    | grep  '>' \
    | cut -d' ' -f1 \
    | cut -c2- \
    | awk '{print $1 "\t" $1}' > tmp.map
mv tmp.map transcripts.gene2tr.map

join <(cat transcripts.gene2tr.map | sort -k1,1b) <(cat $qfna | awk '{if(/>/){printf("\n%s\t",$1)}else{printf("%s",$1)}}' | cut -c2- | sort -k1,1b) | head


####
module load bioinfo-tools
module load trinotate rnammer/1.2 	

/sw/apps/bioinfo/trinotate/3.1.1/rackham/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
#/bioinfo/Trinotate-Trinotate-v3.1.1/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

makeblastdb -in uniprot_sprot.pep -dbtype prot

####

wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

####
#1.2 Determining longest Open Reading Frames (ORF)

module load bioinfo-tools

module load trinotate
TransDecoder.LongOrfs -t Trinity_200.fa
TransDecoder.Predict -t Trinity_200.fa

####
#1.3 BLAST approach
#1.3.1 Perform Blast searches from the command line on Uppmax:

blastx -db ~/RNAseq_assembly_annotation/assembly_annotation/database/uniprot_dmel/uniprot_dmel.fa -query Trinity_200.fa -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > swissprot.blastx.outfmt6
blastp -db ~/RNAseq_assembly_annotation/assembly_annotation/database/uniprot_dmel/uniprot_dmel.fa -query Trinity_200.fa.transdecoder.pep -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > swissprot.blastp.outfmt6

#
blastx -db uniprot_sprot.pep -query trinity_gene_splice_modeler.fasta -max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 -num_threads 68 > transdecoder.pep.swissprot.blastx.outfmt6
blastp -db uniprot_sprot.pep -query ../4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep -max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 -num_threads 68 > transdecoder.pep.swissprot.blastp.outfmt6

#
NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin/"
BLASTDB="/mnt/udda-backup/data/data_backup_2020-Jun/DATABASE/"
#${NCBIBLAST}makeblastdb -in ${BLASTDB}uniprot_sprot.fa -dbtype prot

#blastx -db ${BLASTDB}uniprot_sprot.fa -query ../2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta -max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 -num_threads 68 > swissprot.blastx.outfmt6
blastp -db ${BLASTDB}uniprot_sprot.fa -query ../4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.peptide -max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-5 -num_threads 68 > swissprot.blastp.outfmt6


####
#1.4 Domain/profiles homology searches
#1.4.1 Pfam database
#hmmscan --cpu 5 --domtblout TrinotatePFAM.out ~/RNAseq_assembly_annotation/assembly_annotation/database/trinotate_database/Pfam-A.hmm Trinity_200.fa.transdecoder.pep

#hmmscan --cpu 6 --domtblout TrinotatePFAM.out Pfam-A.hmm ../4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep

hmmscan --cpu 5 --domtblout TrinotatePFAM.out Pfam-A.hmm trinity_gene_splice_modeler.fasta.transdecoder.pep
####
#1.4.2 Signal peptide
export SNIC_TMP=$SNIC_TMP:~/RNAseq_annotation/trinotate
export SNIC_TMP=$SNIC_TMP:~/RNAseq_annotation/bioinfo/Trinotate-Trinotate-v3.1.1/Trinotate

signalp -f short -n signalp.out Trinity_200.fa.transdecoder.pep


/bioinfo/signalp/signalp-4.1/signalp -f short -n transdecoder.pep.signalp.out ../4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep

####
#1.4.3 Transmembrane region
tmhmm --short < Trinity_200.fa.transdecoder.pep > tmhmm.out

/bioinfo/tmhmm-2.0c/bin/tmhmm --short < ../4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep > transdecoder.pep.tmhmm.out


####
#1.4.4 rnammer
rnammer=/bioinfo/rnammer/rnammer
trinotate_dir="/bioinfo/Trinotate-Trinotate-v3.1.1"
$trinotate_dir/util/rnammer_support/RnammerTranscriptome.pl --transcriptome ../2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta --path_to_rnammer $rnammer &>rnammer.log &

# 
module load rnammer/1.2-fasrc01
module load trinotate

/sw/apps/bioinfo/trinotate/3.1.1/rackham/util/rnammer_support/RnammerTranscriptome.pl --transcriptome trinity_gene_splice_modeler.fasta --path_to_rnammer /sw/bioinfo/rnammer/1.2/rackham/rnammer

####
#1.5 Load results into the database

cp ~/RNAseq_assembly_annotation/assembly_annotation/database/trinotate_database/Trinotate.sqlite .

/sw/apps/bioinfo/trinity/2.8.2/rackham/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity_200.fa > Trinity_200.fa.gene_trans_map

Trinotate Trinotate.sqlite init --gene_trans_map Trinity_200.fa.gene_trans_map --transcript_fasta Trinity_200.fa --transdecoder_pep Trinity_200.fa.transdecoder.pep

Trinotate Trinotate.sqlite LOAD_swissprot_blastp transdecoder.pep.swissprot.blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx swissprot.blastx.outfmt6

Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.out

#
#prepare the database
#/bioinfo/Trinotate-Trinotate-v3.2.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

TRINITY_HOME="/bioinfo/trinityrnaseq-Trinity-v2.8.5/"
Trinotate="/bioinfo/Trinotate-Trinotate-v3.2.2/Trinotate"
Trinotate_dir="/bioinfo/Trinotate-Trinotate-v3.2.2"

#${TRINITY_HOME}/util/support_scripts/get_Trinity_gene_to_trans_map.pl ../2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta > trinity_gene_splice_modeler.fasta.gene_trans_map

${Trinotate} Trinotate.sqlite init --gene_trans_map transcripts.gene2tr.map --transcript_fasta ../2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta --transdecoder_pep ../4.1_Trinotate/trinity_gene_splice_modeler.fasta.transdecoder.pep	

#Load to database

${Trinotate} Trinotate.sqlite LOAD_swissprot_blastp transdecoder.pep.swissprot.blastp.outfmt6
${Trinotate} Trinotate.sqlite LOAD_swissprot_blastx transdecoder.pep.swissprot.blastx.outfmt6
${Trinotate} Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
${Trinotate} Trinotate.sqlite LOAD_tmhmm transdecoder.pep.tmhmm.out
${Trinotate} Trinotate.sqlite LOAD_signalp transdecoder.pep.signalp.out

####
#1.6 Output and Annotation Report
#Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls	
# --pfam_cutoff DGC 

${Trinotate} Trinotate.sqlite report > trinotate_annotation_report.xls


#### feature_name_encoding_attributes
${Trinotate_dir}/util/Trinotate_get_feature_name_encoding_attributes.pl trinotate_annotation_report.xls > trinotate_annotation_report_feature_map.txt


####
#GO Annotation (optional)
${Trinotate_dir}/util/extract_GO_assignments_from_Trinotate_xls.pl \
--Trinotate_xls trinotate_annotation_report.xls \
-G --include_ancestral_terms > trinotate_annotation_report_GO.txt

#### run_GOseq

