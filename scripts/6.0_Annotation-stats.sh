#!/bin/bash


#Pars DEG

echo "run script "

#./6.0_Annotation-stats.sh 

###############################################
#all stats
#########
# usefull github repo
#https://github.com/yuragal/fradians-rnaseq/blob/master/2.6.annotation_basic_stats.sh

###############################################
# Total Annotated genes
###############################################

wdir=/home/gala0002/proj/proj_dir
cd $wdir
ty=$wdir/2.0_trinity_out_all_sampleinfo
tns=$wdir/4.1_Trinotate
tnn=$wdir/4.1.2_Trinotate
tnnxml=$tnn/trinotate_annotation_report.xls
bl=$wdir/4.2_blastout


##### Annotated Transcripts: without SignalP and TmHMM 
wc -l $tnn/trinotate_annotation_report_feature_map.txt 
#63809
cat $tnn/trinotate_annotation_report_feature_map.txt | sed 's/\^Tm[0-9]\+//;s/\^sigP//' | awk 'length($1)!=length($2)' | wc -l
#61359
##### Total of 61359-44251 = 17108 TD peptides that were not annotated; 5)sprot_Top_BLASTP_hit, 3)sprot_Top_BLASTX_hit, 7,8) SignalP,TmHMM
cat $tnnxml | sed 1d | awk -F'\t' '$5!="." && ($3!="." || $7!="." || $8!="."){print $5}' | sort -u | wc -l
#44351

##### BLAST annotation
cat $bl/trinity_gene_splice_modeler.fasta.transdecoder.pep_Medicago-truncatula.blastp.outfmt6 | awk '$11<=1e-10 && $3>70 {print $1}' | sort -u | wc -l 
#25606

##### Join annotaions 

##### ORFs- annotation :

wc -l trinity_gene_splice_modeler.fasta.transdecoder.pep*_anno.tsv

34724 trinity_gene_splice_modeler.fasta.transdecoder.pep_Arabidopsis-thaliana.blastp.outfmt6_anno.tsv
41639 trinity_gene_splice_modeler.fasta.transdecoder.pep_Glycine-max.blastp.outfmt6_anno.tsv
42272 trinity_gene_splice_modeler.fasta.transdecoder.pep_Medicago-truncatula.blastp.outfmt6_anno.tsv
39716 trinity_gene_splice_modeler.fasta.transdecoder.pep_Pisum-sativum.blastp.outfmt6_anno.tsv

#
cat *.transdecoder.pep_Medicago-truncatula.blastp.outfmt6_anno.tsv | cut -f1 | cut -d. -f1 | sort -u | wc -l
30005

cat *.transdecoder.pep_Glycine-max.blastp.outfmt6_anno.tsv | cut -f1 | cut -d. -f1 | sort -u | wc -l
29559

cat *.transdecoder.pep_Pisum-sativum.blastp.outfmt6_anno.tsv | cut -f1 | cut -d. -f1 | sort -u | wc -l
27966

cat *.transdecoder.pep_Arabidopsis-thaliana.blastp.outfmt6_anno.tsv | cut -f1 | cut -d. -f1 | sort -u | wc -l
24909


##### All Supertranscripts blast hits:

wc -l *_modeler.fasta_*.blastx.outfmt6_anno.tsv

51248 trinity_gene_splice_modeler.fasta_Arabidopsis-thaliana.blastx.outfmt6_anno.tsv
76030 trinity_gene_splice_modeler.fasta_Glycine-max.blastx.outfmt6_anno.tsv
80438 trinity_gene_splice_modeler.fasta_Medicago-truncatula.blastx.outfmt6_anno.tsv
72102 trinity_gene_splice_modeler.fasta_Pisum-sativum.blastx.outfmt6_anno.tsv
#

##### TF databases 
wc -l *_PlantTFDBv5.0_*.blastx.outfmt6
4269 trinity_gene_splice_modeler.fasta_PlantTFDBv5.0_Arabidopsis-thaliana.blastx.outfmt6
3737 trinity_gene_splice_modeler.fasta_PlantTFDBv5.0_Glycine-max.blastx.outfmt6
5764 trinity_gene_splice_modeler.fasta_PlantTFDBv5.0_Medicago-truncatula.blastx.outfmt6

wc -l *_iTAKv18.12_*.blastx.outfmt6
5827 trinity_gene_splice_modeler.fasta_iTAKv18.12_Arabidopsis-thaliana.blastx.outfmt6
5657 trinity_gene_splice_modeler.fasta_iTAKv18.12_Glycine-max.blastx.outfmt6
8522 trinity_gene_splice_modeler.fasta_iTAKv18.12_Medicago-truncatula.blastx.outfmt6






