#!/bin/bash -l


#
#Transcriptome annotation: Pars reference genome annotaion (GFF) file  

#######################
########

####### 1)  Arabidopsis-thaliana.blastp

for j in Trinity.fasta.transdecoder.pep_Arabidopsis-thaliana.blastp.outfmt6 #trinity_gene_splice_modeler.fasta.transdecoder.pep_Arabidopsis-thaliana.blastp.outfmt6 trinity_gene_splice_modeler.fasta_Arabidopsis-thaliana.blastx.outfmt6
do 

echo ${j}
rm Thaliana_PIDs_anno
for i in $(cat ${j} | awk '{print $2}' | sort -u)
do
#echo $i
LANG=C grep -w -m 1 $i /home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.4_TAIR10.1_genomic.gff | cut -f9 | awk -F';' '{ split($4,a,"=");print a[2],"\t",$3,"\t",$(NF-1)}' >> Thaliana_PIDs_anno
done
awk 'NR==FNR{a[$1]=$0; next} ($2 in a){split(a[$2],b,"\t");print $0,'\t',b[2],"\t",b[3]}' Thaliana_PIDs_anno ${j} > ${j}_anno.tsv

done

####### 2) Glycine-max

for j in Trinity.fasta.transdecoder.pep_Glycine-max.blastp.outfmt6 #trinity_gene_splice_modeler.fasta.transdecoder.pep_Glycine-max.blastp.outfmt6 trinity_gene_splice_modeler.fasta_Glycine-max.blastx.outfmt6
do 
echo ${j}

rm Glycine-max_PIDs_anno
for i in $(cat ${j} | awk '{print $2}' | sort -u)
do
#echo $i
LANG=C grep -w -m 1 $i /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/GCF_000004515.6_Glycine_max_v4.0_genomic.gff | cut -f9 | awk -F';' '{ split($3,a,",");split(a[1],b,"=");split($4,c,"=");print c[2],"\t",b[2],"\t",$(NF-1)}' >> Glycine-max_PIDs_anno
done

awk 'NR==FNR{a[$1]=$0; next} ($2 in a){split(a[$2],b,"\t");print $0,'\t',b[2],"\t",b[3]}' Glycine-max_PIDs_anno ${j} > ${j}_anno.tsv

done

####### 3) Medicago-truncatula


for j in Trinity.fasta.transdecoder.pep_Medicago-truncatula.blastp.outfmt6 #trinity_gene_splice_modeler.fasta.transdecoder.pep_Medicago-truncatula.blastp.outfmt6 trinity_gene_splice_modeler.fasta_Medicago-truncatula.blastx.outfmt6
do 
echo ${j}


rm Medicago-truncatula_PIDs_anno
for i in $(cat ${j} | awk '{print $2}' | sort -u)
do
#echo $i
LANG=C grep -w -m 1 $i /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gff | cut -f9 | awk -F';' '{ split($3,a,",");split(a[1],b,"=");split($4,c,"=");print c[2],"\t",b[2],"\t",$(NF-1)}' >> Medicago-truncatula_PIDs_anno
done

awk 'NR==FNR{a[$1]=$0; next} ($2 in a){split(a[$2],b,"\t");print $0,'\t',b[2],"\t",b[3]}' Medicago-truncatula_PIDs_anno ${j} > ${j}_anno.tsv

done

####### 4) Pisum-sativum


for j in Trinity.fasta.transdecoder.pep_Pisum-sativum.blastp.outfmt6 #trinity_gene_splice_modeler.fasta.transdecoder.pep_Pisum-sativum.blastp.outfmt6 trinity_gene_splice_modeler.fasta_Pisum-sativum.blastx.outfmt6
do 
echo ${j}

####
rm Pisum-sativum_PIDs_anno

for i in $(cat ${j} | awk '{print $2}' | sort -u)
do
#echo $i
grep -w -m 1 $i /home/gala0002/proj/proj_dir/REF/Pisum-sativum/Pisum_sativum_v1a_genes.gff3 | cut -f9 | awk -F';' '{ split($1,a,"=");print a[2],"\t",$2,"\t",$3}' >> Pisum-sativum_PIDs_anno

done
####

awk 'NR==FNR{a[$1]=$0; next} ($2 in a){split(a[$2],b,"\t");print $0,'\t',b[2],"\t",b[3]}' Pisum-sativum_PIDs_anno ${j} > ${j}_anno.tsv

done

##############################################
