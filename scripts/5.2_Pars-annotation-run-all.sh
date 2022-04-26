#!/bin/bash


#Pars DEG

echo "run script "

#./5.2_Pars-annotation-run-all.sh > 5.2_Pars-annotation-run-all.sh.log 2>&1


###############################################
#all samples with stCuSTr-D_tr
#########

workdir="/home/gala0002/proj/proj_dir/"

###
#F-ESI_vs_F-ESII


for i in F-EMI_vs_F-EMII F-EMI_vs_F-EMIII F-EMII_vs_F-EMIII F-ESI_vs_F-ESII T-EMI_vs_T-EMII T-EMI_vs_T-EMIII T-EMII_vs_T-EMIII T-ESI_vs_T-ESII F-EM-I-II-III_vs_T-EM-I-II-III F-ES-I-II_vs_T-ES-I-II F-PS_vs_T-PS

do
echo ${i}

#ii=$(echo $i | tr "-" "_")
#echo ${ii}

python 5.1_Pars-annotation-DET_DEseq2.py ${i}".DESeq2.DE_results.P5e-2_C1.DE.subset" ${i}".DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv" ${workdir}"DESeq2_genes_${i}/"

#python2.7 5.1_Pars-annotation.py ${i}"_vs_${ij}.DESeq2.DE_results" ${i}"_vs_${ij}.DESeq2.DE_results_anno.tsv" ${workdir}"DESeq2_genes_stCuSTr-D_tr_A.DH/"

done

###############################################

echo "Script done...."



###############################################
# Intersection of two time points 

#cat 155-T1_vs_155-T2.DESeq2.DE_results.P5e-2_C0.DE.subset | cut -f1  | grep "HOR" | wc -l

cat DESeq2_genes_*/*.DESeq2.DE_results.P5e-2_C0.DE.subset | cut -f1  | grep "HOR" | sort -u > P17565_4.0_Intersect/DESeq2.DE_results.P5e-2_C0.DE_Intersect.subset


cat DESeq2_genes_*/*.DESeq2.DE_results.P5e-2_C0*T1-UP.subset | cut -f1  | grep "HOR" | sort -u > P17565_4.0_Intersect/DESeq2.DE_results.P5e-2_C0.DE_Intersect-T1-UP.subset

cat DESeq2_genes_*/*.DESeq2.DE_results.P5e-2_C0*T2-UP.subset | cut -f1  | grep "HOR" | sort -u > P17565_4.0_Intersect/DESeq2.DE_results.P5e-2_C0.DE_Intersect-T2-UP.subset



########

cat DESeq2_genes_*/*.DESeq2.DE_results.P5e-2_C0.DE.subset | cut -f1  | grep "HOR" | sort | uniq -c | sort


cat DESeq2_genes_*/*.DESeq2.DE_results.P5e-2_C0*T1-UP.subset | cut -f1  | grep "HOR" | sort | uniq -c | sort 

cat DESeq2_genes_*/*.DESeq2.DE_results.P5e-2_C0*T2-UP.subset | cut -f1  | grep "HOR" | sort | uniq -c | sort 


####
cat DESeq2_genes_*/*.DESeq2.DE_results.P5e-2_C0.DE.subset | grep "HORVU6Hr1G093470.1"



cat DESeq2_genes_*/*.DESeq2.DE_results.P5e-2_C0.DE.subset | cut -f1  | grep "HOR" | grep "HORVU6Hr1G093470.1"

grep -f inter_652.txt ../DESeq2_genes_*/*.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv -w






############## Pars FC


workdir="/home/gala0002/proj/proj_Aakash/"

###
for i in 155 155x199 199 199x155 199x235 224x155 224x199 235 244

#for i in 155
do
echo ${i}

cat ${workdir}"DESeq2_genes_${i}/"${i}"-T1_vs_${i}-T2.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv" | awk -F"\t" '{if ($7 > 1 || $7 < -1 ) print $0}' > ${workdir}"DESeq2_genes_${i}/"${i}"-T1_vs_${i}-T2.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv"

cat ${workdir}"DESeq2_genes_${i}/"${i}"-T1_vs_${i}-T2.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv" | awk -F"\t" '{if ($7 > 2 || $7 < -2 ) print $0}' > ${workdir}"DESeq2_genes_${i}/"${i}"-T1_vs_${i}-T2.DESeq2.DE_results.P5e-2_C2.DE.subset_anno.tsv"


#sed -i '1i TransID	sampleA	sampleB	baseMeanA	baseMeanB	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	155-T1-R1	155-T1-R2	155-T1-R3	155-T2-R1	155-T2-R2	155-T2-R3	155x199-T1-R1	155x199-T1-R2	155x199-T1-R3	155x199-T2-R1	155x199-T2-R2	155x199-T2-R3	199-T1-R1	199-T1-R2	199-T1-R3	199-T2-R1	199-T2-R2	199-T2-R3	199x155-T1-R1	199x155-T1-R2	199x155-T1-R3	199x155-T2-R1	199x155-T2-R2	199x155-T2-R3	199x235-T1-R1	199x235-T1-R2	199x235-T1-R3	199x235-T2-R1	199x235-T2-R2	199x235-T2-R3	224x155-T1-R1	224x155-T1-R2	224x155-T1-R3	224x155-T2-R1	224x155-T2-R2	224x155-T2-R3	224x199-T1-R1	224x199-T1-R2	224x199-T1-R3	224x199-T2-R1	224x199-T2-R2	224x199-T2-R3	235-T1-R1	235-T1-R2	235-T1-R3	235-T2-R1	235-T2-R2	235-T2-R3	244-T1-R1	244-T1-R2	244-T1-R3	244-T2-R1	244-T2-R2	244-T2-R3	chromosome	gene	gene_biotype	Hordeum_vulgare.IBSC_v2_Gene-Description' ${workdir}"DESeq2_genes_${i}/"${i}"-T1_vs_${i}-T2.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv"

done


for i in 155 155x199 199 199x155 199x235 224x155 224x199 235 244

#for i in 155
do

echo ${i}

cp DESeq2_genes_${i}/*.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv DESeq2_genes_All_FC1/
cp DESeq2_genes_${i}/*.DESeq2.DE_results.P5e-2_C2.DE.subset_anno.tsv DESeq2_genes_All_FC2/

done

