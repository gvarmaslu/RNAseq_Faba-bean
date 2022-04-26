#!/bin/bash


#Rna-seq using STAR
#./4.0_DE_CDS_PGSC_v4.03.sh > 4.0_DE_CDS_PGSC_v4.03.sh.log 2>&1


echo "run script for rna-seq-analysis"

############################################
# 2.2. Identifying differentially expressed (DE) transcripts
############################################
#Extracting differentially expressed transcripts and generating heatmaps

Trinity="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5/Trinity"
Trinity_path="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5"
work_dir=/home/gala0002/proj/proj_dir/

cd $work_dir

############################################
#Build Transcript and Gene Expression Matrices
#https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#express-output
# perl $Trinity_path/util/align_and_estimate_abundance.pl \
# --transcripts Trinity.fasta \
# --seqType fq \
# --samples_file Metadata_CDS_Mock-PMTVWT-delta8k_merge.txt \
# --est_method salmon \
# --trinity_mode --prep_reference \
# --output_dir outdir_estimate-ab > estimate-ab.log 2>&1
	
############################################
#Run the DE analysis at the gene level

#DESeq2
#
#F-EMI_vs_F-EMII F-EMI_vs_F-EMIII 

for i in F-EMI_vs_F-EMII F-EMI_vs_F-EMIII F-EMII_vs_F-EMIII F-ESI_vs_F-ESII T-EMI_vs_T-EMII T-EMI_vs_T-EMIII T-EMII_vs_T-EMIII T-ESI_vs_T-ESII F-EM-I-II-III_vs_T-EM-I-II-III F-ES-I-II_vs_T-ES-I-II F-PS_vs_T-PS

do 
echo ${i}

mkdir -p ${work_dir}"DESeq2_genes_"${i}

cd $work_dir

$Trinity_path/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix trinity_GSM.isoform.counts_round.matrix \
--samples_file sampleinfo_all_${i}.txt \
--method DESeq2 \
--output DESeq2_genes_${i} > DGE-run_${i}.log 2>&1

#Extracting differentially expressed transcripts and generating heatmaps
#Extract those differentially expressed (DE) transcripts that are at least 4-fold (C is set to 2^(2) ) differentially expressed at a significance of <= 0.001 (-P 1e-3) in any of the pairwise sample comparisons
#-C 1.0

#cd DESeq2_genes_${i}/
#nice -n 5 $Trinity_path/Analysis/DifferentialExpression/analyze_diff_expr.pl \
#--matrix ../trinity_GSM.isoform.counts_round.matrix \
#--samples ../sampleinfo_all_${i}.txt -P 5e-2 -C 1 > DGE-analyze_${i}.log 2>&1

cd DESeq2_genes_${i}/
nice -n 5 $Trinity_path/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../trinity_GSM.isoform.TMM.EXPR.matrix \
--samples ../sampleinfo_all_${i}.txt -P 5e-2 -C 1 > DGE-analyze_${i}.log 2>&1


rename 's/trinity_GSM.isoform.TMM.EXPR.matrix.//' *

#######



#rm *_vs_P*

#rm *${i}-18_vs_*${i}-36*
#rm *${i}-36_vs_*${i}-18*

#rm diffExpr.P5e-2_C0.matrix.RData
#rm diffExpr.P5e-2_C2.matrix.RData
rm DGE-analyze_${i}.log


done


####

#find DESeq2_genes_P_vs_O -iname "*AP_vs_O*" -exec rename 's/AP_vs_O/P_vs_O/' '{}' \;
#sed -i 's/AP/P/g' *

############################################

echo "Script done...."


