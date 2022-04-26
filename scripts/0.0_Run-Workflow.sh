#!/bin/bash -l


# Run the entire workflow

#bash 0.0_Run-Workflow.sh > 0.0_Run-Workflow.sh.log 2>&1

# Processing the raw FASTQ data
./1.0_Merge_QC.sh
./1.1_rRNA_AdaptRM.sh
./1.2_Merge_QC.sh

# De-novo Assembly 
./2.0_de-novo-assembly_Trinity_all.sh
./2.0_de-novo-assembly_Trinity_all_sampleinfo.sh

# Assembly Annotation 
./3.0_Trinity_annotation_Trinity-raw.sh
./3.1_Trinity_annotation_trinity_GSM.sh
./3.2_Trinity_annotation_trinity_GSM_TF-db.sh
./3.3_Trinity_annotation_Trinotate.sh

# Differential Gene Expression (DEG) analysis 
./4.0_DE_CDS.sh

# Parse annotations and Statistics 
./5.0_Trinity_annotation_Pars_GFF.sh
./5.0_clean_up_anyinput_fasta.py
./5.1_Pars-annotation-DET_DEseq2.py
./5.2_Pars-annotation-run-all.sh
./6.0_Annotation-stats.sh

# GO and KEGG Gene Set Enrichment analysis (GSEA) 
./7.0_KEGG_Pathway-Enrichment_gseKEGG.R
./7.1_KEGG_GSEA.R






