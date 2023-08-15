## Workflow scripts for Transcriptomics analysis of Faba bean 

>
* Note: All these workflow scripts need to be adopted to the local Linux-based working environment.
* Copy this project to your local directory and access the scripts directory

```
git clone https://github.com/gvarmaslu/RNAseq_Faba-bean

cd RNAseq_Faba-bean/scripts/

```

> 
NOTE: Before running the workflow one should install all dependent software/packages and set the necessary environment paths to the working directly. 


* Script to run the entire workflow

```
bash 0.0_Run-Workflow.sh > 0.0_Run-Workflow.sh.log 2>&1

```

#### Usage of workflow scripts and run independently:

> 1.0
* Run script for RNAseq raw data filter with rRNA and Adapter removal steps

```
bash 1.1_rRNA_AdaptRM.sh > 1.1_rRNA_AdaptRM.sh.log 2>&1

```

* Run script to merge technical replicates 

```
bash 1.2_Merge_QC.sh.sh > 1.2_Merge_QC.sh.log 2>&1

```

> 2.0

* De-novo transcriptome assembly with Trinity

```
bash 2.0_de-novo-assembly_Trinity_all.sh > 2.0_de-novo-assembly_Trinity_all.sh.log 2>&1

```

> 3.0

* Assembly annotation

```
bash 3.0_Trinity_annotation_Trinity-raw.sh > 3.0_Trinity_annotation_Trinity-raw.sh.log 2>&1

```

> 4.0 

* Differential Gene or Transcript Expression (DEG/DET) analysis

```
bash 4.0_DE_CDS.sh > 4.0_DE_CDS.sh.log 2>&1

```

> 5.0

* Annotation and parsing of Differential Gene Expression (DEG/DET) output files

```
bash 5.2_Pars-annotation-run-all.sh > 5.2_Pars-annotation-run-all.sh.log 2>&1

--
5.0_Trinity_annotation_Pars_GFF.sh
5.0_clean_up_anyinput_fasta.py
5.1_Pars-annotation-DET_DEseq2.py
5.2_Pars-annotation-run-all.sh

```
> 6.0

* Annotation stats

```
bash 6.0_Annotation-stats.sh > 6.0_Annotation-stats.sh.log 2>&1

--
cd seq-len-dist/

R CMD ExN50-plot.R
python seq-len-dist.py
python seq-len-dist_v2.py
python seq-len-dist_v3.py

```

> 7.0

* Gene Ontology terms and KEGG pathway Enrichment Analysis

```
R CMD 7.1_KEGG_GSEA.R

--
7.0_KEGG_Pathway-Enrichment_gseKEGG.R
7.1_KEGG_GSEA.R
7.1_KEGG_GSEA_v2.R

```


### Citation:

When using this workflow, please cite our publication (manuscript in progress):
- Spatio-Temporal Transcriptome and Storage Compound Profiles of Developing Faba Bean (Vicia faba) Seed Tissues. 2023
- Department of Plant Breeding, Swedish University of Agricultural Sciences (SLU), SE-234 22 Lomma, Sweden

###

