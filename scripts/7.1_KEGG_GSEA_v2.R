#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")
#install.packages("clusterProfiler")

######################

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("KEGGprofile")
#BiocManager::install("enrichKEGG")

library(clusterProfiler)
library(dplyr)
library(tidyverse)

library(enrichKEGG)
library(KEGGprofile)
library(KEGGREST)
library(pathview)
#library(emapplot)

###############

# gseMKEGG or enrichMKEGG

#https://yulab-smu.github.io/clusterProfiler-book/chapter6.html#kegg-gene-set-enrichment-analysis

#Chapter 6 KEGG analysis
#https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#######
#3.1 Input data
#######
#d <- read.csv(your_csv_file)
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_dir/RESULTS_v2/DESeq2_genes_F-EMI_vs_F-EMII/F-EMI_vs_F-EMII.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv")

#######

indir <-"/Volumes/Mac_HD2/proj_dir/RESULTS_v2/"
#basedir <- "/Volumes/Mac_HD2/proj_dir/KEGG_output_sel-pathways_ath"
basedir <- "/Volumes/Mac_HD2/proj_dir/KEGG_output_sel-pathways_mtr"
#basedir <- "/Volumes/Mac_HD2/proj_dir/KEGG_output_sel-pathways_gmx"
setwd(basedir)


#lst1=c("F-EMI_vs_F-EMII","F-EMI_vs_F-EMIII","F-EMII_vs_F-EMIII","F-ESI_vs_F-ESII","T-EMI_vs_T-EMII","T-EMI_vs_T-EMIII","T-EMII_vs_T-EMIII","T-ESI_vs_T-ESII","F-EM-I-II-III_vs_T-EM-I-II-III","F-ES-I-II_vs_T-ES-I-II","F-PS_vs_T-PS")
lst1=c("F-PS_vs_T-PS")

for (i in lst1){
#print(i)

basedir_ID <-  paste(basedir, i, sep="/")

#create new dir 
dir.create(basedir_ID)

setwd(basedir_ID)

indir_ID <- paste(indir,"/DESeq2_genes_",i,"/",i, ".DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv", sep="")

DE.table <- read.delim(indir_ID)
#print(basedir_ID1)

#DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("GeneID_ath")]),]
#DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$GeneID_ath == ".",]


DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("GeneID_mtr")]),]
DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$GeneID_mtr == ".",]

#DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("GeneID_gmx")]),]
#DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$GeneID_gmx == ".",]


# Create a vector of the gene unuiverse
kegg_gene_list <- DE.table_dedup_ids$log2FoldChange

#### AT
#names(kegg_gene_list) <- data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$GeneID_ath), "@", fixed = TRUE)))$X1
names(kegg_gene_list) <- DE.table_dedup_ids$GeneID_mtr
#names(kegg_gene_list) <- DE.table_dedup_ids$GeneID_gmx

# omit any NA values 
kegg_gene_list <- na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

##keyType       = "ncbi-geneid")
#### AT code 1
#kegg_organism = "ath"

####:Medicago truncatula
#kegg_organism = "mtr"

#### Glycine max 
#kegg_organism = "gmx"

#kk2 <- gseKEGG(geneList     = kegg_gene_list,
#               organism     = kegg_organism,
#               nPerm        = 10000,
#               minGSSize    = 3,
#               maxGSSize    = 800,
#               pvalueCutoff = 0.05,
#               pAdjustMethod = "none", keyType       = "ncbi-geneid")

#head(kk2, n = 20L)


# Map
# Produce the native KEGG plot (PNG)
#dme <- pathview(gene.data=kegg_gene_list, pathway.id="ath00620", species = "ath")

#sellist<- c("ath00500","ath00010","ath00061","ath00561","ath01200","ath00710","ath04075")
#sellist<- c("ath01200")

sellist<- c("mtr00500","mtr00010","mtr00061","mtr00561","mtr01200","mtr00710","mtr04075")

#sellist<- c("gmx00500","gmx00010","gmx00061","gmx00561","gmx01200","gmx00710","gmx04075")

#for (ID in kk2$ID)
for (ID in sellist)
{
  print(ID)
  #dme <- pathview(gene.data=kegg_gene_list, pathway.id=ID, species = "ath")
  dme <- pathview(gene.data=kegg_gene_list, pathway.id=ID, species = "mtr")
  #dme <- pathview(gene.data=kegg_gene_list, pathway.id=ID, species = "gmx")
}

#write.table(kk2, file = paste(basedir_ID,"/", i,"_gseKEGG.tsv", sep=""), row.names=T, sep="\t")

# delete xml files
#unlink("*xml")
}



#########
