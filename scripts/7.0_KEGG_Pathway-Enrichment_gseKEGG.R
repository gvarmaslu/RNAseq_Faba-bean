#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")
#install.packages("clusterProfiler")

######################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("KEGGprofile")
BiocManager::install("enrichKEGG")

library(clusterProfiler)
library(dplyr)
library(tidyverse)

library(enrichKEGG)
library(KEGGprofile)
library(KEGGREST)

###############

# gseMKEGG or enrichMKEGG

#https://yulab-smu.github.io/clusterProfiler-book/chapter6.html#kegg-gene-set-enrichment-analysis

#Chapter 6 KEGG analysis
#https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#######
#3.1 Input data
#######
#d <- read.csv(your_csv_file)
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/Reg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv")
DE.table <- read.delim("/Volumes/Mac_HD2/proj_dir/RESULTS_v2/DESeq2_genes_F-EMI_vs_F-EMII/F-EMI_vs_F-EMII.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv")

#########################################################
#Prepare Input
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
# remove duplicates 

#DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("AT_hit")]),]
#DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$AT_hit == "No_AT_hit",]

DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("GeneID_AT")]),]
DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$GeneID_AT == ".",]

# Create a vector of the gene unuiverse
kegg_gene_list <- DE.table_dedup_ids$log2FoldChange

# Name vector with ENTREZ ids
#dedup_ids_ext_GID <- data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$Gene_ID), ":", fixed = TRUE)))
#names(kegg_gene_list) <- dedup_ids_ext_GID$X2
#names(kegg_gene_list) <- DE.table_dedup_ids$GeneID_AT

names(kegg_gene_list) <- data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$GeneID_AT), "@", fixed = TRUE)))$X1
names(kegg_gene_list) <- data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$GeneID_AT), "@", fixed = TRUE)))$X2

# omit any NA values 
kegg_gene_list <- na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

#http://yulab-smu.top/clusterProfiler-book/chapter6.html#kegg-over-representation-test
#6.1 KEGG over-representation test
# gene <-names(kegg_gene_list)[abs(kegg_gene_list) > 2]

gene <- na.omit(data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$GeneID_AT), "@", fixed = TRUE)))$X2)

#gene <- na.omit(DE.table_dedup_ids$GeneID_AT)

# quinoa: cqi 
# Arabidopsis thaliana: ath

# 6.1 KEGG over-representation test
kk <- enrichKEGG(gene         = gene,
                 organism     = 'ath',
                 pvalueCutoff = 0.05)
head(kk)


#Create gseKEGG object
#6.2 KEGG Gene Set Enrichment Analysis

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'ath',
               minGSSize    = 5,
               maxGSSize = 100,
               pvalueCutoff = 0.05, verbose      = FALSE)

head(kk2, n = 20L)

#6.3 KEGG Module over-representation test

mkk <- enrichMKEGG(gene = gene,
                   organism = 'ath')

# 6.4 KEGG Module Gene Set Enrichment Analysis

mkk2 <- gseMKEGG(geneList = kegg_gene_list,
                 organism = 'ath')

#Merge DF

outdat_PM <- kk2
outdat_delta <- kk2

joined_df <- merge(outdat_PM, outdat_delta, by.x = "ID", 
                   by.y = "ID", all.x = T, all.y = T)

###plot

require(DOSE)
#dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
#dotplot(outdat_PM , showCategory = 24, title = "Enriched Pathways in PMTVWT" , split=".sign") + facet_grid(.~.sign)
#dotplot(outdat_delta, showCategory = 24, title = "Enriched Pathways in Delta8K" , split=".sign") + facet_grid(.~.sign)



#dataframe3 <- merge(outdat_PM, outdat_delta, by=c("pathway.code"),all=TRUE)

#write.table(kk, file = "/Volumes/Mac_HD2/proj_dir/RESULTS_2021/Reg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars_gseKEGG_join.tsv", row.names=T, sep="\t")

write.table(kk, file = "/Volumes/Mac_HD2/proj_dir/RESULTS_2021/Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars_gseKEGG_join.tsv", row.names=T, sep="\t")

#########################################################

DE.table <- read.delim("/Volumes/Mac_HD2/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv")

dedup_ids = DE.table[!duplicated(DE.table[c("Gene_ID")]),]

# Create a vector of the gene unuiverse
kegg_gene_list <- dedup_ids$log2FoldChange

# 
dedup_ids_ext_GID <- data.frame(do.call("rbind", strsplit(as.character(dedup_ids$Gene_ID), ":", fixed = TRUE)))

names(kegg_gene_list) <- dedup_ids_ext_GID$X2


## assume that 1st column is ID
## 2nd column is fold change
## feature 1: numeric vector
#geneList <- DE.table[,2]

## feature 2: named vector
#names(geneList) <- as.character(DE.table[,1])

## feature 3: decreasing order
geneList <- sort(kegg_gene_list, decreasing = TRUE)

#######
