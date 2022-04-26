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
basedir <- "/Volumes/Mac_HD2/proj_dir/KEGG_output"
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

DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("GeneID_ath")]),]
DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$GeneID_ath == ".",]

# Create a vector of the gene unuiverse
kegg_gene_list <- DE.table_dedup_ids$log2FoldChange

#### AT
names(kegg_gene_list) <- data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$GeneID_ath), "@", fixed = TRUE)))$X2

# omit any NA values 
kegg_gene_list <- na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

  
#### AT code 1
kegg_organism = "ath"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

head(kk2, n = 20L)


# save boxplot
pdf(file= paste(basedir_ID,"/", i,"_Enriched_Pathways_boxplot.pdf", sep="") )

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

dev.off()

# Map
# Produce the native KEGG plot (PNG)
#dme <- pathview(gene.data=kegg_gene_list, pathway.id="ath00620", species = "ath")

for (ID in kk2$ID)
{
  #print(id)
  dme <- pathview(gene.data=kegg_gene_list, pathway.id=ID, species = "ath")
}

write.table(kk2, file = paste(basedir_ID,"/", i,"_gseKEGG.tsv", sep=""), row.names=T, sep="\t")

# delete xml files
unlink("*xml")
}



#########################################################
#Prepare Input
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
# remove duplicates 

#### AT
DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("GeneID_ath")]),]
DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$GeneID_ath == ".",]
#### MTR
DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("GeneID_mtr")]),]
DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$GeneID_mtr == ".",]
#### GMX
DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("GeneID_gmx")]),]
DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$GeneID_gmx == ".",]

### remove "YP_IDs", # remove word with charecter 
DE.table_dedup_ids <- DE.table_dedup_ids %>% filter(!str_detect(GeneID_mtr, 'YP'))

DE.table_dedup_ids <- DE.table_dedup_ids %>% filter(!str_detect(GeneID_gmx, 'YP'))

# Create a vector of the gene unuiverse
kegg_gene_list <- DE.table_dedup_ids$log2FoldChange

#### AT
names(kegg_gene_list) <- data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$GeneID_ath), "@", fixed = TRUE)))$X1
#names(kegg_gene_list) <- data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$GeneID_ath), "@", fixed = TRUE)))$X2

#### MTR
names(kegg_gene_list) <- DE.table_dedup_ids$GeneID_mtr
#### GMX
names(kegg_gene_list) <- DE.table_dedup_ids$GeneID_gmx

# omit any NA values 
kegg_gene_list <- na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


#http://yulab-smu.top/clusterProfiler-book/chapter6.html#kegg-over-representation-test
#6.1 KEGG over-representation test
# gene <-names(kegg_gene_list)[abs(kegg_gene_list) > 2]
#gene <- na.omit(DE.table_dedup_ids$GeneID_AT)

#### AT
gene <- na.omit(data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$GeneID_ath), "@", fixed = TRUE)))$X2)
#gene <- na.omit(data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$GeneID_ath), "@", fixed = TRUE)))$X1)


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
               maxGSSize = 1000,
               pvalueCutoff = 0.05, verbose      = FALSE)

head(kk2, n = 20L)

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'ath',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none")

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
