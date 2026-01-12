
rm(list = ls())

library(ropls)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggrepel)

kfold = function(x){
  if ( nrow(x) >= 7){
    nn = 7
  }else{
    nn <- nrow(x)
  }
  return(nn)
}

# Load Data -----------------------------------------
load("metabolomics_refQC.Rdata")

Metabolites_info[Metabolites_info$KEGG %in% c("C00956","C00188"),]$SuperClass = "Organic acids and derivatives"
Metabolites_info[Metabolites_info$KEGG %in% c("C00423","C05629"),]$SuperClass = "Phenylpropanoids and polyketides"

# Sample select -------------------------------------
colnames(alldata)=gsub("-",".",colnames(alldata))
alldata=alldata[,c(paste("health",c(1:20),sep = "."),paste("M",c(1:20),sep = "."))]
alldata$health.5=NULL

alldata <- t(alldata)
# Metadata -----------------------------------------
sample_info = data.frame(rownames(alldata))
colnames(sample_info) = "Sample_ID"
sample_info$group = as.vector(unlist(lapply(strsplit(sample_info$Sample_ID,"[.]"),function(x) x[1])))
sample_info$group = ifelse(sample_info$group == "health","Health",ifelse(sample_info$group == "M","MD",ifelse(sample_info$group == "QC",sample_info$group,"No-MD"))) %>% as.factor()

# PCA ----------------------------------------------
print("PCA is analysising")
scaleC = "standard"

alldata.PCA = opls(x = alldata, scaleC = scaleC)
modelDF_PCA = alldata.PCA@modelDF

scoreMN_PCA = as.data.frame(alldata.PCA@scoreMN)
PCA_out = data.frame(ID = rownames(scoreMN_PCA),scoreMN_PCA,check.names=F)

print(paste0("PC1(",round(modelDF_PCA[1,1]*100,2),"%)"))
print(paste0("PC2(",round(modelDF_PCA[2,1]*100,2),"%)"))
print(paste0("PC3(",round(modelDF_PCA[3,1]*100,2),"%)"))

PCA_out = dplyr::left_join(PCA_out,sample_info,by = c("ID"="Sample_ID"))
save(PCA_out, file = "metabolomics_refQC_PCAout.Rdata")

# OPLS-DA -----------------------------------------
alldata.OPLSDA <- opls(x = alldata,
                       y = sample_info$group,
                       predI = 1,
                       orthoI = NA,
                       crossvalI = kfold(alldata),
                       scaleC = scaleC,
                       permI = 200)

save(alldata.OPLSDA, file = "metabolomics_refQC_OPLSDA.Rdata")

# wilcox.test --------------------------------------
df.normdata = data.frame(alldata)
df.normdata$Degree <- rep(c("Health", "MD"), c(19, 20))

type1 <- "MD"
type2 <- "Health"

matab_num <- ncol(df.normdata) - 1
test_results <- data.frame(KEGG = colnames(df.normdata)[1:matab_num])

# Mann-Whitney U test
test_results$pvalue <- apply(df.normdata[, 1:matab_num], 2,
                             function(x) unlist(wilcox.test(as.numeric(x) ~ df.normdata$Degree, data = df.normdata)[3]))
test_results$padj <- p.adjust(test_results$pvalue, method = "BH")

cn <- paste("FC_", type1, "/", type2, sep = "")
test_results[cn] <- apply(df.normdata[,1:matab_num], 2,
                          function(x)
                            mean(as.numeric(x[which(df.normdata$Degree == type1)]))/
                            mean(as.numeric(x[which(df.normdata$Degree == type2)])))

test_results$log2FoldChange <- log2(test_results[,4])
test_results = dplyr::left_join(test_results,Metabolites_info,by = "KEGG")
save(test_results,file = "metabolites_refQC_wilcox.RData")

# KEGG database ------------------------------------
kegg_table_cid = read.xlsx("metabolomic-ORA.kegg_table_cid.xlsx",sheet = 1,rowNames = F,colNames = T,check.names = F)
test_results_change = dplyr::left_join(kegg_table_cid,test_results,by = c("cid"="KEGG")) %>% na.omit()
path_metab = kegg_table_cid[,c(2,3)]
path_intro = kegg_table_cid[,c(2,4)] %>% unique()

# ORA ---------------------------------------------
MET = unique(test_results$KEGG)
MET_sig = test_results[test_results$padj<0.05 & abs(test_results$log2FoldChange) > log2(1.5),]$KEGG
x <- clusterProfiler::enricher(MET_sig,
                               univers = MET,
                               TERM2GENE = path_metab,
                               TERM2NAME = path_intro,
                               minGSSize = 2,
                               pvalueCutoff = 0.1,
                               pAdjustMethod = "BH")
kegg_table <- na.omit(as.data.frame(x))
kegg_table <- kegg_table[order(kegg_table$p.adjust),]

kegg_table$Description = stringr::str_remove(kegg_table$Description," - Homo sapiens \\(human\\)")

path_avelog2FC <- function(x){
  metabs <- unlist(strsplit(x[8], "[/]"))
  mean(abs(test_results[match(metabs, test_results$KEGG),5]))
}
kegg_table$metabs_mean <- apply(kegg_table, 1, path_avelog2FC)

kegg_table$FoldEnrich <- apply(kegg_table, 1,
                               function(x) as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                                 as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                                 as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                                 as.numeric(unlist(strsplit(x[4], '[/]'))[1]))
kegg_table$Group <- "MDvsHC"

MDH_kegg <- kegg_table
save(MDH_kegg,file = "metabolites_refQC_wilcoxMDH_kegg.RData")


