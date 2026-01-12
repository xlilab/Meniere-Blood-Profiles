
rm(list = ls())

library(ropls)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Input Data ------------------------------------
load(file = "LipidIon.uniq_KNN_impute_refQC.Rdata")
colnames(alldata)=gsub("-",".",colnames(alldata))

# Metadata --------------------------------------
sample_info=data.frame(colnames(alldata))
colnames(sample_info) = "Sample_ID"
sample_info$group = as.vector(unlist(lapply(strsplit(sample_info$Sample_ID,"[.]"),function(x) x[1])))
sample_info$group = ifelse(sample_info$group == "health","Health", ifelse(sample_info$group == "M","MD", ifelse(stringr::str_detect(sample_info$group,"QC"),"QC","No-MD"))) %>% as.factor()

# LipidClass ------------------------------------
sample_info = sample_info[sample_info$group %in% c("Health","MD"),]
alldata = alldata[colnames(alldata) %in% sample_info$Sample_ID]

LipidClass = data.frame(alldata) %>% rownames_to_column(var = "LipidIon")
LipidClass$Class = as.vector(unlist(lapply(strsplit(LipidClass$LipidIon,"[()]"),function(x) x[1])))
LipidClass$LipidIon=NULL
LipidClass_percent <- LipidClass %>%
  group_by(Class) %>%
  summarise(across(c(colnames(LipidClass)[-length(LipidClass)]),sum)) %>%
  column_to_rownames(var = "Class")

alldata = alldata %>% t() %>% as.matrix()

# PCA ------------------------------------------
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

library(plotly)
fig <-  plot_ly(PCA_out,
                x = ~p1,
                y = ~p2,
                z = ~p3,
                color = ~group,
                colors = c("#68a9cf","#f58456"))

fig

# wilcox.test -------------------------------
df.normdata = data.frame(alldata)
df.normdata$Degree <- rep(c("Health", "MD"), c(20, 20))

type1 = "MD"
type2 = "Health"

matab_num = ncol(df.normdata) - 1
test_results = data.frame(Metabolites = colnames(df.normdata)[1:matab_num])

# Mann-Whitney U test
test_results$pvalue = apply(df.normdata[, 1:matab_num], 2,
                             function(x) unlist(wilcox.test(as.numeric(x) ~ df.normdata$Degree, data = df.normdata)[3]))
test_results$padj = p.adjust(test_results$pvalue, method = "BH")
cn = paste("FC_", type1, "/", type2, sep = "")
test_results[cn] = apply(df.normdata[,1:matab_num], 2,
                          function(x)
                            mean(as.numeric(x[which(df.normdata$Degree == type1)]))/
                            mean(as.numeric(x[which(df.normdata$Degree == type2)])))
test_results$log2FoldChange = log2(test_results[,4])

save(test_results,file = "lipidomics_refQCwilcoxtest.Rdata")

# Lipid class wilcox.test ------------------------
df.normdata = data.frame(t(LipidClass_percent))
df.normdata$Steroid = df.normdata$ChE+df.normdata$ZyE + df.normdata$StE + df.normdata$SiE
df.normdata[,c("ChE","ZyE","SiE")] = NULL
colnames(df.normdata)

df.normdata$Degree = rep(c("MD", "Health"), c(20, 20))
type1 = "MD"
type2 = "Health"

matab_num = ncol(df.normdata) - 1
test_results = data.frame(Metabolites = colnames(df.normdata)[1:matab_num])

# Mann-Whitney U test
test_results$pvalue = apply(df.normdata[, 1:matab_num], 2,
                             function(x) unlist(wilcox.test(as.numeric(x) ~ df.normdata$Degree, data = df.normdata)[3]))
test_results$padj = p.adjust(test_results$pvalue, method = "BH")
cn = paste("FC_", type1, "/", type2, sep = "")
test_results[cn] = apply(df.normdata[,1:matab_num], 2,
                          function(x)
                            mean(as.numeric(x[which(df.normdata$Degree == type1)]))/
                            mean(as.numeric(x[which(df.normdata$Degree == type2)])))
test_results$log2FoldChange <- log2(test_results[,4])
write.xlsx(test_results,"lipidomicsclass_refQCwilcoxtest.xlsx",sheetName = "variable",rowNames = T,sep = "\t",append = T)

