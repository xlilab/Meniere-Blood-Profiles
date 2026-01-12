
rm(list = ls())

library(openxlsx)
library(tidyverse)

# GENCODE
gencode=rtracklayer::import("/home/zhangliandong/ref/GENCODE/hg38GENCODE/gencode.v26.annotation.gtf",format = "gtf")
gencode=as.data.frame(gencode)
genelist=gencode[gencode$type=="gene",c("gene_id","gene_name")]

# Load data
Meniere_Control_TPM = read.table("/home/zhangliandong/project/meniere/RNA/RSEM/CM/Meniere_Control.rsem_genes_tpm.txt.gz",header = T)
colnames(Meniere_Control_TPM) = as.vector(unlist(lapply(strsplit(colnames(Meniere_Control_TPM),"[_]"),function(x) x[1])))
colnames(Meniere_Control_TPM)[1] = "gene_id"

sample_info = as.data.frame(colnames(Meniere_Control_TPM))
colnames(sample_info)[1] = "Sample_ID"
sample_info$Sample = sample_info$Sample_ID
rownames(sample_info) = sample_info$Sample_ID
sample_info = sample_info[-c(1,2),]
sample_info$condition = stringr::str_extract(sample_info$Sample_ID,"^[A-Z]")
sample_info$condition = ifelse(sample_info$condition == "M","MD","Health") %>% as.factor()

Meniere_Control_TPM = dplyr::left_join(Meniere_Control_TPM,genelist,by = "gene_id")
Meniere_Control_TPM = Meniere_Control_TPM[,-c(1,2)]
Meniere_Control_TPM = Meniere_Control_TPM[,c(length(Meniere_Control_TPM),1:(length(Meniere_Control_TPM)-1))]
data = aggregate(Meniere_Control_TPM[,-1], by=list(Meniere_Control_TPM$gene_name), FUN=mean, na.rm=TRUE)

rownames(data) = data$Group.1
Meniere_Control_MAT = data[-1] %>% as.matrix()
logTPM = log2(Meniere_Control_MAT + 1)

# AIBS --------------------------------------------
write.table(data,"Meniere_RNASEQC_InputTPM.txt",sep = '\t',row.names = F)

shiny::runGitHub("ABIS")

# CIBERSORT -------------------------------------
data = data[-1]
source("/home/zhangliandong/data/paper_data/RNAdeconvolution_2022NC/Scripts/Fun_CIBERSORT.R")
write.CIB = function(data, dir) {
  data <- cbind(rownames(data), data)
  colnames(data)[1] <- "gene_name"
  write.table(data, file = dir, sep = "\t", quote = FALSE, row.names = FALSE)}

load("/home/zhangliandong/scripts/Meniere/RNA/05deconv/CIB/CIBERSORT_LM22/LM22.rda")
write.CIB(LM22,"/home/zhangliandong/project/meniere/RNA/deconv/CIB/LM22.txt")
write.CIB(data,"/home/zhangliandong/project/meniere/RNA/deconv/CIB/Meniere_genename.txt")
#write.table(Meniere_Control_TPM,"/home/zhangliandong/project/meniere/RNA/deconv/CIB/Meniere_genename.txt")
table( rownames(LM22) %in% rownames(data))

# run CIB
res = CIBERSORT(sig_matrix = "/home/zhangliandong/project/meniere/RNA/deconv/CIB/LM22.txt",
                 mixture_file = "/home/zhangliandong/project/meniere/RNA/deconv/CIB/Meniere_genename.txt",
                 perm = 1000,
                 QN = FALSE,
                 absolute=TRUE,
                 abs_method='sig.score')

results = as.data.frame(res)
df.estimatedComp = results[,-c(23:26)]
df.estimatedComp$CD4.memory = df.estimatedComp$`T cells CD4 memory activated` + df.estimatedComp$`T cells CD4 memory resting`

df.estimatedComp = df.estimatedComp %>%
  rownames_to_column("Sample") %>%
  reshape2::melt(id.vars = "Sample", variable.name = "Cell_type", value.name = "Proportion") %>%
  dplyr::left_join(sample_info, by = "Sample")

save(results,df.estimatedComp,file = "/home/zhangliandong/project/meniere/RNA/deconv/CIB/LM22.composit.rda")

# GSVA ssGSEA ----------------------------------------
ImmuneCell = read.xlsx("/home/zhangliandong/scripts/Meniere/Metadata/1-s2.0-S0923181123001421-mmc2.xlsx")
colnames(ImmuneCell) = ImmuneCell[1,]
ImmuneCell = ImmuneCell[-1,c(1:4,11)]

ImmuneCell_list = c()
for(i in 1:ncol(ImmuneCell)){ImmuneCell_list[i] = list(ImmuneCell[ImmuneCell[,i]!="",i])}
ImmuneCell_list = lapply(ImmuneCell_list, function(x) x[!is.na(x)])
names(ImmuneCell_list) = colnames(ImmuneCell)

ImmuneCell_list_genename = unlist(ImmuneCell_list) %>% unique()
table(ImmuneCell_list_genename %in% rownames(data))

gsva_mat = gsva(expr = logTPM,
                method = "ssgsea",
                gset.idx.list = ImmuneCell_list,
                verbose = T,
                parallel.sz = parallel::detectCores())

gsva_ssgsea = gsva_mat %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("Sample") %>%
  reshape2::melt(id.vars = "Sample", variable.name = "Cell_type", value.name = "Proportion") %>%
  dplyr::left_join(sample_info, by = "Sample")

save(gsva_mat,gsva_ssgsea,file = "gsva_ssgsea.rda")

