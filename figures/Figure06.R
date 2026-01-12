
rm(list=ls())

library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(openxlsx)
library(clusterProfiler)

gencode=rtracklayer::import("/home/zhangliandong/ref/GENCODE/hg38GENCODE/gencode.v26.annotation.gtf",format = "gtf")
gencode=as.data.frame(gencode)
genelist=gencode[gencode$type=="gene",c("gene_id","gene_name","gene_type")]

# GSEA -----------------------------------------------------------------------------------------
hallmarks = read.gmt("/home/zhangliandong/ref/GSEA/c5.go.mf.v7.5.1.symbols.gmt") %>%
  mutate(term = gsub("GOMF_", " ", term)) %>%
  mutate(term = stringr::str_trim(term))
hallmarks = hallmarks[grep("G_PROTEIN_COUPLED_RECEPTOR",hallmarks$term),]

df.Meniere_Control_MvsC = fread("DESeq2_RNASeQC.txt",data.table = F) %>% na.omit()
colnames(df.Meniere_Control_MvsC)[1] = "gene_id"
df.Meniere_Control_MvsC = dplyr::left_join(df.Meniere_Control_MvsC,genelist,by = "gene_id")

GPCR = df.Meniere_Control_MvsC[df.Meniere_Control_MvsC$gene_name %in% hallmarks$gene,]
table(GPCR$term)
GPCR_p = GPCR[GPCR$pvalue < 0.01,]
GPCR = dplyr::left_join(GPCR,hallmarks,by = c("gene_name" = "gene"))
length(unique(GPCR$gene_id))

sig_level = 0.01

hallmark_levels <-
  GPCR %>%
  filter(pvalue < sig_level) %>%
  dplyr::select(term, gene_name, log2FoldChange) %>%
  arrange( log2FoldChange) %>%
  pivot_wider(names_from = term, values_from = log2FoldChange) %>%
  rowwise() %>%
  mutate(gene_name = factor(gene_name, levels = gene_name))

to_plot <-
  GPCR %>%
  mutate(padj_label = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.1 ~ "*"
    #TRUE ~ "."
  ))

to_plot$term = gsub("G_PROTEIN_COUPLED_RECEPTOR","GPCR",to_plot$term)

to_plot %>%
  filter(gene_name %in% hallmark_levels$gene_name) %>%
  mutate(gene_name = factor(gene_name, levels = hallmark_levels$gene_name)) %>%
  mutate(term = factor(term) ) %>%
  ggplot(aes(x = term, y = gene_name, label = padj_label, fill = log2FoldChange)) +
  geom_tile() +
  geom_text(size = 7) +
  scale_fill_distiller(na.value = "lightgray", palette = "RdBu", limits = c(-4,1) ) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic() +
  coord_fixed(ratio = 0.25) +
  labs(x = "", y = "",
    fill = expression(log[2]~"Fold change")) +
  theme(plot.subtitle = element_text(hjust = 0),
        axis.text.x = element_text(face = "bold",size = 12),
        axis.text.y = element_text(face = "italic",size = 9)
  )


# Correlation Plot -----------------------------------------------------------------------------
ImmuneCell = read.csv("/home/zhangliandong/scripts/Meniere/Metadata/ImmunCell_GSVA.csv")
ImmuneCell_signature = ImmuneCell %>% tidyr::pivot_longer(cols = everything(), names_to = "celltype", values_to = "marker")
ImmuneCell_signature = ImmuneCell_signature[nzchar(ImmuneCell_signature$marker),]
Th1 = ImmuneCell_signature[grep("Th1.related.signature",ImmuneCell_signature$celltype),]

Meniere_Control_RNA = fread("/home/zhangliandong/project/meniere/RNA/WGCNA/Meniere.expression.bed.gz",data.table = F)
colnames(Meniere_Control_RNA) = as.vector(unlist(lapply(strsplit(colnames(Meniere_Control_RNA),"[_]"),function(x) x[1])))
Meniere_Control_RNA = dplyr::left_join(Meniere_Control_RNA,genelist,by = c("gene"="gene_id"))
Th1_gene = Meniere_Control_RNA[Meniere_Control_RNA$gene_name %in% Th1$marker,-c(1:4,length(Meniere_Control_RNA)-1,length(Meniere_Control_RNA))]
Th1_Signature = apply(Th1_gene,2,sum) %>% t() %>% data.frame()
rownames(Th1_Signature) = "Th1.related.signature"

load("/home/zhangliandong/project/meniere/metabolomics/Impute/refQC/metabolomics_refQC.Rdata")
MET = alldata

#sample info
sample_info = read.table("/home/zhangliandong/scripts/Meniere/Metadata/Sample.txt")
colnames(sample_info) = c("RNA_ID","MET_ID","sex","age")
sample_info$group = as.factor(unlist(lapply(strsplit(sample_info$MET_ID,"[-]"),function(x) x[1])))
sample_info$group = ifelse(sample_info$group == "M","MD","Health")
sample_info$MET_ID = gsub("-",".",sample_info$MET_ID)

col_mapping = sample_info$RNA_ID[match(colnames(MET),sample_info$MET_ID)]
colnames(MET) = col_mapping
MET_match=MET[,unique(colnames(Th1_Signature))]
combind_data = rbind(Th1_Signature,MET_match)
#MET_match=MET[,unique(colnames(gsva_mat))]
#combind_data = rbind(gsva_mat,MET_match)

combind_data = t(combind_data) %>% data.frame()
RNAMET_TEST = combind_data[,c("C01747","Th1.related.signature")]
RNAMET_TEST$RNA_ID = rownames(RNAMET_TEST)

RNAMET_TEST = dplyr::left_join(RNAMET_TEST,sample_info,by = "RNA_ID")
colnames(RNAMET_TEST)[1:2] = c("Psychosine","Th1")

ggplot(data = RNAMET_TEST, aes(x = Psychosine, y = Th1 )) +
  geom_point(size = 5) +
  geom_smooth(method = lm) +
  labs(x = "Psychosine", y = "Th1 Sigature Genes") +
  theme_classic() +
  ggpubr::stat_cor(method = 'spearman', aes(x = Psychosine, y = Th1),size = 7) +
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 17))
