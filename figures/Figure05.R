
rm(list=ls())
set.seed(123456)

library(ggpubr)
library(openxlsx)
library(tidyverse)

# GENCODE
gencode=rtracklayer::import("/home/zhangliandong/ref/GENCODE/hg38GENCODE/gencode.v26.annotation.gtf",format = "gtf")
gencode=as.data.frame(gencode)
genelist=gencode[gencode$type=="gene",c("gene_id","gene_name")]

colpalette <- c("#68a9cf","#f58456")
# Plot CD4 memory T cell -----------------------------------------------------------
# ABIS
df.estimatedComp = read.table("/home/zhangliandong/project/meniere/RNA/deconv/ABIS/Meniere_RNASEQC_InputTPM.txt_deconvolution.txt") %>% t() %>% data.frame()
df.estimatedComp[df.estimatedComp<0] = 0

sample_info = as.data.frame(rownames(df.estimatedComp))

colnames(sample_info)[1] = "Sample"
sample_info$condition = stringr::str_extract(sample_info$Sample,"^[A-Z]") %>% as.factor()
sample_info$condition = ifelse(sample_info$condition == "M","MD","Health")

dat = df.estimatedComp %>%
  tibble::rownames_to_column("Sample") %>%
  reshape2::melt()
colnames(dat) = c("Sample","Cell_type","Proportion")

dat = dplyr::left_join(dat,sample_info,by = "Sample")

p <- ggplot(dat,aes(x = Cell_type,y = Proportion,fill = condition)) +
  geom_boxplot(width = .2,alpha = .7,show.legend = T)+
  scale_fill_manual(values = colpalette)  +
  theme_bw(base_size = 16) +
  labs(y = expression("ABIS Score"),x = NULL) +
  theme(legend.position = 'top',legend.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90))
p + stat_compare_means(label = "p.format",method = "wilcox.test")

dat = dat[dat$Cell_type=="T.CD4.Memory",]
ggplot(dat,aes(x = Cell_type,y = Proportion,fill = condition)) +
  introdataviz::geom_split_violin(alpha = .6, trim = F,color = NA,width = 1) +
  geom_boxplot(width = .2,alpha = .7,show.legend = F)+
  scale_fill_manual(values = colpalette)  +
  theme_bw(base_size = 16) +
  labs(subtitle = "ABIS Score", y = NULL,x = NULL) +
  theme(legend.position = NULL,legend.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90)) +
  stat_compare_means(label = "p.format",method = "wilcox.test")


# CIBERSORT
load("/home/zhangliandong/project/meniere/RNA/deconv/CIB/LM22.composit.rda")

p <- ggplot(df.estimatedComp,aes(x = Cell_type,y = Proportion,fill = condition)) +
  geom_boxplot(width = .2,alpha = .7,show.legend = T)+
  theme_bw(base_size = 16) +
  scale_fill_manual(values = colpalette)  +
  labs(y = expression("CIBERSORT Score"), x = NULL) +
  theme(legend.position = 'top',legend.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90))
p + stat_compare_means(label = "p.format",method = "wilcox.test")


df.estimatedComp = df.estimatedComp[df.estimatedComp$Cell_type=="CD4.memory",]

p <- ggplot(df.estimatedComp,aes(x = Cell_type,y = Proportion,fill = condition)) +
  introdataviz::geom_split_violin(alpha = .6, trim = F,color = NA,width = 1) +
  geom_boxplot(width = .2,alpha = .7,show.legend = F)+
  scale_fill_manual(values = colpalette)  +
  theme_bw(base_size = 16) +
  labs(subtitle = "CIBERSORT Score", y= NULL, x = NULL) +
  theme(legend.position = NULL,legend.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90))

p + stat_compare_means()

# Plot CD4 Th1 cell -----------------------
load("gsva_ssgsea.rda")
gsva_ssgsea = gsva_ssgsea[gsva_ssgsea$Cell_type=="Th1.related.signature",]
ggplot(gsva_ssgsea,aes(x = Cell_type,y = Proportion,fill = condition)) +
  introdataviz::geom_split_violin(alpha = .6, trim = F,color = NA,width = 1) +
  geom_boxplot(width = .2,alpha = .7,show.legend = F)+
  scale_fill_manual(values = colpalette)  +
  theme_bw(base_size = 16) +
  #scale_y_continuous(breaks = c(1.2,1.3,1.4),limits = c(1.2,1.4)) +
  labs(subtitle = "ssGSEA Score", y=NULL, x = NULL) +
  theme(legend.position = NULL,legend.title = element_text(size = 12)) +
  stat_compare_means(label = "p.format",method = "wilcox.test")


# GSEA ImmuneCell --------------------------
ImmuneCell = read.xlsx("/home/zhangliandong/scripts/Meniere/Metadata/1-s2.0-S0923181123001421-mmc2.xlsx")
colnames(ImmuneCell) = ImmuneCell[1,]
ImmuneCell = ImmuneCell[-1,c(1:4,11)]

ImmuneCell_signature = ImmuneCell %>%
  pivot_longer(cols = everything(), names_to = "celltype", values_to = "marker") %>%
  filter(!is.na(marker), nzchar(trimws(marker))) %>%
  filter(grepl("signature", celltype, ignore.case = TRUE))

de_res = fread("DESeq2_RNASeQC.txt",data.table = F) %>% na.omit()
colnames(de_res)[1] = "gene_id"
de_res = dplyr::left_join(de_res,genelist,by = "gene_id")

gene_sets_lfc = de_res %>%
  filter(!duplicated(gene_name) ) %>%
  arrange( desc(log2FoldChange) )
lfc = gene_sets_lfc$log2FoldChange
names(lfc) = gene_sets_lfc$gene_name

hallmark_gsea_res <- as.data.frame(GSEA(lfc, TERM2GENE = ImmuneCell_signature, minGSSize = 5, pvalueCutoff = 0.5))

# Plot
hallmark_gsea_res$ID <- gsub(".related.signature", "", hallmark_gsea_res$ID)

hallmark_levels =
  hallmark_gsea_res %>%
  filter(qvalue < 0.1) %>%
  dplyr::select(ID, NES) %>%
  arrange(NES)

to_plot =
  hallmark_gsea_res %>%
  mutate(padj = p.adjust) %>%
  mutate(padj_label = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.1 ~ "*"
  ))

to_plot$set =  "MD/Health"
to_plot %>%
  ggplot(aes(
    x = set,
    y = ID, label = padj_label, fill = NES)) +
  geom_tile() +
  geom_text(size = 12) +
  scale_fill_distiller(na.value = "lightgray", palette = "RdBu", limits = c(-2.5,2.5) ) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  coord_fixed(ratio = 0.7) +
  labs(x = "", y = "", fill = "NES") +
  theme(plot.subtitle = element_text(hjust = 0),
        axis.text.x = element_text(face = "bold",size = 12),
        axis.text.y = element_text(face = "bold",size = 12))

