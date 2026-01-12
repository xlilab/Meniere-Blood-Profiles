
rm(list = ls())

library(data.table)
library(tximport)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(doMC)
library(clusterProfiler)
library(org.Hs.eg.db)

gencode=rtracklayer::import("/home/zhangliandong/ref/GENCODE/hg38GENCODE/gencode.v26.annotation.gtf",format = "gtf")
gencode=as.data.frame(gencode)
genelist=gencode[gencode$type=="gene",c("gene_id","gene_name","gene_type")]

colpalette <- c("#68a9cf","#f58456")
# Load data -----------------------------
Meniere_Control = read.table("/home/zhangliandong/project/meniere/RNA/rnaseqc/Meniere.gene_reads.gct.gz",skip =2,header = T)
rownames(Meniere_Control) = Meniere_Control[,1]
Meniere_Control = Meniere_Control[,-c(1,2)]
colnames(Meniere_Control) = as.vector(unlist(lapply(strsplit(colnames(Meniere_Control),"[_]"),function(x) x[1])))

Meniere_Control$gene_id = rownames(Meniere_Control)
keep = fread("/home/zhangliandong/project/meniere/RNA/WGCNA/Meniere.expression.bed.gz",select = "gene_id",data.table = F)
Meniere_Control = Meniere_Control[Meniere_Control$gene_id %in% keep$gene_id,]

Meniere_Control$gene_id <- NULL

# Metadata -----------------------------
sample_info = read.table("/home/zhangliandong/scripts/Meniere/Metadata/RNA_sample.txt", col.names = c("Sample_ID","Sex","Age"))
row.names(sample_info) = sample_info$Sample_ID
sample_info$condition = stringr::str_extract(sample_info$Sample_ID,"^[A-Z]") %>% as.factor()
sample_info$condition = ifelse(sample_info$condition == "M","MD","Health") %>% as.factor()
sample_info$Sex = as.factor(sample_info$Sex)

idx = match(colnames(Meniere_Control),sample_info$Sample_ID)
sample_info = sample_info[idx,]

# DESeq2 ---------------------------------------
Meniere_Control_dds = DESeqDataSetFromMatrix(Meniere_Control,colData = sample_info,design = ~ Age + Sex + condition)

# PCA ------------------------------------------
vsd = vst(object = Meniere_Control_dds,blind = T)
plotPCA(vsd,"condition",ntop = 19622)

vsd = vst(object = Meniere_Control_dds,blind = T)
pca = prcomp(t(assay(vsd)))
percentVar = pca$sdev^2 / sum( pca$sdev^2 )
intgroup = c("condition","Sex")
intgroup.df = as.data.frame(colData(vsd)[, drop=FALSE])
d = data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3],intgroup.df, name=colnames(vsd))
ggplot(data=d, aes_string(x="PC1", y="PC2", color="condition",shape = "Sex")) +
  geom_point(size=3.5) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  coord_fixed()+
  theme_bw()+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4,linetype="dashed") +
  geom_hline(yintercept = 0, color = 'gray', size = 0.4,linetype="dashed") +
  scale_fill_manual(values = colpalette )+
  scale_colour_manual(values = colpalette) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1) +
  theme(axis.title.x = element_text(size = 12, hjust = 0.5, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 12, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
        legend.position = "top") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black',
                                        fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        legend.position = "top")


library(plotly)
fig <- plot_ly(d,
               x = ~PC1,
               y = ~PC2,
               z = ~PC3,
               color = ~ condition,
               colors = c("#68a9cf","#f58456")
)

fig

# DEG -----------------------------------------
Meniere_Control_dds$condition=relevel(Meniere_Control_dds$condition,ref = "Health")

Meniere_Control_deseq2 = DESeq(Meniere_Control_dds)
resultsNames(Meniere_Control_deseq2)

Meniere_Control_MvsC = results(Meniere_Control_deseq2,name = "condition_MD_vs_Health")
summary(Meniere_Control_MvsC)
df.Meniere_Control_MvsC <- Meniere_Control_MvsC %>% data.frame()
fwrite(df.Meniere_Control_MvsC,"DESeq2_RNASeQC.txt",row.names = T ,sep = "\t")
