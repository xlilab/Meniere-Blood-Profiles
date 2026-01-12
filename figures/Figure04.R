
rm(list = ls())
set.seed(123456)

library(dplyr)
library(data.table)
library(openxlsx)
library(DESeq2)
library(ggplot2)
library(ggrepel)

gencode=rtracklayer::import("/home/zhangliandong/ref/GENCODE/hg38GENCODE/gencode.v26.annotation.gtf",format = "gtf")
gencode=as.data.frame(gencode)
genelist=gencode[gencode$type=="gene",c("gene_id","gene_name","gene_type")]

# Load data -----------------------------------------------------------------------
df.Meniere_Control_MvsC = fread("DESeq2_RNASeQC.txt",data.table = F) %>% na.omit()

colnames(df.Meniere_Control_MvsC)[1] = "gene_id"
df.Meniere_Control_MvsC = dplyr::left_join(df.Meniere_Control_MvsC,genelist,by = "gene_id")

# volcano plot -----------------------------------------
volcano_data <- df.Meniere_Control_MvsC %>%
  mutate(sig = case_when(
    padj >= 0.05 ~ "non_sig",
    padj < 0.05 & abs(log2FoldChange) < log2(1.5) ~ "sig",
    padj < 0.05 & abs(log2FoldChange) >= log2(1.5) ~ "sig - strong"
  )) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  mutate(class = paste(sig, direction))
table(volcano_data$sig)

volcano_data <- volcano_data %>%
  mutate(sigdirection = case_when(
    padj >= 0.05 ~ "non_sig",
    padj < 0.05 & log2FoldChange < 0 ~ "sigdown",
    padj < 0.05 & log2FoldChange > 0 ~ "sigup"
  ))

volcano_data_sig = volcano_data[!(volcano_data$sigdirection == 'non_sig'),]
table(volcano_data_sig$sigdirection)

write.xlsx(volcano_data_sig,"DESeq2_RNASeQC_significant.xlsx",sheetName = "variable",rowNames = T,sep = "\t",append = T)

title = "PBMC mRNA-seq"
subtitle = "10 MD vs 10 Health"

p_cutoff = data.table(p=volcano_data$pvalue, padj=volcano_data$padj)
p_cutoff = p_cutoff[order(p_cutoff$p),] %>% filter(padj <= 0.05)
adj_cutoff = p_cutoff[nrow(p_cutoff),]$p

volcanoplot = ggplot(volcano_data,aes(x=log2FoldChange, y=-log10(pvalue), col=class)) +
  geom_point(alpha=1, size=3.5) +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1) +
  #geom_hline(yintercept = c(-log10(0.05),-log10(adj_cutoff)), lty = 4) +
  #geom_vline(xintercept = c(-log2(1), log2(1)), lty = 4)+
  #geom_vline(xintercept = c(-log2(1.2), log2(1.2)), lty = 4)+
  scale_colour_manual(values = c("non_sig up" = "gray",
                                 "non_sig down" = "gray",
                                 "sig up" = "#EB7F56",
                                 "sig - strong up" = "#B61927",
                                 "sig down" = "#4F8FC4",
                                 "sig - strong down" = "dodgerblue4"
  )) +
  guides(colour = FALSE) +
  theme(plot.title = element_text(colour="black", size=20, angle=0, hjust=.5, vjust=.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        axis.text.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=1, vjust=0, face="plain"),
        legend.position = NULL
        #legend.title = element_text(size=15),legend.text = element_text(size=12)
  )  +
  xlim(c(-5,5))+
  labs(#title = title,
    #subtitle = subtitle,
    y = expression(-log[10]~P~value),
    x = expression(log[2]~"(Fold change,MD/Health)")) +
  geom_text_repel( fontface = "italic",
                   data = volcano_data %>% filter(padj <= 0.05 & abs(log2FoldChange)  >= 0.55),
                   aes(label = gene_name),
                   size = 3.5,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 10)

volcanoplot

pdf(file = "RNA_volcano.pdf",useDingbats = F,width = 8,height = 5)
volcanoplot
dev.off()



# GO ORA ---------------------------------------
de_res = df.Meniere_Control_MvsC

GO_list = c("/home/zhangliandong/ref/GSEA/c5.go.bp.v7.5.1.symbols.gmt",
            "/home/zhangliandong/ref/GSEA/c5.go.cc.v7.5.1.symbols.gmt",
            "/home/zhangliandong/ref/GSEA/c5.go.mf.v7.5.1.symbols.gmt")
GO_pathway = foreach(fn = GO_list, .combine = rbind)%do%{
  GO_pathway = read.gmt(fn)
  return(GO_pathway)
}
GO_pathway = GO_pathway[GO_pathway$gene %in% de_res$gene_name,]


diffgenelist = de_res[which( (de_res$padj<0.05 & abs(de_res$log2FoldChange) > log2(1.5)) ),]

x = clusterProfiler::enricher(diffgenelist$gene_name,
                              univers = de_res$gene_name,
                              TERM2GENE = GO_pathway,
                              minGSSize = 2,
                              pvalueCutoff = 0.1,
                              pAdjustMethod = "BH")

GO_table = na.omit(as.data.frame(x))
GO_table = GO_table[order(GO_table$p.adjust),]
GO_table = GO_table[GO_table$Count >3 & GO_table$p.adjust < 0.05,]
GO_table$ONTOLOGY = as.vector(unlist(lapply(strsplit(GO_table$Description,"[_]"),function(x) x[1])))
GO_table$Description = sapply(strsplit(GO_table$Description, "_"), function(x) paste(x[-1], collapse = "_"))

GO_PLOT =
  ggplot(GO_table) +
  geom_bar(aes(Count,fct_reorder(Description,ONTOLOGY),fill = ONTOLOGY),stat="identity") +
  theme_bw() +
  theme(text = element_text(size =12),
  axis.text = element_text(
    size =12,
    colour = 'black')) +
  labs(x='Number of genes',y = "") +
  scale_fill_manual(values=c(GOBP = "#79B494",  GOCC = "#848CBD", GOMF = "#D67E56"))

pdf(file = "GO_ORA.pdf",useDingbats = F,width = 8.5,height = 7)
GO_PLOT
dev.off()

# GSEA -----------------------------------------------
gene_sets_lfc = de_res %>%
  filter(!duplicated(gene_name) ) %>%
  arrange( desc(log2FoldChange) )

lfc = gene_sets_lfc$log2FoldChange
names(lfc) = gene_sets_lfc$gene_name

# MsigDB
hallmarks = read.gmt("/home/zhangliandong/ref/GSEA/h.all.v7.5.1.symbols.gmt") %>%
  mutate(term = gsub("HALLMARK|_", " ", term)) %>%
  mutate(term = gsub("INTERFERON", "IFN", term)) %>%
  mutate(term = gsub("EPITHELIAL MESENCHYMAL", "E-M", term)) %>%
  mutate(term = gsub("UNFOLDED PROTEIN RESPONSE", "UPR", term)) %>%
  mutate(term = stringr::str_trim(term))

length(unique(hallmarks$term))
hallmark_gsea_res = as.data.frame(GSEA(lfc, TERM2GENE = hallmarks, minGSSize = 5, pvalueCutoff = 0.05))

hallmark_gsea_res$ID <- gsub("^\\s+", "", hallmark_gsea_res$ID)
sig_level <- 0.01

hallmark_levels <-
  hallmark_gsea_res %>%
  filter(p.adjust < sig_level) %>%
  dplyr::select(ID, NES) %>%
  arrange(NES)

to_plot <-
  hallmark_gsea_res %>%
  mutate(padj = p.adjust) %>%
  mutate(padj_label = dplyr::case_when(
    padj < 0.0001 ~ "***",
    padj < 0.001 ~ "**",
    padj < 0.01 ~ "*",
    TRUE ~ "."
  ))

to_plot$set =  "MD/Health"

MsigDB_GSEA_PLOT =
to_plot %>%
  filter(ID %in% hallmark_levels$ID) %>%
  mutate(ID = factor(ID, levels = hallmark_levels$ID)) %>%
  ggplot(aes(
    x = set,
    y = ID, label = padj_label, fill = NES)) +
  geom_tile() +
  geom_text(size=8) +
  scale_fill_distiller(na.value = "lightgray", palette = "RdBu", limits = c(-3,3) ) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  coord_fixed(ratio = 0.75) +
  labs(x = "", y = "", fill = "NES") +
  theme(plot.subtitle = element_text(hjust = 0),
      axis.text.x = element_text(#face = "bold",
        size = 12),
      axis.text.y = element_text(#face = "bold",
        size = 12))

MsigDB_GSEA_PLOT
pdf(file = "MsigDB_GSEA_PLOT.pdf",useDingbats = F,width = 8.5,height = 7)
MsigDB_GSEA_PLOT
dev.off()

# GO
hallmarks = GO_pathway

length(unique(hallmarks$term))
hallmark_gsea_res <- as.data.frame(GSEA(lfc, TERM2GENE = hallmarks, minGSSize = 5, pvalueCutoff = 0.05))

hallmark_gsea_res$ID = gsub("^\\s+", "", hallmark_gsea_res$ID)
sig_level = 0.01

hallmark_levels <-
  hallmark_gsea_res %>%
  filter(p.adjust < sig_level) %>%
  dplyr::select(ID, NES) %>%
  arrange(NES)

to_plot <-
  hallmark_gsea_res %>%
  mutate(padj = p.adjust) %>%
  mutate(padj_label = dplyr::case_when(
    padj < 0.0001 ~ "***",
    padj < 0.001 ~ "**",
    padj < 0.01 ~ "*",
    TRUE ~ "."
  ))

to_plot$set =  "MD/Health"

GO_GSEA_PLOT =
  to_plot %>%
  filter(ID %in% hallmark_levels$ID) %>%
  mutate(ID = factor(ID, levels = hallmark_levels$ID)) %>%
  ggplot(aes(
    x = set,
    y = ID,
    label = padj_label,
    fill = NES)) +
  geom_tile() +
  geom_text(size= 7) +
  scale_fill_distiller(na.value = "lightgray", palette = "RdBu", limits = c(-3,3) ) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  coord_fixed(ratio = 0.75) +
  labs(x = "", y = "", fill = "NES") +
  theme(plot.subtitle = element_text(hjust = 0),
        axis.text.x = element_text(#face = "bold",
          size = 12),
        axis.text.y = element_text(face = "bold",
          size = 10))

pdf(file = "GO_GSEA_PLOT.pdf",useDingbats = F,width = 8.5,height = 7)
GO_GSEA_PLOT
dev.off()


