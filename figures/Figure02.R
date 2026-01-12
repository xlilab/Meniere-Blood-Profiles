
rm(list = ls())

library(ropls)
library(openxlsx)
library(tidyverse)
library(plotly)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

colpalette <- c("#c77eb5","#BEB8DC","#095f8b","#698498","#F3AA20","#2A445E","#E4DCCF","#E0B4A4","#B2BCBB","#9EA3F3","#FF9840","#E20A4A","#A200FF","#004DFF","#FFC606",
                "#9C9C4E","#192667","#EEAD0E","#ffce7b","#2585a6","#e0861a","#FFD75F","#d71345","#6b473c","#78a355","#fdb933","#5e7c85","#D2BEC1","#58094F","#FE7A15")
# Figure01A.PCA -----------------------------
load("metabolomics_refQC_PCAout.Rdata")

fig = plot_ly(PCA_out,
              x = ~p1,
              y = ~p2,
              z = ~p3,
              color = ~group,
              colors = c("#68a9cf","#f58456"))

fig

# Figure01B.Volcano plot ----------------------
load("metabolites_refQC_wilcox.RData")

volcano_data = test_results %>% mutate(change='Not significants')
volcano_data[which( (volcano_data$padj<0.05 & volcano_data$log2FoldChange > log2(1)) ),]$change = 'Up'
volcano_data[which( (volcano_data$padj<0.05 & (-log2(1) > volcano_data$log2FoldChange) )),]$change = 'Down'
volcano_data$label = ifelse(abs(volcano_data$log2FoldChange) < 1.4,"",volcano_data$Name)
volcano_data$label = ifelse(volcano_data$change %in% "Not significants","",volcano_data$label)
volcano_data[volcano_data$Name %in% c("Psychosine","D-erythro-sphingosine-1-phosphate"),]$label = volcano_data[volcano_data$Name %in% c("Psychosine","D-erythro-sphingosine-1-phosphate"),]$Name
volcano_data[volcano_data$KEGG %in% c("C00025"),]$label = volcano_data[volcano_data$KEGG %in% c("C00025"),]$Name

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.95, size=3.5) +
  scale_color_manual(values=c("#6285b5", "grey","#da3f4c"))+
  geom_vline(xintercept=c(-log2(1), log2(1)), lty = 4,col="grey",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.4) +
  labs(y = expression(-log[10]~"FDR"),
       x = expression(log[2]~"(Fold change,MD/Health)")) +
  theme_bw()+
  theme(aspect.ratio = 1) +
  xlim(c(-5,5))+
  theme(axis.title.x = element_text(size = 16, hjust = 0.5, face = "plain"),
        axis.title.y = element_text(size = 16, face = "plain"),
        axis.text.y = element_text(size = 16, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 16, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm')) +
  geom_text_repel(data = volcano_data, aes(x = log2FoldChange,
                                           y = -log10(padj),
                                           label = label),
                  max.overlaps = 15,
                  size = 4,fontface='bold',
                  point.padding = unit(0.4, "lines"),
                  segment.color = "black",
                  show.legend = FALSE)

# Figure01D.Enrich ----------------------------
load("metabolites_refQC_wilcoxMDH_kegg.RData")

ggplot(MDH_kegg, aes(Group, Description)) +
  geom_point(aes(fill = FoldEnrich,size = -log10(pvalue)), color = "black", shape = 21) +
  scale_size(range = c(1, 3),breaks = c(0,2,4)) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  ggdendro::theme_dendro() +
  theme_bw() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.7) +
  labs(x = "", y = "", title = "", fill = "FoldEnrich", size = expression(-log[10]~"P value")) +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"))

