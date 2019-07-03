library(ggplot2)
library(Seurat)
library(MAST)
library(dplyr)
library(ggrepel)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE_add_Mac2_pos.Rdata")

table(experiment.aggregate@ident)

Cluster_compare <- names(table(experiment.aggregate@ident))

i=10
markers <-FindMarkers(experiment.aggregate, ident.1= Cluster_compare[i],
                      ident.2=Cluster_compare[1],
                      print.bar = TRUE, logfc.threshold = 0,
                      test.use = "MAST")

markers$gene <- rownames(markers)
markers$cluster <- Cluster_compare[i]

df <- data.frame(Gene=markers$gene,
                 Log2FC=markers$avg_logFC, 
                 FDR=markers$p_val_adj)


df$colours <- c("NC")
df$colours[df$Log2FC >= 0.5 & df$FDR <= 0.05] <- c("UP")
df$colours[df$Log2FC <= -0.5 & df$FDR <= 0.05] <- c("DN")
# 
genes_select_mature <- c("Cx3cr1", "Csf1r", "Mafb", "Tmem119")
genes_to_plot_mature <- df[df$Gene %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- c("Mature")

genes_select_immature <- c("Lyz2", "Mmp8", "Mmp9", "Pf4", "Ifit3", "Cdk1", "Lgals3")
genes_to_plot_immature <- df[df$Gene %in% genes_select_immature, ]
genes_to_plot_immature$Cluster <- c("Immature")

genes_to_plot <- rbind(genes_to_plot_mature, genes_to_plot_immature)


my_color <- c("#2B8CBE", "#D7301F", "skyblue","seashell3", "plum1")
my_color_1 <- c("#2B8CBE","Grey", "#D7301F")


dev.off()
ggplot() + 
  geom_point(data=df, aes(x=Log2FC, y=-log10(FDR), colour=colours),
             shape=19, alpha=1, size=1) +
  scale_color_manual(values = my_color_1,
                     name="DEGs",
                     breaks=rev(names(table(df$colours))),
                     labels=rev(names(table(df$colours)))) +
  geom_point(data=genes_to_plot,
             aes(x=Log2FC, y=-log10(FDR)),
             shape=19, alpha=1, size=3) +
  geom_text_repel(data=genes_to_plot,
                  aes(x=Log2FC, y=-log10(FDR), label = Gene), 
                  color="black", fontface = 'bold',size = 5, box.padding = 0.5,
                  point.padding = 0.5, segment.size=0.25, segment.colour="black") +
  ylab("-Log10[FDR]") + xlab("Log2FC[Mac2+ cells\n vs Homeostatic MG1]") +
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=20),
        axis.text.y = element_text(colour = "black", size=20),
        axis.title.x = element_text(colour = "black", size=20, face="bold"),
        axis.title.y = element_text(colour = "black", size=20, face="bold")) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks=seq(-4, 4, 2), limits=c(-4, 4))

ggsave("Volcano_Mac2_pos_genes.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 0.8, width = 7, height = 7, units = c("in"),
       dpi = 600, limitsize = FALSE)
