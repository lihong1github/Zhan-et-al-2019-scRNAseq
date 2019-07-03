library(ggplot2)
library(Seurat)
library(plotly)
library(dplyr)
library(scales)
library(stringr)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")

genes_to_plot <- c("Tmem119", "P2ry12", "Mrc1", "Lgals3", "S100a8", "Nkg7", "Ccr7", "Cd74", "Mki67")

for (i in 1:length(genes_to_plot)) {
  gene = genes_to_plot
cluster_id <- experiment.aggregate@ident
gene.data <- experiment.aggregate@data[rownames(experiment.aggregate@data)==gene]
names(gene.data) <- colnames(experiment.aggregate@data)

df <- as.data.frame(cbind(cluster_id, gene.data))
df <- df[order(df$cluster_id),]
df$cluster_id <- factor(df$cluster_id, level=c(1:9)) 
# add cell numbers which will be used for x-axis
df$cells <- c(1:2100)

my_color <- hue_pal()(9)

y_label <- str_pad(gene, width = max(str_length(genes_to_plot)), side = "left")

g <- ggplot(df) +
  geom_bar(aes(x=cells, y=gene.data, fill=cluster_id), 
           stat="identity", position = "identity") + 
  scale_fill_manual(values=my_color) +
    theme(legend.position="none",
          axis.line = element_blank(),
          panel.border = element_rect(linetype = "solid", color="black", size=1),
          axis.title.y.left = element_text(size=15, angle =0, hjust=0.5, family="mono"),
          axis.title.y.right = element_blank(),
          axis.ticks.y.left = element_blank(),
          axis.text.y.left = element_blank(),
          axis.text.y.right = element_text(size=8, angle =0, hjust=0.5),
          axis.ticks.x  = element_blank(),
          axis.text.x = element_blank(),
          aspect.ratio = 0.05) +
  geom_vline(xintercept = 788, linetype="solid", colour = "black") +
  geom_vline(xintercept = 1233, linetype="solid", colour = "black") +
  geom_vline(xintercept = 1279, linetype="solid", colour = "black") +
  geom_vline(xintercept = 1450, linetype="solid", colour = "black") +
  geom_vline(xintercept = 1579, linetype="solid", colour = "black") +
  geom_vline(xintercept = 1705, linetype="solid", colour = "black") +
  geom_vline(xintercept = 1761, linetype="solid", colour = "black") +
  geom_vline(xintercept = 2058, linetype="solid", colour = "black") +
  scale_x_continuous(limits=c(0, 2100), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0, ceiling(max(df$gene.data))), breaks=c(0, ceiling(max(df$gene.data))),
                     sec.axis = dup_axis(), expand = c(0, 0)) +
 xlab(NULL) + ylab(y_label)

ggsave(paste(i, gene, "LogUMI_barplot.pdf", sep="_"), plot = g, device = "pdf", path = "~/Desktop/",
       scale = 1.2, width = 6, height = 3, units = c("in"),
       dpi = 600, limitsize = F)  

}



