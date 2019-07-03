library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(ggrepel)
library(scales)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")

my_color <- hue_pal()(9)

genes_to_plot <- c("Itgam", "Csf1r", "Cd14", 
                   "Ly6g", "Ly6c2", "Sell", "Klrb1c", "Nkg7",
                   "Ccr7", "Il7r", "Cd209a", "Itgax")

for (i in 1:length(genes_to_plot)) {
gene <- genes_to_plot[i]
  g <- VlnPlot(
    experiment.aggregate, 
    features.plot=gene, point.size.use = 0.1,
    use.raw=F, x.lab.rot=0, do.return=T, size.y.use = 20
  ) + ylab("Log10[nUMI+1]") + xlab("Clusters") +
    scale_fill_manual(values=my_color, guide=F) +
    theme(plot.title = element_text(hjust = 0.5, size=24),
          axis.text = element_text(colour = "black", size=20),
          axis.title.x = element_text(colour = "black", size=20, hjust = 0.5),
          axis.title.y = element_text(colour = "black", size=20, hjust = 0.5),
          aspect.ratio = 0.75,
          axis.ticks.y = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black"))
    
  
    ggsave(paste(gene, "VlnPlot.pdf", sep="_"), plot = g, device = "pdf", path = "~/Desktop/Violin_plots",
         scale = 0.5, width = 10, height = 8, units = c("in"),
         dpi = 300, limitsize = FALSE)
  
}
  

