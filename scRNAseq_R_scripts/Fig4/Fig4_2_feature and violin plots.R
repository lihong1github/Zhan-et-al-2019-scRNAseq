library(RColorBrewer)
library(ggplot2)
library(Seurat)
# install.packages('ggrepel')
# devtools::install_github('VPetukhov/ggrastr')
library(Cairo)
library(ggrastr)
library(shadowtext)
library(ggrepel)
library(scales)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")
# makes nice plots for genes

dev.off()

# calculating the center cooridnate for each cluster on TSNE plot
# retrive TSNE coordiantes from @dr$tsne@cell.embeddings slot
df <- experiment.aggregate@dr$tsne@cell.embeddings

clu <- c(1:9)

# find center coordinates of each cluster in TSNE plot
cluster_cord_summary <- data.frame()
for (i in 1:length(clu)) {
cell_names <- names(experiment.aggregate@ident)[which(experiment.aggregate@ident==clu[i])]
clu_cord <- df[rownames(df) %in% cell_names,]
clu_cord_mean <- data.frame(x=mean(clu_cord[,1]), y=mean(clu_cord[,2]), name=i)
cluster_cord_summary <- rbind(cluster_cord_summary, clu_cord_mean)
}


# put raw data in imputed slot
# use feature plot to highlight certain marker genes
# impute data: normalized raw UMI
# my_color <- hue_pal()(11)
dev.off()

gene <- c("Lgals3")


g <- FeaturePlot(
  experiment.aggregate, 
  features.plot=gene, # desired marker gene to plot
  pt.size=1,
  vector.friendly=F,
  no.legend = F, reduction.use = "tsne", use.imputed=T, do.return=T) 
# n = length(experiment.aggregate@imputed)
lapply(g, function(x){x + 
    scale_colour_gradientn(name ="Log10[nUMI+1]", colours = rev(brewer.pal(11,"RdYlBu"))) +
    geom_label_repel(data=cluster_cord_summary, aes(x=x, y=y, label=name), 
                    color="grey30", fill="white", fontface = 'bold',size = 5, box.padding = 0.3,
                    point.padding = 0.5, segment.size=0.5, label.r=0.2, segment.colour="gold") +
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_blank(),
          panel.border = element_rect(linetype = "solid", color="black", size=1),
          axis.text=element_text(size=20, colour="black"),
          axis.title=element_text(size=20,face="bold"),
          plot.title=element_text(size=30,face="bold"),
          legend.title=element_text(size=15,vjust=1, angle=90),
          legend.direction ="vertical",
          legend.position="right") +
    theme(aspect.ratio = 1)})

ggsave(paste(gene, "TSNE.pdf", sep="_"), plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 0.6, width = 10, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)

# make a violin plot for the gene of interest
g <- VlnPlot(
  experiment.aggregate, 
  features.plot=gene, point.size.use = 0.1,
  use.raw=F, x.lab.rot=0, do.return=T, size.y.use = 20
) + ylab("Log10[nUMI+1]") + xlab("Clusters") +
  scale_fill_manual(values=hue_pal()(9), guide=F) +
  theme(plot.title = element_text(hjust = 0.5, size=24),
        axis.text = element_text(colour = "black", size=20),
        axis.title.x = element_text(colour = "black", size=20, hjust = 0.5),
        axis.title.y = element_text(colour = "black", size=20, hjust = 0.5),
        aspect.ratio = 0.75,
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"))

ggsave(paste(gene, "VlnPlot.pdf", sep="_"), plot = g, device = "pdf", path = "~/Desktop/",
       scale = 0.5, width = 10, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)


