library(Seurat)
# install.packages('ggrepel')
# devtools::install_github('VPetukhov/ggrastr')
library(Cairo)
library(ggrastr)
library(shadowtext)
library(ggrepel)
library(scales)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")

# calculating the center cooridnate for each cluster on TSNE plot
# retrive TSNE coordiantes from @dr$tsne@cell.embeddings slot


# put raw data in imputed slot
# use feature plot to highlight certain marker genes
# impute data: normalized raw UMI
df <- experiment.aggregate@dr$tsne@cell.embeddings

n <- dim(table(experiment.aggregate@ident))
clu <- 1:n

# find center coordinates of each cluster in TSNE plot
cluster_cord_summary <- data.frame()
for (i in 1:length(clu)) {
  cell_names <- names(experiment.aggregate@ident)[which(experiment.aggregate@ident==clu[i])]
  clu_cord <- df[rownames(df) %in% cell_names,]
  clu_cord_mean <- data.frame(x=mean(clu_cord[,1]), y=mean(clu_cord[,2]),
                              name=i)
  cluster_cord_summary <- rbind(cluster_cord_summary, clu_cord_mean)
}


TSNEPlot(object = experiment.aggregate, 
         #group.by=c("orig.ident"),
         pt.size=0.5,
         do.label = F,
         vector.friendly=F) + 
  geom_text(data=cluster_cord_summary, aes(x=x, y=y, label=name),
                    color="white", fontface = 'bold',size = 6) +
  theme(panel.background = element_rect(fill="black"),
        axis.line = element_blank(),
        panel.border = element_rect(linetype = "solid", color="black", size=1),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=30,face="bold"),
        legend.title=element_blank(),
        legend.direction ="vertical",
        legend.position="right") +
  theme(aspect.ratio = 1)

ggsave("TSNE_labeled.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/TSNE_plots",
       scale = 0.6, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)

ggsave("TSNE_labeled.png", plot = last_plot(), device = "png", path = "~/Desktop/TSNE_plots",
       scale = 0.6, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)
## make TSNE plot with samples

sample_color = c("azure4", "Deeppink1", "royalblue1")

TSNEPlot(object = experiment.aggregate, 
         group.by=c("orig.ident"), colors.use=sample_color, 
         pt.size=0.5,
         do.label = F,
         vector.friendly=F) +
  geom_text(data=cluster_cord_summary, aes(x=x, y=y, label=name),
            color="white", fontface = 'bold',size = 6) +
  # geom_label_repel(data=cluster_cord_summary, aes(x=x, y=y, label=name),
  #                 color="black", fill="white", fontface = 'bold',size = 5, box.padding = 0.3,
  #                 point.padding = 0.5, segment.size=0.5, label.r=0.2, segment.colour="white") +
  theme(panel.background = element_rect(fill="black"),
        axis.line = element_blank(),
        panel.border = element_rect(linetype = "solid", color="black", size=1),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=30,face="bold"),
        legend.title=element_blank(),
        legend.direction ="vertical",
        legend.position="right") +
  theme(aspect.ratio = 1)

ggsave("TSNE_by_samples.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/TSNE_plots",
       scale = 0.6, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)

ggsave("TSNE_by_samples.png", plot = last_plot(), device = "png", path = "~/Desktop/TSNE_plots",
       scale = 0.6, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)