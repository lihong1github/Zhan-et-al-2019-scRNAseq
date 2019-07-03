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



load(file="~/desktop/PLX_Rdata/Data_RunTSNE_add_Mac2_pos_ctrl.Rdata")

table(experiment.aggregate@ident)

dev.off()

# calculating the center cooridnate for each cluster on TSNE plot
# retrive TSNE coordiantes from @dr$tsne@cell.embeddings slot
df <- experiment.aggregate@dr$tsne@cell.embeddings

n <- dim(table(experiment.aggregate@ident))[[1]]
clu <- c(1:n)

# find center coordinates of each cluster in TSNE plot
cluster_cord_summary <- data.frame()
for (i in 1:length(clu)) {
  cell_names <- names(experiment.aggregate@ident)[which(experiment.aggregate@ident==clu[i])]
  clu_cord <- df[rownames(df) %in% cell_names,]
  clu_cord_mean <- data.frame(x=mean(clu_cord[,1]), y=mean(clu_cord[,2]), name=i)
  cluster_cord_summary <- rbind(cluster_cord_summary, clu_cord_mean)
}

  gene = "Lgals3"
  RidgePlot(
    experiment.aggregate,
    use.scaled=F, cols.use = c(rep("grey88", 9), "orange"),
    features.plot=gene, do.return=T) + 
    ylab("Clusters") + xlab("Log10[UMI+1]") + 
    geom_vline(xintercept =  2.172241, linetype="dashed", colour = "grey36") +
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.text = element_text(colour = "black", size=16),
          axis.text.y = element_text(colour = "black", size=14, vjust = 0.5),
          axis.title.x = element_text(colour = "black", size=16, hjust = 0.5),
          axis.title.y = element_text(colour = "black", size=16, hjust=0.5, vjust=0.5),
          aspect.ratio = 1,
          axis.ticks.y = element_line(colour = "white"),
          axis.ticks.x = element_line(colour = "black"))
  
  ggsave(paste(gene, "RidgePlot.pdf", sep="_"), plot = last_plot(), device = "pdf", path = "~/Desktop/",
         scale = 0.4, width = 10, height = 8, units = c("in"),
         dpi = 300, limitsize = FALSE)


TSNEPlot(object = experiment.aggregate, 
         pt.size=0.5,
         do.label = F, colors.use =  c(rep("grey", 9), "deeppink"),
         vector.friendly=F) + 
  geom_text(data=cluster_cord_summary, aes(x=x, y=y, label=name),
            color="black", fontface = 'bold',size = 6) +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_blank(),
        panel.border = element_rect(linetype = "solid", color="black", size=1),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=30,face="bold"),
        legend.title=element_blank(),
        legend.direction ="vertical",
        legend.position="right") +
  theme(aspect.ratio = 1)

ggsave("TSNE_labeled_Mac2_pos_ctrl.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 0.6, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)



