library(RColorBrewer)
library(ggplot2)
library(Seurat)
load(file="~/desktop/PLX_Rdata/pca_sample_corrected.RData")

use.pcs = 1:22
# pick and choose pc
# use.pc2 <- c(1:2, 4:35)

experiment.aggregate <- FindClusters(
  object = experiment.aggregate, 
  reduction.type = "pca", 
  dims.use = use.pcs, 
  resolution = seq(0.1, 1, by=0.1), 
  print.output = FALSE, 
  save.SNN = TRUE,
  force.recalc = T
)


PrintFindClustersParams(object = experiment.aggregate)

?FindClusters()


# Lets first investigate how many clusters each resolution produces 
sapply(grep("^res",colnames(experiment.aggregate@meta.data),value = TRUE),
       function(x) length(unique(experiment.aggregate@meta.data[,x])))

# set the cluster based on the resolution and 
# stored in the experiment.aggregate@ident slot
experiment.aggregate <- SetAllIdent(experiment.aggregate, id = "res.0.3")

# finnaly lets produce a table of cluster to sample assignments.
table(experiment.aggregate@ident,experiment.aggregate@meta.data$orig.ident)
table(experiment.aggregate@ident)

# tSNE dimensionality reduction plots are then used to visualise clustering results. 
# As input to the tSNE, you should use the same PCs as input to the clustering analysis.

experiment.aggregate <- RunTSNE(
  object = experiment.aggregate,
  reduction.use = "pca",
  dims.use = use.pcs,
  do.fast = TRUE)


# visualize TSNE plots
dev.off()
TSNEPlot(object = experiment.aggregate, 
         pt.size=0.5,
         group.by="res.0.3",
         do.label = TRUE,
         vector.friendly=T) +
  theme(aspect.ratio = 1)

ggsave("TSNE_plot_ori.png", plot = last_plot(), device = "png", path = "~/Desktop/TSNE_plots",
       scale = 0.6, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)

## Plot TSNE coloring by the slot ‘orig.ident’ (sample names).

TSNEPlot(object = experiment.aggregate, 
         group.by="orig.ident", colors.use=c("azure4", "Deeppink1", "royalblue1"), 
         pt.size=0.5, do.return=T) + theme(aspect.ratio = 1)

ggsave("TSNE_sample_ori.png", plot = last_plot(), device = "png", path = "~/Desktop/TSNE_plots",
       scale = 0.6, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)
# Building a tree relating the ‘average’ cell from each cluster.
experiment.aggregate <- BuildClusterTree(
  experiment.aggregate,
  pcs.use = use.pcs,
  do.reorder = F,
  reorder.numeric = F,
  do.plot=F)

experiment.aggregate@cluster.tree
PlotClusterTree(experiment.aggregate)

dev.copy(png, file="~/desktop/R_plots/Cluster_trees.png", 
         width = 2000, height = 1600, res=300 )
dev.off()

# rename clusters, put 0, 1, 7 at first
n <- dim(table(experiment.aggregate@ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
# switch 3 and 8
new.cluster.ids[c(3, 8)] <- new.cluster.ids[c(8, 3)]

experiment.aggregate@ident <- plyr::mapvalues(x = experiment.aggregate@ident, from = current.cluster.ids, to = new.cluster.ids)

experiment.aggregate@ident <- factor(experiment.aggregate@ident, levels=1:n)
table(experiment.aggregate@ident)

dev.off()
TSNEPlot(object = experiment.aggregate, 
         pt.size=0.5,
         do.label = TRUE,
         vector.friendly=T)

save(experiment.aggregate, file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")


