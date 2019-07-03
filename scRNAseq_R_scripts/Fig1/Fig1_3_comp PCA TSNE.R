library(Seurat)
load(file="~/desktop/PLX_Rdata/pre_sample_corrected.RData")

# regress with nUMI and percent.mito, then scale all genes
experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  do.scale = T,
  do.center = T,
  vars.to.regress = c("nUMI", "percent.mito")) 
# vars.to.regress use any metadata that we used to correct for batch effects, etcs

# @rawdata is raw umi
# @imputed is saved from @data which was normalized data
# @scale.data is regressed data with orig.ident and mito and scaled


## Dimensionality reduction with PCA 
# Next we perform PCA on the scaled data. By default, 
# the genes in object@var.genes are used as input, but can be alternatively defined. 
# Running dimensionality reduction on highly variable genes can improve performance. 

experiment.aggregate <- RunPCA(
  object = experiment.aggregate,
  pc.genes = experiment.aggregate@var.genes,
  do.print = TRUE,
  pcs.print = 1:5,
  genes.print = 5,
  pcs.compute = 40,
  maxit = 500)

PrintPCAParams(experiment.aggregate)

# make PCA plot based on pc1 and pc2
sample_color = c("azure4", "Deeppink1", "royalblue1")

PCAPlot(
  object = experiment.aggregate,
  dim.1 = 1,
  dim.2 = 2, do.return=T, cols.use=sample_color ) + 
  theme(panel.background = element_rect(fill="white"),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20,face="bold")) +
  theme(aspect.ratio = 1)

ggsave("PCA_plot.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/R_plots",
       scale = 0.6, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)

# Visualize top genes associated with principal components
VizPCA(
  object = experiment.aggregate,
  pcs.use=1:2
)

dev.copy(png, file="~/desktop/R_plots/PC12_genes.png", 
         width = 2000, height = 1600, res=300 )
dev.off()

# Draws a heatmap focusing on a principal component. 
# Both cells and genes are sorted by their principal component scores. 
# Allows for nice visualization of sources of heterogeneity in the dataset.
?PCHeatmap()
dev.off()
PCHeatmap(
  object = experiment.aggregate, 
  pc.use = 1:6, 
  cells.use = 500, 
  do.balanced = TRUE, 
  label.columns = FALSE,
  use.full = FALSE,
  remove.key = FALSE)

dev.copy(png, file="~/desktop/R_plots/PC_heatmap.png", 
         width = 2000, height = 1600, res=300 )
dev.off()

# Selecting which PCs to use
# PCElbowPlot plots the standard deviations of the principle components for easy identification of an elbow in the graph. 
# This elbow often corresponds well with the significant PCs and is much faster to run.

PCElbowPlot(
  experiment.aggregate,
  num.pc = 40) + ggtitle("PC Elbow Plot") +
  theme(panel.background = element_rect(fill="lightsteelblue1"),
        axis.line = element_blank(),
        panel.grid.major = element_line(size=0.5, linetype = "solid", color="white"),
        panel.border = element_rect(linetype = "solid", color="black", size=1),
                       axis.text=element_text(size=20, colour="black"),
                       axis.title=element_text(size=20,face="bold"),
                       plot.title=element_text(size=30,face="bold"),
                       legend.title=element_blank()) +
  theme(aspect.ratio = 0.5) + 
  scale_x_continuous(expand = c(0, 0), breaks=seq(0, 40, by=10), limits=c(0, 40)) + 
  scale_y_continuous(expand = c(0, 0), breaks=seq(0, 10, by=2), limits=c(0, 10))

ggsave("PC_bow_plot.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/R_plots",
       scale = 0.6, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)

dev.off()


# The JackStraw function randomly permutes a subset of data, 
# and calculates projected PCA scores for these ‘random’ genes. 
# Then compares the PCA scores for the ‘random’ genes with the observed 
# PCA scores to determine statistical signifance. End result is a p-value 
# for each gene’s association with each principal component. 
# We identify significant PCs as those who have a strong enrichment of low p-value genes.

experiment.aggregate <- JackStraw(
  experiment.aggregate, 
  num.replicate = 100, 
  num.pc = 40,
  display.progress = TRUE
  )

# find pc number that has two lines significantly different from each other
JackStrawPlot(object = experiment.aggregate, PCs = 1:40, nCol = 5)

ggsave("JackStrawPlot.png", plot = last_plot(), device = "png", path = "~/Desktop/R_plots",
       scale = 1.2, width = 12, height = 8, units = c("in"),
       dpi = 300, limitsize = FALSE)

save(experiment.aggregate, file="~/desktop/PLX_Rdata/pca_sample_corrected.RData")
