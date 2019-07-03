library(Seurat)

load(file="~/desktop/PLX_Rdata/original_seurat_object.RData")

table(experiment.aggregate@meta.data$orig.ident)

experiment.aggregate@ident
head(experiment.aggregate@meta.data)
# Show 5% qunatiles for number of genes per cell per sample
do.call("cbind", tapply(experiment.aggregate@meta.data$nGene,
                        experiment.aggregate@ident,quantile,probs=seq(0,1,0.05)))

# Show 5% qunatiles for number of UMI per cell per sample
do.call("cbind", tapply(experiment.aggregate@meta.data$nUMI,
                        experiment.aggregate@ident,quantile,probs=seq(0,1,0.05)))

# Show 5% qunatiles for number of mitochondrial percentage per cell per sample
round(do.call("cbind", tapply(experiment.aggregate@meta.data$percent.mito,
                              experiment.aggregate@ident,quantile,probs=seq(0,1,0.05))), digits = 3)

# Show 5% qunatiles for number of ribosomal percentage per cell per sample
round(do.call("cbind", tapply(experiment.aggregate@meta.data$percent.ribo,
                              experiment.aggregate@ident,quantile,probs=seq(0,1,0.05))), digits = 3)

# Plot the number of cells with each gene is represented by

plot(sort(Matrix::rowSums(experiment.aggregate@data>=2)) , 
     xlab="gene rank", ylab="number of cells", main="Histogram: gene distribution in cells")

dev.copy(png, file="~/desktop/R_plots/histogram_genes.png", width = 2000, height = 1600, res=300 )

dev.off()


violin_plot_item <- c("nGene", "nUMI","percent.mito")

# Violin plot of 1) number of genes, 2) number of UMI and 3) percent mitochondrial genes
sample_color = c("azure4", "Deeppink1", "royalblue1")

# removed cells that showed less than 200 genes
VlnPlot(experiment.aggregate, 
        c("nGene"),
        nCol = 1, do.return = T, adjust.use = 4, point.size.use = 0.25, cols.use = sample_color) +
  xlab(NULL) +
  geom_hline(yintercept = 200, linetype="dashed", colour = "black") +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_blank(),
        panel.border = element_rect(linetype = "solid", color="black", size=1),
        axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=30,face="bold"),
        legend.title=element_blank()) +
  scale_y_continuous(breaks=seq(0, 8000, 2000), limits=c(0, 8000), expand=c(0,0)) +
  theme(aspect.ratio = 1)

ggsave(paste("QC_violin_plot_", "nGene", ".pdf", sep=""), plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 300, limitsize = FALSE)



#### remove cells that showed UMI less 500 or above 2000
VlnPlot(experiment.aggregate, 
        c("nUMI"),
        nCol = 1, do.return = T, adjust.use = 4, point.size.use = 0.25, cols.use = sample_color) +
  xlab(NULL) +
  geom_hline(yintercept = 500, linetype="dashed", colour = "black") +
  geom_hline(yintercept = 20000, linetype="dashed", colour = "black") +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_blank(),
        panel.border = element_rect(linetype = "solid", color="black", size=1),
        axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=30,face="bold"),
        legend.title=element_blank()) +
  scale_y_continuous(breaks=seq(0, 60000, 20000), limits=c(0, 60000), expand=c(0,0)) +
  theme(aspect.ratio = 1)

ggsave(paste("QC_violin_plot_", "nUMI", ".pdf", sep=""), plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 300, limitsize = FALSE)


#### remove cells that showed mito above 10%
VlnPlot(experiment.aggregate, 
        c("percent.mito"),
        nCol = 1, do.return = T, adjust.use = 4, point.size.use = 0.25, cols.use = sample_color) +
  xlab(NULL) +
  geom_hline(yintercept = 0.1, linetype="dashed", colour = "black") +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_blank(),
        panel.border = element_rect(linetype = "solid", color="black", size=1),
        axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=30,face="bold"),
        legend.title=element_blank()) +
  scale_y_continuous(breaks=seq(0, 0.3, 0.1), limits=c(0, 0.3), expand=c(0,0)) +
  theme(aspect.ratio = 1)

ggsave(paste("QC_violin_plot_", "percent.mito", ".pdf", sep=""), plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 1, width = 12, height = 6, units = c("in"),
       dpi = 300, limitsize = FALSE)



# Gene Plot, scatter plot of gene expression across cells, (colored by sample)

sample_color = c("azure4", "Deeppink1", "royalblue1")
GenePlot(
  experiment.aggregate, "nUMI", "nGene", cex.use = 0.5, pch.use = 16, col.use = sample_color)

dev.copy(png, file="~/desktop/R_plots/nGene_nUMI_coorelation.png", width = 2000, height = 1600, res=300 )

dev.off()
# We use the information above to filter out cells. 
# Here we choose those that have percent mitochondrial genes max of 10% 
# and unique UMI counts under 20,000 or greater than 500, 
# Note that low.thresholds and high.thresholds are used to define a ‘gate’ 
# -Inf and Inf should be used if you don’t want a lower or upper threshold.

table(experiment.aggregate@meta.data$orig.ident)
# remove cells that have over 10% mitochondrial genes
experiment.aggregate <- FilterCells(
  object = experiment.aggregate,
  subset.names = c("percent.mito"),
  low.thresholds = c(-Inf),
  high.thresholds = c(0.1))

table(experiment.aggregate@meta.data$orig.ident)

# remove cells that have UMI less than 500 or UMI greater than 20000
experiment.aggregate <- FilterCells(
  object = experiment.aggregate,
  subset.names = c("nUMI"),
  low.thresholds = c(500),
  high.thresholds = c(20000))

table(experiment.aggregate@meta.data$orig.ident)

# additional filtering [optional]
FilterGenes <- 
  function (object, min.value=1, min.cells = 0, genes = NULL) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FilterGenes"))]
    object <- Seurat:::SetCalcParams(object = object, calculation = "FilterGenes", ... = parameters.to.store)
    genes.use <- rownames(object@data)
    
    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- object@data[genes.use, ]
      return(object)
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(object@data > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      object@data <- object@data[genes.use, ]
      return(object)
    } else {
      return(object)
    }
  }
# filters genes requiring a min.value (log-normalized) in at least min.cells, 
# here expression of 1 in at least 100 cells.
# experiment.aggregate <- FilterGenes(object = experiment.aggregate, min.value = 1, min.cells = 100)


#  normalizes the gene expression measurements for each cell by the total expression, 
#  multiplies this by a scale factor (10,000 by default), 
#  and then log-transforms the data.
experiment.aggregate <- NormalizeData(
  object = experiment.aggregate,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

?NormalizeData()

# put the normalized data to imputed slot
norm.data <- experiment.aggregate@data
normalized_data <- as.data.frame(as.matrix(norm.data))
experiment.aggregate@imputed <- normalized_data
table(experiment.aggregate@meta.data$orig.ident)

library(Seurat)

# Identify variable genes

?FindVariableGenes()
experiment.aggregate <- FindVariableGenes(
  object = experiment.aggregate,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  x.low.cutoff = 0.125,
  x.high.cutoff = 4,
  y.cutoff = 0.5, do.plot=T)

dev.copy(png, file="~/desktop/R_plots/Find_variable_genes.png", width = 2000, height = 1600, res=300 )
dev.off()

?FindVariableGenes()

experiment.aggregate@meta.data$orig.ident
length(experiment.aggregate@var.genes)

save(experiment.aggregate, file="~/desktop/PLX_Rdata/pre_sample_corrected.RData")
