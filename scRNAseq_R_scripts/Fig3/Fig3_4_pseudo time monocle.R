# source("http://bioconductor.org/biocLite.R")
# biocLite("monocle")

library(monocle)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")

head(experiment.aggregate@meta.data)

# add cluster_id into meta.data
cluster_id <- experiment.aggregate@ident

experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = cluster_id,
  col.name = "cluster_id")

head(experiment.aggregate@meta.data)
table(experiment.aggregate@meta.data$cluster_id)
class(experiment.aggregate@meta.data$cluster_id)

seurat_import <-importCDS(experiment.aggregate)

mono_object <- newCellDataSet(exprs(seurat_import),
                              phenoData = new("AnnotatedDataFrame", data = pData(seurat_import)),
                              featureData = new("AnnotatedDataFrame", data = fData(seurat_import)),
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())


mono_object <- estimateSizeFactors(mono_object)
mono_object <- estimateDispersions(mono_object)

mono_object <- detectGenes(mono_object, min_expr = 0.1)
print(head(fData(mono_object)))
print(head(pData(mono_object)))

expressed_genes <- row.names(subset(fData(mono_object),
                                    num_cells_expressed >= 10))

pData(mono_object)$Total_mRNAs <- Matrix::colSums(exprs(mono_object))

mono_object <- mono_object[,pData(mono_object)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(mono_object)$Total_mRNAs)) +
                     2*sd(log10(pData(mono_object)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(mono_object)$Total_mRNAs)) -
                     2*sd(log10(pData(mono_object)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(mono_object), color=cluster_id, geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

library(reshape2)
L <- log(exprs(mono_object[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

# generate cluster without markers
disp_table <- dispersionTable(mono_object)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
mono_object <- setOrderingFilter(mono_object, unsup_clustering_genes$gene_id)
plot_ordering_genes(mono_object)

plot_pc_variance_explained(mono_object, return_all = F) # norm_method='log'

mono_object <- reduceDimension(mono_object, max_components = 2, num_dim = 22,
                        reduction_method = 'tSNE', verbose = T)

mono_object <- clusterCells(mono_object, num_clusters = 12)
mono_object@phenoData
plot_cell_clusters(mono_object, 1,2, color="cluster_id")

plot_cell_clusters(mono_object, 1,2, color="orig.ident")
### making pseudo time
diff_test_res <- differentialGeneTest(mono_object[expressed_genes,],
                                      fullModelFormulaStr = "~cluster_id")

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

mono_object <- setOrderingFilter(mono_object, ordering_genes)
plot_ordering_genes(mono_object)

mono_object <- reduceDimension(mono_object, max_components = 2,
                            method = 'DDRTree')


mono_object <- orderCells(mono_object)

g <- plot_cell_trajectory(mono_object, color_by = "Pseudotime", cell_size = 0.5) +
  theme(aspect.ratio = 2, legend.position = "bottom") 

ggsave("pseduotime_plot_legend.pdf", plot = g, device = "pdf", path = "~/Desktop/pseudotime/",
       scale = 0.8, width = 14, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)


g <- plot_cell_trajectory(mono_object, color_by = c("cluster_id"), cell_size = 0.5) + 
  theme(aspect.ratio = 2) +
  facet_wrap(~cluster_id, nrow = 1) + 
  theme(legend.position="none") 

ggsave("Clusters_pseduotime.pdf", plot = g, device = "pdf", path = "~/Desktop/pseudotime/",
       scale = 0.8, width = 14, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)

