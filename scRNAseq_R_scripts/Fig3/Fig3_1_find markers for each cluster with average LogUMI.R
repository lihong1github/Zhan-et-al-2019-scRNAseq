library(dplyr)
library(Seurat)
# install.packages("BiocManager") # Needed to install all Bioconductor packages
# BiocManager::install("MAST")
library(MAST)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")


# find unique genes in each cluster
table(experiment.aggregate@ident)

cluster_id <- experiment.aggregate@ident

experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = cluster_id,
  col.name = "cluster_id")


Cluster_compare <-  as.character(c(1:9))

df <- as.data.frame(as.matrix(experiment.aggregate@data))

total_mean <- data.frame(gene=rownames(df), mean=rowMeans(df))
colnames(total_mean) <- c("gene","avg(logUMI)_cluster_all")

dir <- "~/desktop/marker_each_cluster/"

for (i in 1:length(Cluster_compare)) {
  markers <-FindMarkers(experiment.aggregate, ident.1= Cluster_compare[i],
                        ident.2=Cluster_compare[-i],
                        print.bar = TRUE, logfc.threshold = 0.5,
                        test.use = "MAST")
  markers$gene <- rownames(markers)
  cell_names <- rownames(experiment.aggregate@meta.data)[experiment.aggregate@meta.data$cluster_id==i]
  mean_umi_cluster <- data.frame(gene=rownames(df), mean=rowMeans(df[, colnames(df) %in% cell_names]))
  colnames(mean_umi_cluster) <- c("gene", paste("avg(logUMI)_cluster_", i, sep=""))
  markers <- dplyr::left_join(markers, mean_umi_cluster, by="gene")     
  markers <- dplyr::left_join(markers, total_mean, by="gene")
  markers$cluster <- i
  markers <- markers[order(-markers$avg_logFC), ]
  markers_2 <- markers[, c(6, 2, 7, 8, 1, 5, 9)]
  write.csv(markers_2, file=paste(dir, "markers_UMI_cluster-", i, ".csv", sep=""), row.names=FALSE)
}




