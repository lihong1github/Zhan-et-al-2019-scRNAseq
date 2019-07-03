library(dplyr)
library(Seurat)
# install.packages("BiocManager") # Needed to install all Bioconductor packages
# BiocManager::install("MAST")
library(MAST)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE_add_Mac2_neg.Rdata")
table(experiment.aggregate@ident)

cluster_id <- experiment.aggregate@ident
experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = cluster_id,
  col.name = "cluster_id")

Cluster_compare <-  as.character(c(1:9, "Mac2_neg"))

df <- as.data.frame(as.matrix(experiment.aggregate@data))



  i=10
  markers <-FindMarkers(experiment.aggregate, ident.1= Cluster_compare[i],
                        ident.2=Cluster_compare[1],
                        print.bar = TRUE, logfc.threshold = 0.5,
                        test.use = "MAST")
  
  markers$gene <- rownames(markers)
  cell_names <- rownames(experiment.aggregate@meta.data)[experiment.aggregate@meta.data$cluster_id==Cluster_compare[i]]
  mean_umi_cluster <- data.frame(gene=rownames(df), mean=rowMeans(df[, colnames(df) %in% cell_names]))
  colnames(mean_umi_cluster) <- c("gene", paste("avg(logUMI)_cluster_", Cluster_compare[i], sep=""))
  
  cell_names_1 <- rownames(experiment.aggregate@meta.data)[experiment.aggregate@meta.data$cluster_id==1]
  mean_umi_cluster_1 <- data.frame(gene=rownames(df), mean=rowMeans(df[, colnames(df) %in% cell_names_1]))
  colnames(mean_umi_cluster_1) <- c("gene", paste("avg(logUMI)_cluster_", 1, sep=""))
  
  markers <- dplyr::left_join(markers, mean_umi_cluster, by="gene")     
  markers <- dplyr::left_join(markers, mean_umi_cluster_1, by="gene")
  markers$cluster <- Cluster_compare[i]
  markers <- markers[order(-markers$avg_logFC), ]
  markers_2 <- markers[, c(6, 2, 7, 8, 1, 5, 9)]
  write.csv(markers_2, file=paste("~/desktop/", "markers_UMI_cluster-", Cluster_compare[i], ".csv", sep=""), row.names=FALSE)

# load(file="~/desktop/PLX_Rdata/mac2_dev_genes.Rdata")
# 
# 
# colnames(mac2_genes_devstages) <- c("gene", "Dev stage", "Exp Direction")
# markers_2_dev <- dplyr::left_join(markers_2, mac2_genes_devstages, by="gene")
# 
# write.csv(markers_2_dev, file=paste("~/desktop/", "markers_UMI_cluster-", Cluster_compare[i], "_Dev.csv", sep=""), row.names=FALSE)


