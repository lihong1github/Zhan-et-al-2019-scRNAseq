library(dplyr)
library(Seurat)
# install.packages("BiocManager") # Needed to install all Bioconductor packages
# BiocManager::install("MAST")
library(MAST)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE_add_Mac2_pos_ctrl.Rdata")
table(experiment.aggregate@ident)

cluster_id <- experiment.aggregate@ident
experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = cluster_id,
  col.name = "cluster_id")

Cluster_compare <-  as.character(c(1:9, "Mac2_pos_ctrl"))

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

  # add microglia developmental gene information
  
  df <- read.csv(file="~/desktop/MG_dev.csv", stringsAsFactors = F, header = T)
  df_mac2 <- markers_2
  
  df_mac2_up <- df_mac2$gene[df_mac2$avg_logFC>0]
  df_mac2_dn <- df_mac2$gene[df_mac2$avg_logFC<0]
  
  YS_genes <- df$Gene[df$Cluster==1]
  E1_genes <- df$Gene[df$Cluster==2]
  E2_genes <- df$Gene[df$Cluster==3]
  P1_genes <- df$Gene[df$Cluster==4]
  P2_genes <- df$Gene[df$Cluster==5]
  A1_genes <- df$Gene[df$Cluster==6]
  A2_genes <- df$Gene[df$Cluster==7]
  
  ##################
  # combine all Mac2_pos Dev genes in one table
  current.clusternames <- as.character(c(1:7))
  new.clusternames <- c("YS", "E1", "E2", "P1", "P2", "A1", "A2")
  
  mac2up_genes_devstages <- df[df$Gene %in% df_mac2_up,c(1,2)]
  mac2up_genes_devstages$Cluster <- plyr::mapvalues(x = mac2up_genes_devstages$Cluster, from = current.clusternames, to = new.clusternames)
  mac2up_genes_devstages$Exp.direction <- rep("UP", length(mac2up_genes_devstages$Gene))
  
  mac2dn_genes_devstages <- df[df$Gene %in% df_mac2_dn,c(1,2)]
  mac2dn_genes_devstages$Cluster <- plyr::mapvalues(x = mac2dn_genes_devstages$Cluster, from = current.clusternames, to = new.clusternames)
  mac2dn_genes_devstages$Exp.direction <- rep("DN", length(mac2dn_genes_devstages$Gene))
  
  mac2_genes_devstages <- rbind(mac2up_genes_devstages, mac2dn_genes_devstages)
colnames(mac2_genes_devstages) <- c("gene", "Dev stage", "Exp Direction")
markers_2_dev <- dplyr::left_join(markers_2, mac2_genes_devstages, by="gene")

write.csv(markers_2_dev, file=paste("~/desktop/", "markers_UMI_cluster-", Cluster_compare[i], "_Dev.csv", sep=""), row.names=FALSE)


