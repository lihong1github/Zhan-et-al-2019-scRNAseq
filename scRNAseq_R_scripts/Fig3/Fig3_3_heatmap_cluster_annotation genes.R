library(Seurat)
library(tools)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")

cluster_1 <- c("Trem2", "Tmem119", "P2ry12", "P2ry13", "Hexb")
cluster_2 <- c("Ccl3", "Ccl4", "Tnf", "Cx3cr1", "Csf1r")
cluster_3 <- c("Pf4", "Apoe", "Ccl7", "Mrc1", "Dab2")
cluster_4 <- c("Plac8", "Lyz2", "Chil3", "Ear2", "Ifi27l2a")
cluster_5 <- c("S100a8", "S100a9", "Ngp", "Retnlg", "Camp")
cluster_6 <- c("Gzma", "Nkg7", "Gzmb", "AW112010", "Xcl1")
cluster_7 <- tolower(c("CD40", "CD83", "STAT1", "IL12B", "CCR7"))
cluster_7 <- toTitleCase(cluster_7)
cluster_8 <- c("H2-Ab1", "Cd74", "H2-Eb1", "H2-Aa", "Cd209a")
cluster_9 <- c("Birc5", "Top2a", "Ccna2", "Ccnb2", "Ube2c")

all_genes <- c(cluster_1, cluster_2, cluster_3, cluster_4,
               cluster_5, cluster_6, cluster_7, cluster_8,
               cluster_9)


g <- DoHeatmap(
  object = experiment.aggregate, 
  genes.use = all_genes, 
  #use.scaled=F, # if false, it uses data slot which contains normalized data (log10(UMI+1))
  data.use= experiment.aggregate@scale.data,
  remove.key = F,
  group.label.rot=F,
  rotate.key=T,
  slim.col.label=T,
  draw.line = F,
  # col.low="gray90",
  # col.mid="grey70",
  # col.high="red",
  group.cex=0,
  cex.row=8,
  do.plot=F) 

ggsave(paste("Cluster_annotation", "heatmap.pdf", sep="_"), plot = g, device = "pdf", path = "~/Desktop/",
       scale = 0.6, width = 12, height = 10, units = c("in"),
       dpi = 300, limitsize = F)


