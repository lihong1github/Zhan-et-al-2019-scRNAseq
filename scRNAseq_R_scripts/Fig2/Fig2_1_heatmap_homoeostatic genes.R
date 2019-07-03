library(Seurat)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")

MG_homeostasis <- c("Trem2", "Tmem119", "P2ry12", "P2ry13", "Mafb", "Selplg",
                    "Aif1", "Cx3cr1", "Csf1r", "Hexb")

DoHeatmap(
  object = experiment.aggregate, 
  genes.use = MG_homeostasis, 
  #use.scaled=F, # if false, it uses data slot which contains normalized data (log10(UMI+1))
  data.use= experiment.aggregate@scale.data,
  remove.key = F,
  group.label.rot=F,
  rotate.key=T,
  slim.col.label=T,
  draw.line = T,
  # col.low="gray90",
  # col.mid="grey70",
  # col.high="red",
  group.cex=0,
  cex.row=10,
  do.plot=F) 

ggsave(paste("MG_homeostasis", "heatmap.pdf", sep="_"), plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 0.6, width = 12, height = 3, units = c("in"),
       dpi = 300, limitsize = F)


