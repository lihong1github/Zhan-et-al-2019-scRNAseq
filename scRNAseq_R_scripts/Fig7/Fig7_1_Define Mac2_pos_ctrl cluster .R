# install.packages('devtools')
# Replace '2.3.0' with your desired version
# install.packages("Hmisc")
# install.packages("ggplot2")
# devtools::install_version(package = 'Seurat', version = package_version('2.3.0'))

library(Seurat)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")

table(experiment.aggregate@meta.data$orig.ident)

cluster_id <- experiment.aggregate@ident

experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = cluster_id,
  col.name = "cluster_id")

class(experiment.aggregate@ident)

# find the lgals3 data and put it in the object
lgals3.data <- experiment.aggregate@data[rownames(experiment.aggregate@data)=="Lgals3"]
names(lgals3.data) <- colnames(experiment.aggregate@data)


experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = lgals3.data,
  col.name = "lgals3.data")
head(experiment.aggregate@meta.data)


lgals3.cutoff <- mean(lgals3.data) + 1*sd(lgals3.data)
table(experiment.aggregate@meta.data$lgals3.data>lgals3.cutoff)

cell_index <- which(experiment.aggregate@meta.data$lgals3.data>lgals3.cutoff & experiment.aggregate@meta.data$orig.ident=="Ctrl") 
cell_name <- rownames(experiment.aggregate@meta.data)[cell_index]

# add a level to ident 
levels(experiment.aggregate@ident) <- c(levels(experiment.aggregate@ident), "Mac2_pos_ctrl")
experiment.aggregate@ident[names(experiment.aggregate@ident) %in% cell_name] <- c("Mac2_pos_ctrl")
table(experiment.aggregate@ident)

save(experiment.aggregate, file="~/desktop/PLX_Rdata/Data_RunTSNE_add_Mac2_pos_ctrl.Rdata")





