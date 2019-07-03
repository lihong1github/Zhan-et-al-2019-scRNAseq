# install.packages("Seurat")
library(Seurat)
## Dataset for analysis

# load data in seurat, change directory name if needed
experiment.data <- Read10X(data.dir = "~/desktop/10x_results")

## create a seurate object
experiment.aggregate <- CreateSeuratObject(
  raw.data = experiment.data,
  project = "PLX", 
  min.cells = 10,  # to keep a gene, that gene must have at least in 10 cells
  min.genes = 200, # a cell must have at least 200 genes
  names.field = 2, # names are stored in the 2nd field
  names.delim = "\\-") 

### Calc mitocondrial content
# Calculate percent mitochondrial genes per cell. 

mito.genes <- grep("^mt-", rownames(experiment.aggregate@data), value = T)
percent.mito <- Matrix::colSums(experiment.aggregate@raw.data[mito.genes, ]) / 
  Matrix::colSums(experiment.aggregate@raw.data)

# add mito percentage to metadata
experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = percent.mito,
  col.name= "percent.mito")

### Lets fix the sample names, reassign names with more meaningful factors
current.samplenames = c("1", "2", "3")
new.samplenames = c("Ctrl", "D0", "D2")
library(dplyr)

# change sample name in ident slot
experiment.aggregate@ident <- plyr::mapvalues(x = experiment.aggregate@ident, from = current.samplenames, to = new.samplenames)
table(experiment.aggregate@ident)

# change sample name in meta.data$orig.ident
experiment.aggregate@meta.data$orig.ident <- plyr::mapvalues(x = experiment.aggregate@meta.data$orig.ident, from = current.samplenames, to = new.samplenames)
table(experiment.aggregate@meta.data$orig.ident)

# save the data as a R object
save(experiment.aggregate,file="~/desktop/PLX_Rdata/original_seurat_object.RData")


