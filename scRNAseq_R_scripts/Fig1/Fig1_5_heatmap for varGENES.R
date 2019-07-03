# source("https://bioconductor.org/biocLite.R")
# biocLite("ComplexHeatmap")
# install.packages("dendextend")
# install.packages("circlize")
library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(Seurat)
library(scales)

# make a heatmap that shows all cells and all genes in clusters
load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")

# @rawdata is raw umi
# @data is now nomralized regressed data with orig.ident and mito
# @imputed is saved from @data
# @data.scaled is regressed to var.genes

df <- experiment.aggregate@imputed

clu <- c(1:9)

df_reorder_summary <- data.frame(genes = rownames(df))
df_cluster_length <- data.frame()
for (i in 1:length(clu)) {
  cell_names <- names(experiment.aggregate@ident)[which(experiment.aggregate@ident==clu[i])]
  cluster_len <- data.frame(cluster= clu[i], n=length(cell_names))
  df_cluster_length <- rbind(df_cluster_length, cluster_len)
  df_reorder <- df[ ,colnames(df) %in% cell_names]
  df_reorder_summary <- cbind(df_reorder_summary, df_reorder)
}

df_reorder_select <- df_reorder_summary[rownames(df_reorder_summary) %in% experiment.aggregate@var.genes,]
df_reorder_select <- df_reorder_select[, -1]

my_color <- hue_pal()(9)

type <- character()
for (i in 1:length(clu)){
type <-  append(type, c(rep(as.character(i), df_cluster_length$n[i])))
}



ha = HeatmapAnnotation(df = data.frame(Clusters=type),
                       col= list(Clusters = c("1" = my_color[1],
                                            "2" = my_color[2],
                                            "3" = my_color[3],
                                            "4" = my_color[4],
                                            "5" = my_color[5],
                                            "6" = my_color[6],
                                            "7" = my_color[7],
                                            "8" = my_color[8],
                                            "9" = my_color[9]
                                        
                       )))
                        

draw(ha, 1:2100)


Heatmap(df_reorder_select, 
        col = circlize::colorRamp2(c(0, 4, 6), c("snow1", "orange", "red")),
        heatmap_legend_param = list(color_bar = "continuous"), 
        name="Normed Exp",
        top_annotation = ha, top_annotation_height = unit(3, "mm"),
        na_col = "grey", # color for NA
        color_space = "LAB",
        rect_gp = gpar(col = NA),
        row_title_gp = gpar(fontsize = 14),
        column_title = character(0),
        column_title_side = c("bottom"),
        column_title_gp = gpar(fontsize = 14),
        # column_title_rot = 0, # rotation of column title
        cluster_rows = T,
        row_order = NULL,
        #split = 4,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        #row_dend_side = c("left", "right"),
        #row_dend_width = unit(10, "mm"),
        show_row_dend = FALSE,
        row_dend_reorder = TRUE,
        #row_dend_gp = gpar(),
        # methods for clustering "ward.D", "ward.D2", "single", "complete", "average"
        
        cluster_columns = F,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2",
        #"complete" on columns using logCPM
        #"ward.D2" on columns using CPM
        column_dend_side = c("top"),
        column_dend_height = unit(20, "mm"),
        show_column_dend = F,
        #column_dend_gp = gpar(),
        column_dend_reorder = NULL,
        
        #row_order = NULL,
        #column_order = rev(column_order),
        #row_names_side = c("right", "left"),
        show_row_names = FALSE,
        #row_names_max_width = default_row_names_max_width(),
        #row_names_gp = gpar(fontsize = 12),
        column_names_side = c("bottom"),
        show_column_names = F,
        #column_names_max_height = default_column_names_max_height(),
        # column_names_gp = gpar(fontsize = 12),
        # bottom_annotation = ha1, 
        
        # top_annotation = new("HeatmapAnnotation"),
        # top_annotation_height = top_annotation@size,
        # bottom_annotation = new("HeatmapAnnotation"),
        # bottom_annotation_height = bottom_annotation@size,
        
        # # provide kmeans clustering
        # km = 4,
        # km_title = "%i",
        split = NULL,
        gap = unit(0.5, "mm"),
        #combined_name_fun = function(x) paste(x, collapse = "/"),
        width = NULL,
        #show_heatmap_legend = TRUE,
        #heatmap_legend_param = list(title = "Heatmap"),
        use_raster = TRUE
        # raster_device = c("png", "jpeg", "tiff", "CairoPNG", "CairoJPEG", "CairoTIFF"),
        # raster_quality = 2,
        # raster_device_param = list()
)


