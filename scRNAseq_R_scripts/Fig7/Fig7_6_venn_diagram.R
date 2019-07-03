library(VennDiagram)

df_all <- read.csv(file="~/desktop/markers_UMI_cluster-Mac2_pos.csv", header = T, stringsAsFactors = F)
df_ctrl <- read.csv(file="~/desktop/markers_UMI_cluster-Mac2_pos_ctrl.csv", header = T, stringsAsFactors = F)

df_all_up <- df_all$gene[df_all$avg_logFC>0]
df_all_dn <- df_all$gene[df_all$avg_logFC<0]

df_ctrl_up <- df_ctrl$gene[df_ctrl$avg_logFC>0]
df_ctrl_dn <- df_ctrl$gene[df_ctrl$avg_logFC<0]

df_up_overlap <- df_all_up[df_all_up %in% df_ctrl_up]
df_dn_overlap <- df_all_dn[df_all_dn %in% df_ctrl_dn]

## cat.pos - position of category titles, represented by degree from the
## middle of the circle

## cat.dist - distance of the category titles from the edge of the circle


dev.off()
# draw venn with upregulated genes
g <- draw.pairwise.venn(area1 = length(df_all_up), area2 = length(df_ctrl_up), cross.area = length(df_up_overlap), 
                        category = c("Mac2+ all", "Mac2+ ctrl"),fontfamily="sans", lwd=c(5,5),
                        fill = c("maroon1", "maroon1"), col = c("white", "white"), 
                        alpha = c(0.5, 0.5), cat.pos = c(0, 0), cat.dist = c(0.02, 0.02))

ggsave("venn_up_mac2.pdf", plot = g, device = "pdf", path = "~/desktop/",
       scale = 1, width = 10, height = 10, units = c("in"),
       dpi = 300, limitsize = FALSE)

#############

adjustcolor( "dodgerblue", alpha.f = 1)
dev.off()
g <- draw.pairwise.venn(area1 = length(df_all_dn), area2 = length(df_ctrl_dn), cross.area = length(df_dn_overlap), 
                        category = c("Mac2+ all", "Mac2+ ctrl"),fontfamily="sans", lwd=c(5,5),
                        fill = c("dodgerblue", "dodgerblue"), col = c("white", "white"), 
                        alpha = c(0.5, 0.5), cat.pos = c(0, 0), cat.dist = c(0.02, 0.02))


ggsave("venn_dn_mac2.pdf", plot = g, device = "pdf", path = "~/desktop/",
       scale = 1, width = 10, height = 10, units = c("in"),
       dpi = 300, limitsize = FALSE)
