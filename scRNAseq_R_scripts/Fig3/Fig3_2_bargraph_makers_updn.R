library(ggplot2)
library(reshape2)
library(dplyr)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")

Cluster_compare <-  as.character(c(1:9))

markers_all <- data.frame()
for (i in 1:length(Cluster_compare)) {
  markers <-FindMarkers(experiment.aggregate, ident.1= Cluster_compare[i],
                        ident.2=Cluster_compare[-i],
                        print.bar = TRUE, logfc.threshold = 0.5,
                        test.use = "MAST")
  markers$gene <- rownames(markers)
  markers$cluster <- i
  markers_all <- dplyr::bind_rows(markers_all, markers)
}


# filter genes with FDR<0.05 and order them by logFC
markers_all <-  markers_all %>% dplyr::filter(p_val_adj<0.05) %>% dplyr::arrange(desc(avg_logFC))
markers_all <- markers_all[order(markers_all$cluster, -markers_all$avg_logFC), ]


# count number of upregulated genes in each cluster
up_genes <- subset(markers_all, markers_all$avg_logFC>0)
dn_genes <- subset(markers_all, markers_all$avg_logFC<0)

# write.csv(up_genes, file="~/desktop/marker_up_genes.csv")
# write.csv(dn_genes, file="~/desktop/marker_dn_genes.csv")

df_1 <- as.data.frame(table(up_genes$cluster))
colnames(df_1) <- c("Clusters", "up")
df_2 <- as.data.frame(table(dn_genes$cluster))
colnames(df_2) <- c("Clusters", "dn")
df_3 <- dplyr::full_join(df_1, df_2, by="Clusters")


df_long <- melt(df_3, na.rm=F, id.vars=c("Clusters"))


g <- ggplot(data=df_long, aes(x=Clusters, y=value, fill=variable)) +
  geom_bar(stat="Identity", position="dodge", width=0.5, alpha=1) +
  scale_fill_manual(values=c("maroon1", "dodgerblue"),
                    name="DEGs", labels=c("Up", "Dn")
  ) +
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(size = 1),
        axis.ticks.length=unit(0.1, "cm"),
        panel.border = element_rect(colour = "black", size = 1), 
        plot.title = element_text(hjust = 0.5, vjust = 0, size=15),
        axis.text.x = element_text(colour = "black", size=15, angle=0, hjust=0.5),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15)) +
  ylab("Number of Genes") + xlab("Clusters") +
  theme(aspect.ratio = 0.5)
# scale_y_continuous(limits = c(0, 2000),
#                    expand = c(0, 0), position = "top") 
dev.off()
ggsave("marker_genes_each_clu_summary.pdf", plot = g, device = "pdf", path = "~/Desktop/",
       scale = 0.6, width = 8, height = 4, units = c("in"),
       dpi = 300, limitsize = FALSE)
