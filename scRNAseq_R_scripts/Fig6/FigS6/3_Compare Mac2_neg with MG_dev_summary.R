library(dplyr)
library(ggplot2)
library(RColorBrewer)

df <- read.csv(file="~/desktop/MG_dev.csv", stringsAsFactors = F, header = T)
df_mac2 <- read.csv(file="~/desktop/markers_UMI_cluster-Mac2_neg.csv", stringsAsFactors = F, header = T)

df_mac2_up <- df_mac2$gene[df_mac2$avg_logFC>0]
df_mac2_dn <- df_mac2$gene[df_mac2$avg_logFC<0]

length(df_mac2_up)
length(df_mac2_dn)

dev_stage <- c("YS", "E1", "E2", "P1", "P2", "A1", "A2")

percent_summary <- data.frame()
percent_each <- numeric()
for (i in 1:length(dev_stage)) {
  dev_gene <- df$Gene[df$Cluster==i]
  dev_gene_up <- dev_gene[dev_gene %in% df_mac2_up]
  percent_up <- length(dev_gene_up)/length(dev_gene)*100
  dev_gene_dn <- dev_gene[dev_gene %in% df_mac2_dn]
  percent_dn <- length(dev_gene_dn)/length(dev_gene)*100
  percent_each <- c(percent_up, percent_dn)
  percent_summary <- rbind(percent_summary,percent_each)
}
colnames(percent_summary) <- c("Percent_UP", "Percent_DN")
percent_summary$dev_stage <- dev_stage
write.csv(percent_summary, file="~/desktop/percent_summary_each_stage.csv", row.names = F)

# calculates % of genes in mac2_up that belong to each dev_stage
percent_summary_mac2 <- data.frame()
percent_each_mac2 <- numeric()
for (i in 1:length(dev_stage)) {
  dev_gene <- df$Gene[df$Cluster==i]
  dev_gene_up <- df_mac2_up[df_mac2_up %in% dev_gene]
  percent_up <- length(dev_gene_up)/length(df_mac2_up)*100
  dev_gene_dn <- df_mac2_dn[df_mac2_dn %in% dev_gene]
  percent_dn <- length(dev_gene_dn)/length(df_mac2_dn)*100
  percent_each_mac2 <- c(percent_up, percent_dn)
  percent_summary_mac2 <- rbind(percent_summary_mac2,percent_each_mac2)
}
unmatched <- c((100-colSums(percent_summary_mac2[1])), (100-colSums(percent_summary_mac2[2])))
percent_summary_mac2 <- rbind(percent_summary_mac2, unmatched)

colnames(percent_summary_mac2) <- c("Percent_UP_MAC2_pos", "Percent_DN_MAC2_pos")
percent_summary_mac2$dev_stage <- c(dev_stage, "Others")
percent_summary_mac2$dev_stage <- factor(percent_summary_mac2$dev_stage, levels = percent_summary_mac2$dev_stage)


write.csv(percent_summary_mac2, file="~/desktop/percent_summary_each_stage_MAC2_pos.csv", row.names = F)

### make graph
# my_color <- c("#E74C3C", 
#               "#2980B9", "#2980B999", 
#               "#3498DB", "#3498DB99", 
#               "#2C3E50", "#2C3E5099",
#               "#ECF0F1")

my_color <- rev(brewer.pal(8,"OrRd"))[1:8]
# adjustcolor( "#2C3E50", alpha.f = 0.6)
my_labels <- paste(percent_summary_mac2$dev_stage, sprintf("%.2f",round(percent_summary_mac2$Percent_UP_MAC2_pos, 2)), sep=": ")
my_labels <- paste(my_labels, "%", sep="")


ggplot(percent_summary_mac2, aes(x=2, y=percent_summary_mac2$Percent_UP_MAC2_pos, fill=percent_summary_mac2$dev_stage)) +
  geom_bar(width = 1, stat = "identity", color="black", size=0.25) + coord_polar("y", start=0, direction = -1) +
  xlim(0.5, 2.5) + scale_fill_manual(values = my_color,
                    name="MG Dev Genes", labels=my_labels) +
  ggtitle("Mac2- UP DEGs") +
  geom_text(aes(x=0.5, y=0), label="UP\n191", size=8) +
  ylab(NULL) + xlab(NULL) +
  theme(aspect.ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=14),
        legend.position='right',
        legend.direction = "vertical",
        axis.text.x = element_text(colour = "black", size=0),
        panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank()) 

ggsave("Mac2_neg_UP_devstages.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 0.8, width = 8, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)


#################

percent_summary_mac2_reorder <- percent_summary_mac2[c(7:1, 8),]
percent_summary_mac2_reorder$dev_stage <- factor(percent_summary_mac2_reorder$dev_stage, levels=percent_summary_mac2_reorder$dev_stage)
my_color <- rev(brewer.pal(8,"GnBu"))

my_labels <- paste(percent_summary_mac2_reorder$dev_stage, sprintf("%.2f",round(percent_summary_mac2_reorder$Percent_DN_MAC2_pos, 2)), sep=": ")
my_labels <- paste(my_labels, "%", sep="")

ggplot(percent_summary_mac2_reorder, aes(x=2, y=percent_summary_mac2_reorder$Percent_DN_MAC2_pos, fill=percent_summary_mac2_reorder$dev_stage)) +
  geom_bar(width = 1, stat = "identity", color="black", size=0.25) + coord_polar("y", start=0, direction = -1) +
  xlim(0.5, 2.5) + scale_fill_manual(values = my_color, breaks=percent_summary_mac2_reorder$dev_stage,
                                     name="MG Dev Genes", label=my_labels) +
  ggtitle("Mac2- DN DEGs") +
  geom_text(aes(x=0.5, y=0), label="DN\n150", size=8) +
  ylab(NULL) + xlab(NULL) +
  theme(aspect.ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=14),
        legend.position='right',
        legend.direction = "vertical",
        axis.text.x = element_text(colour = "black", size=0),
        panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())

ggsave("Mac2_neg_DN_devstages.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 0.8, width = 8, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)
