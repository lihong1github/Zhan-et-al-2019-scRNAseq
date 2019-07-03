library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringr)

##read in the data - note that it should be in .csv format

df_input <- read.csv(file="~/desktop/markers_UMI_cluster-4.csv", header=T, stringsAsFactors = F)


df <- data.frame(gene=df_input$gene,
                 Log2FC=df_input$avg_logFC, 
                 FDR=df_input$p_val_adj, stringsAsFactors=FALSE)

df$colours <- c("NC")

df$colours[df$Log2FC >= 1 & df$FDR <= 0.05] <- c("UP")
df$colours[df$Log2FC <= -1 & df$FDR <= 0.05] <- c("DN")


# label particular genes
genes_to_plot_1 <- df_input$gene[1:11]
genes_to_plot_1 <- genes_to_plot_1[-c(4, 6)]

df$colours[df$gene %in% genes_to_plot_1] <- c("GOI")
df$select_gene[df$gene %in% genes_to_plot_1] <- df$gene[df$gene %in% genes_to_plot_1]


table(df$colours)


df$colours <- factor(df$colours, levels= c("NC", "UP","DN","GOI"))
my_color <- c("gray60","deeppink", "royalblue1","black")


##draw the plot
dev.off()

ggplot(data=df, aes(x=Log2FC, y=-log(FDR, base = 10), colour=colours)) + 
  geom_point(shape=19, alpha=1, size=1) + 
  scale_colour_manual(values = my_color,
                      name="Legend",
                      breaks=c("NC", "UP","DN", "GOI"),
                      labels=c("NC", "UP", "DN", "GOI"))+
  geom_vline(xintercept =  1,linetype="dashed", colour = "grey36") +
  geom_vline(xintercept = -1,linetype="dashed", colour = "grey36") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", colour = "grey36") +
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=20),
        axis.text.y = element_text(colour = "black", size=20),
        axis.title.x = element_text(colour = "black", size=20),
        axis.title.y = element_text(colour = "black", size=20)) +
  ylab("-Log10[FDR]") + xlab("Log2FC[Cluster-4 vs Other Clusters]") +
  geom_point(data=df[df$colours=="GOI",],
             aes(x=`Log2FC`, y=-log(FDR, base = 10)),
             shape=21, alpha=1, size=2, fill="black",color="white") +
  geom_text_repel(aes(x=`Log2FC`, y=-log10(FDR), label = select_gene), 
                  color="black", fontface = 'bold',size = 5, box.padding = 0.4,
                  point.padding = 0.5, segment.size=0.5, segment.colour="gold") +
  theme(aspect.ratio = 1) 
  scale_y_continuous(breaks=seq(0, 400, 100), limits=c(0, 400)) +
  scale_x_continuous(breaks=seq(-4, 4, 2), limits=c(-4, 4))

ggsave("Volcano_Cluster-4.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/",
       scale = 0.8, width = 7, height = 7, units = c("in"),
       dpi = 600, limitsize = FALSE)
