library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(scales)

load(file="~/desktop/PLX_Rdata/Data_RunTSNE.Rdata")


cluster_number <- table(experiment.aggregate@ident,experiment.aggregate@meta.data$orig.ident)
cluster_number <- as.data.frame.matrix(t(cluster_number))


cluster_percent_summary <- data.frame()
for (i in 1:9) {
cluster_percent <- rbind(cluster_number[,i]/rowSums(cluster_number)*100)
cluster_percent_summary <- rbind(cluster_percent_summary, cluster_percent)}

cluster_percent_relative_summary <- data.frame()
for (i in 1:9) {
  cluster_percent_relative <- rbind(cluster_percent_summary[i,]/rowSums(cluster_percent_summary)[i]*100)
  cluster_percent_relative_summary <- rbind(cluster_percent_relative_summary, cluster_percent_relative)}

df <- as.data.frame(t(cluster_percent_relative_summary))
clu <- paste("Cluster", c(1:9), sep="-")

sample_color = c("azure4", "Deeppink1", "royalblue1")
my_color <- hue_pal()(9)


# get the cell number for each cluster
cluster_cellnum <- as.data.frame(table(experiment.aggregate@ident))

df$Sample <- factor(c("Ctrl", "D0", "D2"), levels=c("Ctrl", "D0", "D2"))

class(df$Sample)
for (i in 1:9) {
  my_labels <- paste(df$Sample, sprintf("%.2f",round(df[, i], 2)), sep=": ")
  my_labels <- paste(my_labels, "%", sep="")
  ggplot(df) +
    geom_bar(aes(x="", y=df[,i], fill=Sample), width = 1, stat = "identity", color="white", size=1) + 
    coord_polar("y", start=0) +
    scale_fill_manual(values = sample_color,
                      name="Relative Freq",
                      breaks=c("Ctrl", "D0", "D2"),
                      labels=my_labels) +
    geom_text(aes(x=0, y=0), label=cluster_cellnum$Freq[i], size=8) +
    theme_minimal() +
    ggtitle(clu[i]) +
    ylab(NULL) + xlab(NULL) +
    theme(aspect.ratio = 1) +
    theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"),
          legend.title = element_text(colour="black", size=14, hjust=0.5),
          legend.text = element_text(colour="black", size=14),
          legend.position='bottom',
          legend.direction = "vertical",
          axis.text.x = element_text(colour = "black", size=0),
          panel.grid.major.x  = element_line(size = 2, colour=my_color[i]),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank())
  
  ggsave(paste(clu[i], "pie.pdf", sep="_"), plot = last_plot(), device = "pdf", path = "~/Desktop/pie_charts",
         scale = 0.8, width = 4, height = 4, units = c("in"),
         dpi = 600, limitsize = FALSE)
  
}
