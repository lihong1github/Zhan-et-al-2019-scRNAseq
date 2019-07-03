# Author: Lihong Zhan
# Date: 06/05/2018
# The following script is used to make bar graph based on GESA data.

library(stringr)
library(dplyr)
library(ggplot2)


# set the directory where your csv file is located
setwd("~/desktop/gsea_input")
files <- list.files()

data_conslidate <- data.frame()

for (i in 1:length(files)) {
  i=2
  df <- read.table(files[i], sep="\t", header=F, fill=T, stringsAsFactors = FALSE)
  df <- data.frame(df[4:13, 1:7])
  colnames(df) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR")
  
  df$FDR <- as.numeric(df$FDR)
  
  BP_1 <- dplyr::mutate(df, minus_logp=0-log10(FDR))
  BP_1 <- dplyr::mutate(BP_1, enrich=paste(BP_1$genes_in_data, "/", BP_1$genes_in_gsea, sep=" ") )
  BP_2 <- BP_1[order(-BP_1$minus_logp),]
  BP_2$Name <-  gsub("HALLMARK_", "", BP_2$name)
  BP_2$Name <- gsub("*_", " ", BP_2$Name)
  BP_2$Name <-  factor(BP_2$Name, levels=rev(BP_2$Name))
  
  df_summary <- data.frame(name=BP_2$Name, cluster=files[i])
  data_conslidate <- rbind(data_conslidate, df_summary)
  
  ggplot(data=BP_2, aes(x=Name  , y=minus_logp)) +
    theme_bw()+
    theme(panel.grid.major.x  = element_line(size =0.2, colour="grey88"),
          panel.grid.major.y  = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(colour = "black", size=12),
          axis.text.y = element_text(colour = "black", size=12),
          axis.title.x = element_text(colour = "black", size=14),
          axis.title.y = element_text(colour = "black", size=14)) +
    ylab("-Log(FDR)") + xlab(NULL) +
    # change Title to your desired title, \n is used to here to return the line
    ggtitle(paste(files[i], "\nGESA:Hallmark")) +
    geom_bar(stat="Identity",  width=0.7, fill="maroon1", alpha=0.5) +
    # maroon1 dodgerblue
    geom_text(aes(label=enrich), vjust=0.4, 
              hjust=0, size=4, color="black", 
              stat="Identity", y=0.05*max(BP_2$minus_logp)) +
    coord_flip() + 
    theme(aspect.ratio = 1.5)+
    scale_y_continuous(limits = c(0,max(BP_2$minus_logp)+1),
                       expand = c(0, 0), position = "top") 
  
  ggsave(paste(files[i], ".pdf", sep=""), plot = last_plot(), device = "pdf", path = "~/Desktop/",
         scale = 0.5, width = 12, height = 6, units = c("in"),
         dpi = 600, limitsize = FALSE)

}


