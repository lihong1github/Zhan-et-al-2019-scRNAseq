source("https://bioconductor.org/BiocManager.R")
# biocLite("GSEABase")
# biocLite("rWikiPathways")
BiocManager::install(c("GSEABase", "rWikiPathways"))
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"))

library(rWikiPathways)
library(XML)
# install.packages("tidyverse")
library(tidyr)
library(dplyr)
library(clusterProfiler)
library(GSEABase)
library(org.Mm.eg.db)

# download and prepare pathway data from wikipathway
wp.mm.gmt <- downloadPathwayArchive(organism="Mus musculus", format = "gmt")

# convert entries to GMT format
wp2gene <- clusterProfiler::read.gmt(wp.mm.gmt)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")

wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME

# prepare my data
df <- read.csv(file="~/desktop/markers_UMI_cluster-Mac2_pos.csv", header=T, stringsAsFactors = F)

# select on upregulated genes

  query_genes <- df$gene[df$avg_logFC>0]
  query_genes_1 <- bitr(query_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  wp <- clusterProfiler::enricher(
    query_genes_1[[2]],
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05, #p.adjust cutoff
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
  
  # head(wp, 15)
  
  # barplot(wp, showCategory = 20)
  
  BP_2 <- head(wp,10)
  BP_2$Minus_logFDR <- -log(BP_2$p.adjust, 10)
  BP_2 <- BP_2 %>%
    arrange(Minus_logFDR)
  BP_2$Description <- factor(BP_2$Description, levels = BP_2$Description)
  
  
  title = "Mac2+ UP DEGs\nWikipathway"
  
  ggplot(data=BP_2, aes(x=Description, y=Minus_logFDR)) +
    geom_point(stat="Identity",aes(size=Count), alpha=0.8,
               colour="white",pch=21, fill="maroon1") +
    # maroon1 dodgerblue
    theme_bw()+
    theme(panel.background = element_rect(fill="gray40"),
          panel.grid.major.x  = element_line(size =0.25, colour="grey88"),
          panel.grid.major.y  = element_line(size =0.25, colour="grey88"),
          panel.grid.minor = element_line(size =0.1, colour="grey88"),
          plot.title = element_text(hjust = 0.5, size=14),
          axis.text.x = element_text(colour = "black", size=12),
          axis.text.y = element_text(colour = "black", size=12),
          axis.title.x = element_text(colour = "black", size=14),
          axis.title.y = element_text(colour = "black", size=14),
          legend.position="none") +
    ylab("-Log(FDR)") + xlab(NULL) + ggtitle(title) +
    geom_text_repel(aes(x=Description, y=Minus_logFDR, label = GeneRatio),
                    color="white", fontface = 'italic',size = 4, box.padding = 1,
                    point.padding = 1, segment.size=0.2, segment.colour="white") +
    coord_flip() + 
    theme(aspect.ratio = 1.5) 
  
  ggsave("Mac2+UP_WK.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/",
         scale = 0.5, width = 12, height = 6, units = c("in"),
         dpi = 600, limitsize = FALSE)
  

