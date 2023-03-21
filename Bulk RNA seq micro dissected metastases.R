#code for the analysis of bulk transcriptomic profiling of microdissected AKPS mets grown in Plxnb2-overexpressing or NT livers

library(rlang)
library(devtools)
library(BiocManager)
library(biomaRt)
library("dplyr")
library(ggplot2)
library(RColorBrewer) 
library(ggraph)
library(wesanderson)
library(fgsea)
library(msigdbr)
library(ggrepel)
library(cowplot)
library(magrittr)
library(data.table)
library(stringr)
library(tidyr)
library("edgeR")

load(metsinAAV.Rds)
metsinAAV <- all_samples_mat[,c(2,8,16,18,20,22,25)]
colnames(metsinAAV) <- c("AMO25713", "AMO25715", "AMO25297", "AMO25712",  "AMO25294", "AMO25856", "AMO25934")
#Create DGEList object
group <- factor(c(2,2,1,2,1,2,1)) #1 is NT 2 is OE
y <- DGEList(counts=metsinAAvno710,group=group)
y$samples
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y), las=2, main = "Barplot of library sizes")

# Filter reads by counts: Most of the samples should have at least 10 reads, normalize the library and estimate dispersion
#keep <- filterByExpr(y, min.count = 10)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
design <- model.matrix(~group)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1)
qlf.2vs1$table$Gene <- rownames(qlf.2vs1$table)
View(qlf.2vs1$table)


# Export all genes
all_qlf.2vs1 = topTags(qlf.2vs1, n = Inf)
dim(qlf.2vs1)
head(qlf.2vs1$table)
qlf.2vs1$table$Gene <- rownames(qlf.2vs1$table)
View(qlf.2vs1$table)
write.csv(qlf.2vs1$table, file = "PLXNB2_VS_NT.csv")


#volcano plot
DEGs_volcano <- function(DEGs, p_treshold, FC_treshold, title, color, upperylim, xlim) {
  DEGs$Significant <- ifelse((DEGs$PValue < p_treshold & abs(DEGs$logFC) > FC_treshold ), "Significant", "Not Sig")
  DEGs$Gene <- ifelse((DEGs$PValue < p_treshold & abs(DEGs$logFC) > FC_treshold ), rownames(DEGs), NA)
  
  plot <-  ggplot(data=DEGs, aes(x=logFC, y=-log10(PValue), fill=factor(Significant), label = DEGs$Gene) ) + 
    theme_bw()+ 
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank())+
    geom_point(shape = 21,color=color)+
    scale_fill_manual(values = c(color, "black")) +
    geom_text_repel(size=3, max.overlaps = Inf)+
    xlab("log2 fold change") + ylim(0,upperylim)+xlim(-xlim,xlim)+
    ylab("-log10 adjusted p-value") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=12),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}
DEGs<-qlf.2vs1$table
DEGs_volcano(DEGs, 0.02, 3, "Plxnb2 vs NT", "grey",6, 8) +ggsave("mets_volcano.pdf", width =8,  height=8)


#gene set enrichment analysis
msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(Hallmarks)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
REACTOME <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(REACTOME)
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(BP)

preranked_BP <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=log10(PValue)/sign(logFC))
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  
  BP_x <- fgsea(pathways = BP, 
                stats = ranks,
                minSize=10,
                maxSize=500,
                nperm=1000000)
  
  BP_x$pathway<-gsub("GOBP_","",BP_x$pathway)
  BP_x$pathway<-gsub("_"," ",BP_x$pathway)
  return(BP_x)
}

preranked_REACTOME <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=log10(PValue)/sign(logFC))
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  
  REACTOME_x <- fgsea(pathways = REACTOME, 
                      stats = ranks,
                      minSize=10,
                      maxSize=500,
                      nperm=1000000)
  
  REACTOME_x$pathway<-gsub("REACTOME_","",REACTOME_x$pathway)
  REACTOME_x$pathway<-gsub("_"," ",REACTOME_x$pathway)
  return(REACTOME_x)
}
preranked_Hallmark <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=log10(PValue)/sign(logFC))
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  
  Hallmarks_x <- fgsea(pathways = Hallmarks, 
                       stats = ranks,
                       minSize=10,
                       maxSize=500,
                       nperm=1000000)
  
  Hallmarks_x$pathway<-gsub("HALLMARK_","",Hallmarks_x$pathway)
  Hallmarks_x$pathway<-gsub("_"," ",Hallmarks_x$pathway)
  return(Hallmarks_x)
}

BP_enrichment <- preranked_BP(qlf.2vs1$table)
View(BP_enrichment)
REACTOME_enrichment <- preranked_REACTOME(qlf.2vs1$table)
View(REACTOME_enrichment)
Hallmarks_enrichment <- preranked_Hallmark(qlf.2vs1$table)
View(Hallmarks_enrichment)
KEGG_enrichment <- preranked_KEGG(qlf.2vs1$table)
View(KEGG_enrichment)


ggplot(BP_enrichment %>% filter(abs(NES)>2.1 & pval<0.01) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=pval)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=" ") + 
  theme_classic()+
  labs(y="")+   
  theme(axis.text.y=element_text(size=10)) + ggsave("mets_BP.pdf", width =10,  height=8)



