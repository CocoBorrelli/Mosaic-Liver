Annotation_mouse <- function(zUMI_output, prefix) {
  ensembl<-useEnsembl(biomart="ensembl")
  list<-listDatasets(ensembl)
  mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",version = 95)
  attributes<-listAttributes(mart)
  gene_ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart)
  dataframe<-as.data.frame(as.matrix(zUMI_output$umicount$exon$all))%>%
    as.matrix(.)
  colnames(dataframe)<-paste(colnames(dataframe), prefix,sep = "_") 
  dataframe<-mutate(as.data.frame(dataframe),ensembl_gene_id_version=rownames(dataframe))
  join<-dataframe%>%
    left_join(dplyr::select(gene_ids,1:2))
  length(unique(join$external_gene_name))
  join<-join[!duplicated(join$external_gene_name),]
  join[is.na(join)]<-0 #make all empty value to zero
  rownames(join)<-join$external_gene_name
  join<-dplyr::select(join,-ensembl_gene_id_version,-external_gene_name)
  return(join)
}

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

#read in dgecounts.rds file 
poolA <- readRDS("/media/Coco/MOSAIC LIVER/Experiments/shSema4s in AKPS/bulk_orgs/poolA/zUMIs_output/expression/bulk_orgs_round2.dgecounts.rds")
colnames(poolA$umicount$exon$all)
poolA_ann <- Annotation_mouse(poolA, "poolA")
colnames(poolA_ann) <- c("mets in AMO25713", "AKPS sLPmCherry 3","KPC parental", 
                         "AKPS shRNA Sema4s 2", "KPC shRNA EV 1", "KPC shRNA Sema4s 2", 
                         "mets in AMO25715", "KPC shRNA EV 2", "KPC shRNA EV 3")
poolA_ann$Gene <- rownames(poolA_ann)

poolB <- readRDS("/media/Coco/MOSAIC LIVER/Experiments/shSema4s in AKPS/bulk_orgs/poolB/zUMIs_output/expression/bulk_orgs_round2.dgecounts.rds")
colnames(poolB$umicount$exon$all)
poolB_ann <- Annotation_mouse(poolB, "poolB")
colnames(poolB_ann) <- c("AKPS sLPmCherry 2", 
                         "AKPS shRNA EV 3", 
                         "KPC shRNA Sema4s 1",
                         "AKPS shRNA EV 1", 
                         "AKPS shRNA Sema4s 1", 
                         "mets in AMO25297", 
                         "AKPS shRNA Sema4s 3", 
                         "mets in AMO25712", 
                         "mets in AMO25710")      
poolB_ann$Gene <- rownames(poolB_ann) 


poolC <- readRDS("/media/Coco/MOSAIC LIVER/Experiments/shSema4s in AKPS/bulk_orgs/poolC/zUMIs_output/expression/bulk_orgs_round2.dgecounts.rds")
colnames(poolC$umicount$exon$all)
poolC_ann <- Annotation_mouse(poolC, "poolC")
View(poolC_ann)
colnames(poolC_ann) <- c("mets in AMO25294" , "AKPS shRNA EV 2", "mets in AMO25856",
                                        "mets derived org line 1",  "mets derived org line 2",
                                        "mets in AMO25934" ,"AKPS sLPmCherry 1")
poolC_ann$Gene <- rownames(poolC_ann)


#create sample count matrix
poolAB <- merge(poolA_ann, poolB_ann,by="Gene",  stringsAsFactors = TRUE)
colnames(poolAB)
poolABC <- merge(poolAB, poolC_ann,by="Gene",  stringsAsFactors = TRUE)
colnames(poolABC)

all_samples_mat <- as.matrix(poolABC[2:26])
rownames(all_samples_mat)<- poolABC$Gene
metsinAAv <- all_samples_mat[,c(2,8,16,18,19,20,22,25)]
colnames(metsinAAv) <- c("AMO25713", "AMO25715", "AMO25297", "AMO25712", "AMO25710", "AMO25294", "AMO25856", "AMO25934")
metsinAAvno710 <- all_samples_mat[,c(2,8,16,18,20,22,25)]
colnames(metsinAAvno710) <- c("AMO25713", "AMO25715", "AMO25297", "AMO25712",  "AMO25294", "AMO25856", "AMO25934")


#barcode to sample assignment
#barcode          sample                group
# "ACAGTG_poolA"  mets in AMO25713          1
#"AGTTCC_poolA"   AKPS sLPmCherry 3
#"ATCACG_poolA"   KPC shRNA EV 1
#"ATTCCT_poolA"   KPC parental
# "CACTCA_poolA"  AKPS shRNA Sema4s 2       3
# "CTATAC_poolA"  KPC shRNA Sema4s 2
#"GCCAAT_poolA"   mets in AMO25715          1
#"TACAGC_poolA"   KPC shRNA EV 2
#"TATAAT_poolA"   KPC shRNA EV 3
#"AGTCAA_poolB"   AKPS sLPmCherry 2
#"CCAACA_poolB"   AKPS shRNA EV 3           4
#"CTAGCT_poolB"   KPC shRNA Sema4s 1
#"GTAGAG_poolB"   AKPS shRNA EV 1           4
#"TAATCG_poolB"   AKPS shRNA Sema4s 1       3
#"TCATTC_poolB"   mets in AMO25297          2
#"TCGAAG_poolB"   AKPS shRNA Sema4s 3       3
# "TGACCA_poolB"  mets in AMO25712          1
# "TTAGGC_poolB"  mets in AMO25710          2
# "ATCACG_poolC"  mets in AMO25294          2
# "ATGAGC_poolC"  AKPS shRNA EV 2           4
#"CAGATC_poolC"   mets in AMO25856          1
#"CCGTCC_poolC"   mets derived org line 1
#"GACGAC_poolC"   mets derived org line 2
#"GAGTGG_poolC"   mets in AMO25934          2
#"GTGAAA_poolC"   AKPS sLPmCherry 1

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
#y <- estimateDisp(y)
y <- estimateDisp(y,design)
plotMDS(y)
logcpm <- cpm(y, log=TRUE)
View(logcpm)
#write.table(logcpm, file="logcpm_OE_NT.txt")

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


#functions
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
DEGs_volcano(DEGs, 0.02, 3, "Plxnb2 vs NT", "grey",6, 8) #+ggsave("mets_volcano.pdf", width =8,  height=8)


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



#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.zora.uzh.ch/id/eprint/120377/1/2403.full.pdf
#https://www.nature.com/articles/nrm1740



###IMPORT DATASETS ####
msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(Hallmarks)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
#CGP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#head(CGP)
#m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
#KEGG <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#head(KEGG)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
REACTOME <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(REACTOME)
#m_df<- msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT")
#TFT <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#head(TFT)
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(BP)
#m_df<- msigdbr(species = "Mus musculus", category = "C6")
#oncSig <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#head(oncSig)



