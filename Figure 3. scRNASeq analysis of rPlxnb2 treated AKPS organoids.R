######## LOAD PACKAGES AND FUNCTIONS ##############
library(dplyr)
library(Seurat)
library(sctransform)
library(fgsea)
library(msigdbr)
library(rstatix)
library(stringr)
library("Nebulosa")
library(MetBrewer)
library(MoMAColors)


msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(Hallmarks)
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(BP)
preranked_BP <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=avg_log2FC)
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


#######LOAD 10X FILES AND PREPROCESS SEURAT OBJECTS####
AKPS_rPlxnb2 <- Read10X(data.dir = "/rPlxnb2/filtered_feature_bc_matrix/")
AKPS_ctrl <- Read10X(data.dir = "/ctrl/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
AKPS_ctrl <- CreateSeuratObject(counts = AKPS_ctrl, project = "AKPS_ctrl", min.cells = 3, min.features = 200)
AKPS_rPlxnb2 <- CreateSeuratObject(counts = AKPS_rPlxnb2, project = "AKPS_rPlxnb2", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
AKPS_ctrl[["percent.mt"]] <- PercentageFeatureSet(AKPS_ctrl, pattern = "^mt-")
AKPS_rPlxnb2[["percent.mt"]] <- PercentageFeatureSet(AKPS_rPlxnb2, pattern = "^mt-")

VlnPlot(AKPS_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AKPS_rPlxnb2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

AKPS <- merge(AKPS_ctrl, AKPS_rPlxnb2, add.cell.ids = c("ctrl", "rPlxnb2"))
plot1 <- FeatureScatter(AKPS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AKPS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(AKPS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

AKPS <- subset(AKPS, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 25)

AKPS <- NormalizeData(AKPS, normalization.method = "LogNormalize", scale.factor = 10000)
AKPS <- FindVariableFeatures(AKPS, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AKPS), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AKPS)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(AKPS)
AKPS <- ScaleData(AKPS, features = all.genes, vars.to.regress = "percent.mt")
AKPS <- RunPCA(AKPS, features = VariableFeatures(object = AKPS))
DimPlot(AKPS, reduction = "pca")
ElbowPlot(AKPS)
AKPS <- FindNeighbors(AKPS, dims = 1:10)
AKPS <- FindClusters(AKPS, resolution = 0.2)
AKPS <- RunUMAP(AKPS, dims = 1:10)
Idents(AKPS) <- "seurat_clusters"
DimPlot(AKPS, split.by = "orig.ident")
Idents(AKPS) <- "orig.ident"
DimPlot(AKPS)
DEGs <- FindMarkers(AKPS, ident.1 = "AKPS_rPlxnb2", min.pct = 0.25)
View(DEGs)
BP_rPlxnb2_vs_ctrl <- preranked_BP(DEGs)
View(BP_rPlxnb2_vs_ctrl)
Hallmarks_rPlxnb2_vs_ctrl <- preranked_Hallmark(DEGs)
View(Hallmarks_rPlxnb2_vs_ctrl)

Idents(AKPS) <- "seurat_clusters"
AKPS.markers <- FindAllMarkers(AKPS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(AKPS.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC))

#remove cluster 3 and 4 because of mito and contamination, respectively
AKPS_sub <- subset(AKPS, idents = c(0, 1, 2))
all.genes <- rownames(AKPS_sub)
AKPS_sub <- ScaleData(AKPS_sub, features = all.genes, vars.to.regress = "percent.mt")
AKPS_sub <- RunPCA(AKPS_sub, features = VariableFeatures(object = AKPS_sub))
DimPlot(AKPS_sub, reduction = "pca")
ElbowPlot(AKPS_sub)
AKPS_sub <- FindNeighbors(AKPS_sub, dims = 1:10)
AKPS_sub <- FindClusters(AKPS_sub, resolution = 0.15)
AKPS_sub <- RunUMAP(AKPS_sub, dims = 1:10)
Idents(AKPS_sub) <- "seurat_clusters"
current.cluster.ids <- c(0, 1, 2)
DimPlot(AKPS_sub, cols=met.brewer("Isfahan2", 3)) + ggsave("Clusterdimplot.pdf", width =6,  height=5)
Idents(AKPS_sub) <- "orig.ident"
DimPlot(AKPS_sub, reduction = "umap", group.by = "orig.ident", cols=met.brewer("Egypt", 2)) + ggsave("Sampledimplot.pdf", width =6,  height=5)



####COMPOSITIONAL ANALYSIS####
#frequencies per cluster
numberofcells         <- table(AKPS_sub$orig.ident, AKPS_sub$seurat_clusters)
numberofcells
totalcellssample   <- c(sum(numberofcells[1,]), sum(numberofcells[2,]))
a                     <- cbind(numberofcells,totalcellssample)
a
totalcellspercluster  <- c(sum(a[,1]), sum(a[,2]), sum(a[,3]), sum(a[,4]))
b                     <- rbind(a, totalcellspercluster)
b
b[1:3,1]
c0 <- (b[1:2,1]/totalcellssample)*100
c1 <- (b[1:2,2]/totalcellssample)*100
c2 <- (b[1:2,3]/totalcellssample)*100

c <- rbind(c0,c1,c2)
colSums(c)

rownames(c) =  rev(c("cluster 0", "cluster 1", "cluster 2"))
c

#plot
par(mar=c(6,8,2,14))
pdf(file="Clusterbreakdown.pdf")
barplot(c, horiz=TRUE,
        legend = T, border=NA,
        args.legend=list(bty = "n",x=130, cex=.8),
        main = "Cluster breakdown per sample", 
        las = 1, 
        col= met.brewer("Isfahan2", 3))
dev.off()
display.all.moma(n=3, override_order = T, direction = 1)
library(MetBrewer)
display.all.met(n=3, override_order = T, direction = 1)




######MARKERS and GO ANALYSIS######
Idents(AKPS_sub) <- "seurat_clusters"


DEGs_C0 <- FindMarkers(AKPS_sub, ident.1 = 0, min.pct = 0.25, features = genes.use)
View(DEGs_C0)
BP_C0 <- preranked_BP(DEGs_C0)
View(BP_C0)

ggplot(BP_C0 %>% filter(padj<0.05) %>% head(n= 100), aes(reorder(tolower(pathway), NES), NES)) +
  geom_point(aes(size=size, colour=padj)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(y="Normalized Enrichment Score", x= " ", title="Biological processes in Cluster 0") + 
  theme_classic()+
  theme(axis.text.y=element_text(size=10))# + ggsave("PLxnb2OEvsNT_Hallmarks.pdf", width =7,  height=5)

DEGs_C1 <- FindMarkers(AKPS_sub, ident.1 = 1, min.pct = 0.25, features = genes.use)
View(DEGs_C1)
BP_C1 <- preranked_BP(DEGs_C1)
View(BP_C1)

ggplot(BP_C1 %>% filter(padj<0.05) %>% head(n= 100), aes(reorder(tolower(pathway), NES), NES)) +
  geom_point(aes(size=size, colour=padj)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(y="Normalized Enrichment Score", x= " ", title="Biological processes in Cluster 1") + 
  theme_classic()+
  theme(axis.text.y=element_text(size=10))

DEGs_C2 <- FindMarkers(AKPS_sub, ident.1 = 2, min.pct = 0.25, features = genes.use)
View(DEGs_C2)
BP_C2 <- preranked_BP(DEGs_C2)
View(BP_C2%>% filter(padj<0.05))


########### EMT SCORES ##########
Idents(AKPS_sub) <-"orig.ident"
DotPlot(AKPS_sub, features=c("Cdh1", "Vim", "Zeb1", "Klf4"))
DotPlot(AKPS_sub, features = c("Epcam", "Cdh1", 
                               "Lgr5", "Axin2", "Mki67", 
                               "Zeb1", "Klf4", "Id3", "Foxp2",
                               "Ovol1", "Grhl2", "Elf3"))

stemness_list<- list(c("Lgr5", "Axin2",  "Ascl2",  "Tnfrsf19"))
AKPS_sub <-AddModuleScore(AKPS_sub, features= stemness_list,name = "Stem")
names(x = AKPS_sub[[]])

plot_density(AKPS_sub, features = "Stem1")  +labs(title = "Stemness (Lgr5, Axin2, Ascl2, Tnfrsf19)") 

apoptosis.score <- list(c(BP$GOBP_APOPTOTIC_SIGNALING_PATHWAY))
AKPS_sub <-AddModuleScore(AKPS_sub, features= apoptosis.score,name = "Apoptosis")
names(x = AKPS_sub[[]])


Idents(AKPS_sub) <- "orig.ident"
AKPS_sub_Pl <- subset(AKPS_sub, idents = "AKPS_rPlxnb2")
AKPS_sub_ctrl <- subset(AKPS_sub, idents = "AKPS_ctrl")
p<-wilcox.test(AKPS_sub_Pl$EMTmoor1, AKPS_sub_ctrl$EMTmoor1, alternative = "two.sided") #p-value < 2.2e-16
adp <- p.adjust(p$p.value, method = "fdr", n = length(p))
adp


























