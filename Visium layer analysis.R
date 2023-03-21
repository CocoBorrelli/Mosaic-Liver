### This code by Kristina Handler identifies liver layers neighbouring metastatic foci, here applied to a sample of AKPS hepatic lesions

### load libraries 
library(devtools)
library(STutility)
library(spdep)
library(kableExtra)
library(NNLM)
library(RColorBrewer)

### set working directory 
setwd("./Mets6M1")

### load data 
samples <- "./filtered_feature_bc_matrix"
spotfiles <- "./spatial/tissue_positions_list.csv"
imgs <- "./spatial/tissue_hires_image.png"
json <- "./spatial/scalefactors_json.json"
df_paths6M1 <- data.frame(samples,spotfiles,imgs,json)
se6M1 <- InputFromTable(infotable = df_paths6M1)

### load H&E image to overlay
se6M1 <- LoadImages(se6M1, time.resolve = FALSE, verbose = TRUE)
ImagePlot(se6M1, method = "raster", type = "raw")

### Normalization 
se6M1 <- SCTransform(se6M1)

### apply nFeature cutoff  
se6M1_sub <- SubsetSTData(se6M1, expression = nFeature_RNA > 200)

### Matrix factorization 
se6M1_2 <- RunNMF(se6M1_sub) # Specificy nfactors to choose the number of factors, default=20.

FeatureOverlay(se6M1_2, features = c("Gpx2"), 
               sampleids = 1,
               cols = c("lightgray", "mistyrose", "red", "darkred", "black"),
               pt.size = 3, 
               add.alpha = TRUE,
               ncol = 2)

### Clustering 
se6M1_2 <- FindNeighbors(object = se6M1_2, verbose = FALSE, reduction = "NMF", dims = 1:10)
se6M1_2 <- FindClusters(object = se6M1_2, verbose = FALSE, resolution = 0.3)
ST.FeaturePlot(object = se6M1_2, features = "seurat_clusters", pt.size = 3, ncol = 2)
# Run UMAP
se6M1_2 <- RunUMAP(se6M1_2, reduction = "NMF", dims = 1:10, n.neighbors = 10)
DimPlot(se6M1_2, reduction = "umap", label = TRUE)

### Identify the metastasis cluster = 0
FeaturePlot(se6M1_2, features = "Gpx2")

### subset mets positive regions and go spot by spot 
##first layer
se6M1_3 <- SetIdent(se6M1_2, value = "seurat_clusters")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "0", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_0", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
ST.FeaturePlot(object = se6M1_3, features = "nbs_0", pt.size = 3, ncol = 2)

#combine 0 and nbs_0 and call it 1 wich increase the mets neighbouring area and replace NA to 2 which is the not neighbour areas in our case 
current.cluster.ids <- c("0","nbs_0")
new.cluster.ids <- c(1,1)
se6M1_3$layer <- plyr::mapvalues(x = se6M1_3$nbs_0, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer"
FeatureOverlay(se6M1_3, features = "layer", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 2
ST.FeaturePlot(object = se6M1_3, features = "layer", pt.size = 3, ncol = 2)

##second layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "1", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_1", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("1","nbs_1")
new.cluster.ids <- c(2,2)
se6M1_3$layer2 <- plyr::mapvalues(x = se6M1_3$nbs_1, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer2"
FeatureOverlay(se6M1_3, features = "layer2", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 3
ST.FeaturePlot(object = se6M1_3, features = "layer2", pt.size = 3, ncol = 2)

##third layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer2")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "2", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_2", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("2","nbs_2")
new.cluster.ids <- c(3,3)
se6M1_3$layer3 <- plyr::mapvalues(x = se6M1_3$nbs_2, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer3"
FeatureOverlay(se6M1_3, features = "layer3", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 4
ST.FeaturePlot(object = se6M1_3, features = "layer3", pt.size = 3, ncol = 2)

##forth layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer3")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "3", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_3", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("3","nbs_3")
new.cluster.ids <- c(4,4)
se6M1_3$layer4 <- plyr::mapvalues(x = se6M1_3$nbs_3, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer4"
FeatureOverlay(se6M1_3, features = "layer4", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 5
ST.FeaturePlot(object = se6M1_3, features = "layer4", pt.size = 3, ncol = 2)

##fifth layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer4")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "4", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_4", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("4","nbs_4")
new.cluster.ids <- c(5,5)
se6M1_3$layer5 <- plyr::mapvalues(x = se6M1_3$nbs_4, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer5"
FeatureOverlay(se6M1_3, features = "layer5", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 6
ST.FeaturePlot(object = se6M1_3, features = "layer5", pt.size = 3, ncol = 2)

##sixth layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer5")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "5", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_5", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("5","nbs_5")
new.cluster.ids <- c(6,6)
se6M1_3$layer6 <- plyr::mapvalues(x = se6M1_3$nbs_5, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer6"
FeatureOverlay(se6M1_3, features = "layer6", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 7
ST.FeaturePlot(object = se6M1_3, features = "layer6", pt.size = 3, ncol = 2)

##seventh layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer6")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "6", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_6", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("6","nbs_6")
new.cluster.ids <- c(7,7)
se6M1_3$layer7 <- plyr::mapvalues(x = se6M1_3$nbs_6, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer7"
FeatureOverlay(se6M1_3, features = "layer7", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 8
ST.FeaturePlot(object = se6M1_3, features = "layer7", pt.size = 3, ncol = 2)

##8th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer7")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "7", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_7", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("7","nbs_7")
new.cluster.ids <- c(8,8)
se6M1_3$layer8 <- plyr::mapvalues(x = se6M1_3$nbs_7, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer8"
FeatureOverlay(se6M1_3, features = "layer8", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 9
ST.FeaturePlot(object = se6M1_3, features = "layer8", pt.size = 3, ncol = 2)

##9th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer8")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "8", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_8", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("8","nbs_8")
new.cluster.ids <- c(9,9)
se6M1_3$layer9 <- plyr::mapvalues(x = se6M1_3$nbs_8, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer9"
FeatureOverlay(se6M1_3, features = "layer9", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 10
ST.FeaturePlot(object = se6M1_3, features = "layer9", pt.size = 3, ncol = 2)

##10th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer9")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "9", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_9", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("9","nbs_9")
new.cluster.ids <- c(10,10)
se6M1_3$layer10 <- plyr::mapvalues(x = se6M1_3$nbs_9, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer10"
FeatureOverlay(se6M1_3, features = "layer10", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 11
ST.FeaturePlot(object = se6M1_3, features = "layer10", pt.size = 3, ncol = 2)

##11th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer10")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "10", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_10", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("10","nbs_10")
new.cluster.ids <- c(11,11)
se6M1_3$layer11 <- plyr::mapvalues(x = se6M1_3$nbs_10, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer11"
FeatureOverlay(se6M1_3, features = "layer11", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 12
ST.FeaturePlot(object = se6M1_3, features = "layer11", pt.size = 3, ncol = 2)

##12th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer11")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "11", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_11", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("11","nbs_11")
new.cluster.ids <- c(12,12)
se6M1_3$layer12 <- plyr::mapvalues(x = se6M1_3$nbs_11, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer12"
FeatureOverlay(se6M1_3, features = "layer12", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 13
ST.FeaturePlot(object = se6M1_3, features = "layer12", pt.size = 3, ncol = 2)

##13th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer12")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "12", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_12", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("12","nbs_12")
new.cluster.ids <- c(13,13)
se6M1_3$layer13 <- plyr::mapvalues(x = se6M1_3$nbs_12, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer13"
FeatureOverlay(se6M1_3, features = "layer13", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 14
ST.FeaturePlot(object = se6M1_3, features = "layer13", pt.size = 3, ncol = 2)

##14th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer13")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "13", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_13", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("13","nbs_13")
new.cluster.ids <- c(14,14)
se6M1_3$layer14 <- plyr::mapvalues(x = se6M1_3$nbs_13, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer14"
FeatureOverlay(se6M1_3, features = "layer14", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 15
ST.FeaturePlot(object = se6M1_3, features = "layer14", pt.size = 3, ncol = 2)

##15th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer14")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "14", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_14", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("14","nbs_14")
new.cluster.ids <- c(15,15)
se6M1_3$layer15 <- plyr::mapvalues(x = se6M1_3$nbs_14, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer15"
FeatureOverlay(se6M1_3, features = "layer15", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 16
ST.FeaturePlot(object = se6M1_3, features = "layer15", pt.size = 3, ncol = 2)

##16th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer15")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "15", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_15", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("15","nbs_15")
new.cluster.ids <- c(16,16)
se6M1_3$layer16 <- plyr::mapvalues(x = se6M1_3$nbs_15, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer16"
FeatureOverlay(se6M1_3, features = "layer16", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 17
ST.FeaturePlot(object = se6M1_3, features = "layer16", pt.size = 3, ncol = 2)

##17th layer
se6M1_3 <- SetIdent(se6M1_3, value = "layer16")
se6M1_3 <- RegionNeighbours(se6M1_3, id = "16", verbose = TRUE,keep.within.id = T)
FeatureOverlay(se6M1_3, features = "nbs_16", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)

current.cluster.ids <- c("16","nbs_16")
new.cluster.ids <- c(17,17)
se6M1_3$layer17 <- plyr::mapvalues(x = se6M1_3$nbs_16, from = current.cluster.ids, to = new.cluster.ids)
Idents(se6M1_3) <- "layer17"
FeatureOverlay(se6M1_3, features = "layer17", ncols = 2, sampleids = 1:2, cols = c("red", "lightgray"), pt.size = 2)
se6M1_3@meta.data[is.na(se6M1_3@meta.data)] <- 18
ST.FeaturePlot(object = se6M1_3, features = "layer17", pt.size = 3, ncol = 2)

##add column in meta.data to define mets, no mets level 
#extract the spots BC for each layer 
spots_layer0 <- row.names(se6M1_3@meta.data)[which(se6M1_3$seurat_clusters ==0)]
spots_layer1 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_0 =="nbs_0")]
spots_layer2 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_1 =="nbs_1")]
spots_layer3 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_2 =="nbs_2")]
spots_layer4 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_3 =="nbs_3")]
spots_layer5 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_4 =="nbs_4")]
spots_layer6 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_5 =="nbs_5")]
spots_layer7 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_6 =="nbs_6")]
spots_layer8 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_7 =="nbs_7")]
spots_layer9 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_8 =="nbs_8")]
spots_layer10 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_9 =="nbs_9")]
spots_layer11 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_10 =="nbs_10")]
spots_layer12 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_11 =="nbs_11")]
spots_layer13 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_12 =="nbs_12")]
spots_layer14 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_13 =="nbs_13")]
spots_layer15 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_14 =="nbs_14")]
spots_layer16 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_15 =="nbs_15")]
spots_layer17 <- row.names(se6M1_3@meta.data)[which(se6M1_3$nbs_16 =="nbs_16")]

se6M1_3$Mets_layer <- NA
for(i in 1:nrow(se6M1_3@meta.data)) {
  if (rownames(se6M1_3@meta.data)[i] %in% spots_layer0) {
    se6M1_3$Mets_layer[i] <- "layer0" 
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer1) {
    se6M1_3$Mets_layer[i] <- "layer1"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer2) {
    se6M1_3$Mets_layer[i] <- "layer2"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer3) {
    se6M1_3$Mets_layer[i] <- "layer3"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer4) {
    se6M1_3$Mets_layer[i] <- "layer4"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer5) {
    se6M1_3$Mets_layer[i] <- "layer5"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer6) {
    se6M1_3$Mets_layer[i] <- "layer6"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer7) {
    se6M1_3$Mets_layer[i] <- "layer7"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer8) {
    se6M1_3$Mets_layer[i] <- "layer8"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer9) {
    se6M1_3$Mets_layer[i] <- "layer9"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer10) {
    se6M1_3$Mets_layer[i] <- "layer10"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer11) {
    se6M1_3$Mets_layer[i] <- "layer11"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer12) {
    se6M1_3$Mets_layer[i] <- "layer12"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer13) {
    se6M1_3$Mets_layer[i] <- "layer13"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer14) {
    se6M1_3$Mets_layer[i] <- "layer14"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer15) {
    se6M1_3$Mets_layer[i] <- "layer15"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer16) {
    se6M1_3$Mets_layer[i] <- "layer16"
  }else if (rownames(se6M1_3@meta.data)[i] %in% spots_layer17) {
    se6M1_3$Mets_layer[i] <- "layer17"
  }else{
    se6M1_3$Mets_layer[i] <- "negative"
  }
}

#remove negatives 
Idents(se6M1_3) <- "Mets_layer"
se6M1_3 <- subset(se6M1_3, idents = c("layer0","layer1","layer2","layer3","layer4","layer5","layer6","layer7","layer8","layer9","layer10",
                                      "layer11","layer12","layer13","layer14","layer15","layer16","layer17"))

### plot the layers; layer 0 is the metastasis 
ST.FeaturePlot(object = se6M1_3, features = "Mets_layer", pt.size = 3, ncol = 2)

### save objects with layers 
saveRDS(se6M1_3,"./se6M1_3_layers.Rds")

###load R objects 
#layers of distance from metastasis are stored in seurat_object$Mets_layer (layer by layer, layer 0 = witin metastasis, layer 1 = first layer after metastasis, ...)
se6M1_3 <- readRDS("/media/Coco/MOSAIC LIVER/Experiments/Visium/Visium_Kristina/se6M1_3_layers.Rds")
SpatialDimPlot(se6M1_3$Mets_layer)
ST.FeaturePlot(object = se6M1_3, features = "Mets_layer", pt.size = 3, ncol = 2)
Idents(se6M1_3) <- "Mets_layer"
VlnPlot(se6M1_3, features = "App", pt.size=0, log=T)
VlnPlot(se6M1_3, features = "Saa1", pt.size=0, log=T)

#seurat_object$Mets_zone2 combines layers into three zones (mets = within mets, layer 1 = within 3 spots away from mets, rest is layer 2)
ST.FeaturePlot(object = se6M1_3, features = "Mets_zone2", pt.size = 3, ncol = 2)

#visualize genes of interest in an VlnPlot from differnet zones 
Idents(se6M1_3) <- "Mets_zone2"
VlnPlot(se6M1_3, features = "App")
VlnPlot(se6M1_3, features = "Saa1")


#set colors
colfunc <- colorRampPalette(c("#D12429", "#82C341"))
colors <- colfunc(18)

# Define an order of cluster identities
my_levels <- c("layer0", "layer1", "layer2", "layer3", "layer4",
               "layer5", "layer6", "layer7","layer8", "layer9", 
               "layer10", "layer11", "layer12", "layer13", "layer14",
               "layer15", "layer16", "layer17")

# Relevel object@ident
se6M1_3$Mets_layer <- factor(x = se6M1_3$Mets_layer, levels = my_levels)
Idents(se6M1_3) <- "Mets_layer"

FeatureOverlay(object = se6M1_3, features = "Mets_layer") +ggsave("/mets_layers_visium6M1.pdf", width = 8, height =8)
FeatureOverlay(se6M1_3, features = c("App", "Saa1", "Plxnb2", "Psen1"), 
               sampleids = 1,
               cols = c("lightgray", "mistyrose", "red", "darkred", "black"),
               pt.size = 3, 
               add.alpha = F,
               ncol = 4) +ggsave("visiumgenes.pdf", width = 32, height =8)

a <- VlnPlot(se6M1_3, features = "App", pt.size=0, log=T)+ theme(legend.position="none")
b <- VlnPlot(se6M1_3, features = "Saa1", pt.size=0, log=F)+ theme(legend.position="none")
c <- VlnPlot(se6M1_3, features = "Plxnb2", pt.size=0, log=F)+ theme(legend.position="none")
d <- VlnPlot(se6M1_3, features = "Psen1", pt.size=0, log=F)+ theme(legend.position="none")
ggarrange(a,b,c,d, ncol=2, nrow=2)+ggsave("mets_layers_visium.pdf", width = 12, height =20)

#check expression of genes across samples (se6M1 and se4M1 have mets, se6M3 doesn't)
se6M1$sample <- "AKPS_1"
se4M1$sample <- "AKPS_2"
se6M3$sample <- "no_mets"

se.merged <- MergeSTData(se6M1, c(se6M3,se4M1))
se.merged <- SCTransform(se.merged)
VlnPlot(se.merged, features = c("Saa1", "App", "Psen1", "Plxnb2"),ncol=2, group.by = "sample", pt.size=0, log=T)
VlnPlot(se.merged, features = c("Saa1"),ncol=2, group.by = "sample", pt.size=0, log=T)+ggsave("Saa1_visium.pdf", width = 5, height =5)

se.merged <- AddModuleScore(se.merged, features = "Saa1", name = "Saa1")
names(x = se.merged[[]])

Idents(se.merged) <- "sample"
AKPS_1 <- subset(se.merged, idents = c("AKPS_1"))
AKPS_2 <- subset(se.merged, idents = c("AKPS_2"))
no_mets <- subset(se.merged, idents = c("no_mets"))
wilcox.test(AKPS_1$Saa11, no_mets$Saa11, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(AKPS_2$Saa11, no_mets$Saa11, alternative = "two.sided") #p-value < 2.2e-16


