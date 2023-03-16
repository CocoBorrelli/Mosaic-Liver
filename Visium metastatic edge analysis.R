library(STutility)

load("Met1_seuratObj_V2.RData")
load("Met1rep_seuratObj_V2.RData")
FeatureOverlay(Met1, features = "seurat_predicted.id", sampleids = 1, ncols.samples = 1, dark.theme = T)


#STUtility
infoTable <- data.frame(samples=NA, spotfiles=NA, imgs=NA, json=NA, condition=NA)
infoTable <- rbind(infoTable, data.frame(samples="B1_FRB33CB70_Met/outs/filtered_feature_bc_matrix.h5", 
                                         spotfiles="B1_FRB33CB70_Met/outs/spatial/tissue_positions_list.csv", 
                                         imgs="B1_FRB33CB70_Met/outs/spatial/tissue_hires_image.png", 
                                         json="B1_FRB33CB70_Met/outs/spatial/scalefactors_json.json", condition="Met1"))
dim(infoTable)
infoTable <- infoTable[-1,]

STU <- InputFromTable(infotable = infoTable, 
                      min.gene.count = 0, 
                      min.gene.spots = 0,
                      min.spot.count = 0,
                      platform =  "Visium")

Met1.rep.STU <- STU
Met1.rep.temp <- Met1.rep
rownames(Met1.rep.temp@meta.data) <- paste(rownames(Met1.rep.temp@meta.data), "_1", sep = "")

#rownames(Met1.rep.STU@meta.data) <- sapply(rownames(Met1.rep.STU@meta.data), function(x) unlist(strsplit(x, split = "_"))[1])
Met1.rep.STU <- AddMetaData(Met1.rep.STU, Met1.rep.temp@meta.data[rownames(Met1.rep.STU@meta.data), "seurat_predicted.id"], "seurat_predicted.id")
Met1.rep.STU <- SCTransform(Met1.rep.STU)

#Met1.rep.STU@assays$RNA@counts <- Met1.rep@assays$Spatial@counts
#Met1.rep.STU@assays$RNA@data <- Met1.rep@assays$Spatial@data
#Met1.rep.STU@assays$RNA@scale.data <- Met1.rep@assays$Spatial@scale.data
Met1.rep.STU <- LoadImages(Met1.rep.STU, time.resolve = F, verbose = T)

Met1.rep.STU <- SetIdent(Met1.rep.STU, value = "seurat_predicted.id")
FeatureOverlay(Met1.rep.STU, features = "seurat_predicted.id", pt.size = 2, dark.theme = F)

Met1.rep.STU <- RegionNeighbours(Met1.rep.STU, id = "Hepatocytes", verbose = TRUE)
FeatureOverlay(Met1.rep.STU, features = "nbs_Hepatocytes", cols = c("red", "lightgray"), pt.size = 2, dark.theme = F)
Hepa.nb.tumor <- rep(NA, nrow(Met1.rep.STU@meta.data))
Hepa.nb.tumor[which(Met1.rep.STU@meta.data$seurat_predicted.id == "Tumor_Carcinoma" & Met1.rep.STU@meta.data$nbs_Hepatocytes == "nbs_Hepatocytes")] <- "Tumor_nb"
Hepa.nb.tumor[which(Met1.rep.STU@meta.data$seurat_predicted.id == "Hepatocytes" & Met1.rep.STU@meta.data$nbs_Hepatocytes == "Hepatocytes")] <- "Hepatocytes"
Met1.rep.STU <- AddMetaData(Met1.rep.STU, Hepa.nb.tumor, "nbs_Hepatocytes_Tumors")
FeatureOverlay(Met1.rep.STU, features = "nbs_Hepatocytes_Tumors", cols = c("red", "yellow"), pt.size = 2, dark.theme = F)
Hepa.tumor.iner.edge <- rep(NA, nrow(Met1.rep.STU@meta.data))
for (i in 1:nrow(Met1.rep.STU@meta.data)){
  if (Met1.rep.STU@meta.data[i,"seurat_predicted.id"] == "Hepatocytes" & is.na(Met1.rep.STU@meta.data[i,"nbs_Hepatocytes_Tumors"])) {Hepa.tumor.iner.edge[i] <- "Hepatocyte_central"}
  else if (Met1.rep.STU@meta.data[i,"seurat_predicted.id"] == "Hepatocytes" & Met1.rep.STU@meta.data[i,"nbs_Hepatocytes_Tumors"] == "Hepatocytes") {Hepa.tumor.iner.edge[i] <- "Hepatocyte_TumorNeighbour"}
  if (Met1.rep.STU@meta.data[i,"seurat_predicted.id"] == "Tumor_Carcinoma" & is.na(Met1.rep.STU@meta.data[i,"nbs_Hepatocytes_Tumors"])) {Hepa.tumor.iner.edge[i] <- "Tumor_central"}
  else if (Met1.rep.STU@meta.data[i,"seurat_predicted.id"] == "Tumor_Carcinoma" & Met1.rep.STU@meta.data[i,"nbs_Hepatocytes_Tumors"] == "Tumor_nb") {Hepa.tumor.iner.edge[i] <- "Tumor_HepatoNeighbour"}
}
Met1.rep.STU <- AddMetaData(Met1.rep.STU, Hepa.tumor.iner.edge, "Hepa_tumor_iner_edge")

FeatureOverlay(Met1.rep.STU, features = "Hepa_tumor_iner_edge", cols = c("green", "yellow", "red", "blue","grey"), pt.size = 2, dark.theme = F)


library(magrittr)
library(dplyr)

Met1.rep.STU <- SetIdent(Met1.rep.STU, value = "Hepa_tumor_iner_edge")
#Met1.rep.STU@assays$RNA@scale.data <- Met1.rep@assays$SCT@scale.data
#Met1.rep.STU@assays$RNA@scale.data <- as.matrix(Met1.rep.STU@assays$RNA@scale.data)
hepato.edge.central.markers.Met1.rep <- FindMarkers(Met1.rep.STU, ident.1 = "Hepatocyte_TumorNeighbour", ident.2 = "Hepatocyte_central")
Tumor.edge.central.markers.Met1.rep <- FindMarkers(Met1.rep.STU, ident.1 = "Tumor_HepatoNeighbour", ident.2 = "Tumor_central")
save(hepato.edge.central.markers.Met1.rep,file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/DEGs_HepatoNeighbouringTumor_VS_HepatoCentral_Met1rep_V2.RData")
save(Tumor.edge.central.markers.Met1.rep,file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/DEGs_TumorNeighbouringHepato_VS_TumorCentral_Met1rep_V2.RData")
#hepato.edge.vs.central.DEGs <- hepato.edge.central.markers$gene
#Tumor.edge.vs.central.DEGs <- Tumor.edge.central.markers$gene
#save(hepato.edge.vs.central.DEGs,file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/DEGs_HepatoNeighbouringTumor_VS_HepatoCentral_genes.RData")
#save(Tumor.edge.vs.central.DEGs,file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/DEGs_TumorNeighbouringHepato_VS_TumorCentral_genes.RData")

hepato.edge.central.markers.Met1.rep$gene <- rownames(hepato.edge.central.markers.Met1.rep)
hepato.edge.central.Met1.rep.STU.subset <- SubsetSTData(Met1.rep.STU, expression = Hepa_tumor_iner_edge %in% c("Hepatocyte_TumorNeighbour", "Hepatocyte_central"))
hepato.edge.central.Met1.rep.sorted.marks <- hepato.edge.central.markers.Met1.rep %>% top_n(n = 40, wt = abs(avg_log2FC))
hepato.edge.central.Met1.rep.sorted.marks <- hepato.edge.central.Met1.rep.sorted.marks[order(hepato.edge.central.Met1.rep.sorted.marks$avg_log2FC, decreasing = T), ]
DoHeatmap(hepato.edge.central.Met1.rep.STU.subset, features = hepato.edge.central.Met1.rep.sorted.marks$gene, group.colors = c("red", "lightgray"), disp.min = -2, disp.max = 2) 

Tumor.edge.central.markers.Met1.rep$gene <- rownames(Tumor.edge.central.markers.Met1.rep)
Tumor.edge.central.Met1.rep.STU.subset <- SubsetSTData(Met1.rep.STU, expression = Hepa_tumor_iner_edge %in% c("Tumor_HepatoNeighbour", "Tumor_central"))
Tumor.edge.central.Met1.rep.sorted.marks <- Tumor.edge.central.markers.Met1.rep %>% top_n(n = 40, wt = abs(avg_log2FC))
Tumor.edge.central.Met1.rep.sorted.marks <- Tumor.edge.central.Met1.rep.sorted.marks[order(Tumor.edge.central.Met1.rep.sorted.marks$avg_log2FC, decreasing = T), ]
DoHeatmap(Tumor.edge.central.Met1.rep.STU.subset, features = Tumor.edge.central.Met1.rep.sorted.marks$gene, group.colors = c("red", "lightgray"), disp.min = -2, disp.max = 2) 


#Nichenet on neighbouring cells
library(nichenetr)
library(sparseMatrixStats)
Met1.rep.nichenet <- Met1.rep.STU
Met1.rep.nichenet@meta.data %>% head()
Idents(Met1.rep.nichenet) <- "Hepa_tumor_iner_edge"
Met1.rep.nichenet <- RunPCA(Met1.rep.nichenet, verbose = FALSE)  %>% RunUMAP(dims = 1:30)
DimPlot(Met1.rep.nichenet, reduction = "umap")

ligand_target_matrix = readRDS("/Users/ati/Documents/Projects/Visium/colorectal_cancer/nichenet/ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
dim(ligand_target_matrix)

expression = t(as.matrix(Met1.rep.nichenet@assays$RNA@counts))
sample_info = Met1.rep.nichenet@meta.data

# we wanna check between Hepatocytes: sender, and Tumor: reciever
Hepato_ids = rownames(sample_info[which(sample_info$Hepa_tumor_iner_edge == "Hepatocyte_TumorNeighbour"),])
Tumor_ids = rownames(sample_info[which(sample_info$Hepa_tumor_iner_edge == "Tumor_HepatoNeighbour"),])

expressed_genes_sender = colnames(expression[Hepato_ids,])[which(colSums2(expression[Hepato_ids,]) > length(Hepato_ids)*0.1)]
expressed_genes_receiver = colnames(expression[Tumor_ids,])[which(colSums2(expression[Tumor_ids,]) > length(Tumor_ids)*0.1)]
length(expressed_genes_sender)
length(expressed_genes_receiver)

#Step 2: Define the gene set of interest and a background of genes
#genes of which the expression is possibly affected due to communication with other cells.
#Now are you interested in hepato genes or Tumor genes? Upregulated ones or downregulated ones or all
geneset_oi = rownames(Tumor.edge.central.markers)#[which(hepato.edge.central.markers$avg_log2FC > 0),])
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
length(geneset_oi)
length(background_expressed_genes)
length(intersect(geneset_oi, background_expressed_genes))

#Step 3: Define a set of potential ligands
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
save(lr_network_expressed, file= "/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/Hepato_Tumor_all_expressed_ligandReceptorNetwork_NoFiltering.RData")
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

#Now perform the ligand activity analysis: in this analysis, we will calculate the ligand activity of each ligand, or
# in other words, we will assess how well each Hepato-ligand can predict the Tumor gene set compared to the background of expressed genes 
#(predict whether a gene belongs to the Tumor program or not).
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson)
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
save(ligand_activities, file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/Hepato_Tumor_active_ligand_activities_weighted_all_V2.RData")
save(active_ligand_target_links_df, file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/Hepato_Tumor_active_ligand_target_links_df_V2.RData")
save(vis_ligand_target, file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/Hepato_Tumor_vis_ligand_target_V2.RData")

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Hepatocytes-ligands","Tumor_Tumor genes in Metastatic cells", color = "red",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network +theme(text = element_text(size=15, face = "bold"),axis.text.x = element_text(size=15,face="bold"))




##Nichenet on opposite interaction: Tumor as sender and Hepato as Reciever
# we wanna check between Hepatocytes: reciever, and Tumor: sender
Hepato_ids = rownames(sample_info[which(sample_info$Hepa_tumor_iner_edge == "Hepatocyte_TumorNeighbour"),])
Tumor_ids = rownames(sample_info[which(sample_info$Hepa_tumor_iner_edge == "Tumor_HepatoNeighbour"),])

expressed_genes_sender = colnames(expression[Tumor_ids,])[which(colSums2(expression[Tumor_ids,]) > length(Tumor_ids)*0.1)]
expressed_genes_receiver = colnames(expression[Hepato_ids,])[which(colSums2(expression[Hepato_ids,]) > length(Hepato_ids)*0.1)]
length(expressed_genes_sender)
length(expressed_genes_receiver)

#Step 2: Define the gene set of interest and a background of genes
#genes of which the expression is possibly affected due to communication with other cells.
#Now are you interested in hepato genes or Tumor genes? Upregulated ones or downregulated ones or all
geneset_oi = rownames(hepato.edge.central.markers)#[which(hepato.edge.central.markers$avg_log2FC > 0),])
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
length(geneset_oi)
length(background_expressed_genes)
length(intersect(geneset_oi, background_expressed_genes))

#Step 3: Define a set of potential ligands
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
save(lr_network_expressed, file= "/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/Tumor_Hepato_all_expressed_ligandReceptorNetwork.RData")
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

#Now perform the ligand activity analysis: in this analysis, we will calculate the ligand activity of each ligand, or
# in other words, we will assess how well each Hepato-ligand can predict the Tumor gene set compared to the background of expressed genes 
#(predict whether a gene belongs to the Tumor program or not).
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson)
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
save(ligand_activities, file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/Tumor_Hepato_Tumor_active_ligand_activities_weighted_all_V2.RData")
save(active_ligand_target_links_df, file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/Tumor_Hepato_active_ligand_target_links_df_V2.RData")
save(vis_ligand_target, file="/Users/ati/Documents/Projects/Visium/colorectal_cancer/coco/Tumor_Hepato_vis_ligand_target_V2.RData")

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Tumor-ligands","Hepatocytes genes in Metastatic cells", color = "red",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network +theme(text = element_text(size=15, face = "bold"),axis.text.x = element_text(size=8,face="bold"))
