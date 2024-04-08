library(Seurat)
library(ggplot2)
library(MetBrewer)
Resolve_seurat_anno <- readRDS("../NAS/Kristina/Processed_sequencing_data/Sphere_seq_paper_data/
                               Highly_multiplexed_FISH_analysis/data_files_generated/Resolve_seurat_anno.rds")
Resolve_seurat_anno
  DimPlot(Resolve_seurat_anno)
DimPlot(Resolve_seurat_anno, group.by = "samples")
Idents(Resolve_seurat_anno) <- "samples"
no_mets <- subset(Resolve_seurat_anno, idents = "no_mets")
mets <- subset(Resolve_seurat_anno, idents = "mets")
FeaturePlot(no_mets, features = "Plxnb2", order=T)
Idents(Resolve_seurat_anno) <- "annotation"
hepatocytes <- subset(Resolve_seurat_anno, idents = c("Hepatocytes_CV", "Hepatocytes_PV"))
VlnPlot(hepatocytes, features="Plxnb2", cols=met.brewer("Archambault", 2), pt.size = 0)+
  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + ylim(0.01,4)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = "", y = "Semaphorin expression", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") + 
  ggsave("Plxnb2_PV_CV-8.337577e-164.pdf", width =5,  height=5)

dePlxnb2 <- FindMarkers(hepatocytes, ident.1= "Hepatocytes_PV", ident.2="Hepatocytes_CV", logfc.threshold=0)
View(dePlxnb2)

DotPlot(no_mets, features = "Plxnb2", group.by = "annotation")
VlnPlot(no_mets, features = "Plxnb2", pt.size=1, group.by = "annotation") +ylim(0.1, 4)+geom_boxplot()
VlnPlot(mets, features = "Plxnb2", pt.size=1, group.by = "annotation") +ylim(0.1, 4)

VlnPlot(no_mets, features = "Plxnb2", pt.size=1, log=T,  group.by = "vein")
VlnPlot(no_mets, features = "Plxnb2", pt.size=1, log=T,  group.by = "vein")
VlnPlot(no_mets, features = "Plxnb2", pt.size=0,  group.by = "vein") +ylim(0.1, 4)

VlnPlot(hepatocytes, features = "Plxnb2", pt.size=0, group.by = "Mets_distance") + ylim(0.1,5)
VlnPlot(mets, features = "Plxnb2", pt.size=0, group.by = "annotation") 
Idents(mets) <- "annotation"
hepatocytes_mets <- subset(mets, idents =c("Hepatocytes_CV", "Hepatocytes_PV"))
Idents(hepatocytes_mets) <- "Mets_distance"

VlnPlot(hepatocytes_mets, features="Plxnb2", cols=met.brewer("Archambault", 2), pt.size = 0)+geom_boxplot()+
  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + ylim(0.01,4)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = "", y = "Semaphorin expression", x="") + theme(legend.position="right") #+  
 # ggsave("Plxnb2_proximal_distal-8.337577e-164.pdf", width =5,  height=5)

Idents(mets) <- "Mets_distance"
demets <- FindMarkers(mets, ident.1= "proximal", ident.2="distal", logfc.threshold=0)
View(demets)

VlnPlot(mets, features = c("Sema4a", "Sema4c", "Sema4d", "Sema4g", "Plxnb2"), cols = c("yellow", "red", "green", "violet", "blue"), stack = TRUE, sort = TRUE,  flip = TRUE) + ylim(0.001,7) +
  theme(legend.position = "none") +ggsave("Semaexpression.pdf", width=6, height = 5)

Idents(mets) <- "annotation"
AverageExpression <- AverageExpression(mets, features = c("Sema4a", "Sema4c", "Sema4d", "Sema4g"))
library(pheatmap)
pheatmap(AverageExpression$originalexp, scale = 'none', cluster_rows = F, cluster_cols = F, border_color="white", 
         fontsize_number = 0.5, cellwidth =30, cellheight = 20, treeheight_row=0, treeheight_col=0,	
         angle_col = "45",color = colorRampPalette(rev(brewer.pal(n = 5, name =  "RdBu")))(100), filename = "heatmap.pdf")
