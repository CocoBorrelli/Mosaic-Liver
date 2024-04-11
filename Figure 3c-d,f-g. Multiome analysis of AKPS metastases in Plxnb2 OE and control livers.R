#This code performs the analysis of multiomic dataset


######CREATE UNIFIED SET OF PEAKS BEFORE MERGING OBJECT######
#To create a unified set of peaks we can use functions from the GenomicRanges package. 
#The reduce function from GenomicRanges will merge all intersecting peaks.

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("GenomicRanges")

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(magrittr)
library(dplyr)

gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(20, 70, 300), end = c(120, 200, 400)))
gr

# read in peak sets
peaks.OE <- read.table(
  file = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/OE_outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.NT <- read.table(
  file = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/NT_outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)

peaks.OE.2 <- read.table(
  file = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/OE2/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.NT.2 <- read.table(
  file = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/NT2/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
peaks.OE <- makeGRangesFromDataFrame(peaks.OE)
peaks.NT <- makeGRangesFromDataFrame(peaks.NT)
peaks.OE.2 <- makeGRangesFromDataFrame(peaks.OE.2)
peaks.NT.2 <- makeGRangesFromDataFrame(peaks.NT.2)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(peaks.OE, peaks.NT, peaks.OE.2, peaks.NT.2))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
md.OE <- read.table(
  file = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/OE_outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.NT <- read.table(
  file = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/NT_outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.OE.2 <- read.table(
  file = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/OE2/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.NT.2 <- read.table(
  file = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/NT2/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
#md.OE <- md.OE[md.OE$passed_filters > 500, ]
#md.NT <- md.NT[md.NT$passed_filters > 500, ]


######CREATE SEPARATE OBJECT WITH RNA, ATAC AND CALLED PEAKS######
# create fragment objects
frags.OE <- CreateFragmentObject(
  path = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/OE_outs/atac_fragments.tsv.gz",
  cells = rownames(md.OE),
  validate.fragments = F
)

frags.NT <- CreateFragmentObject(
  path = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/NT_outs/atac_fragments.tsv.gz",
  cells = rownames(md.NT),
  validate.fragments = F
)

frags.OE.2 <- CreateFragmentObject(
  path = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/OE2/outs/atac_fragments.tsv.gz",
  cells = rownames(md.OE.2),
  validate.fragments = F
)

frags.NT.2 <- CreateFragmentObject(
  path = "~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/NT2/outs/atac_fragments.tsv.gz",
  cells = rownames(md.NT.2),
  validate.fragments = F
)

#Quantify peaks in each dataset
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

OE.counts <- FeatureMatrix(
  fragments = frags.OE,
  features = combined.peaks,
  cells = rownames(md.OE)
)

NT.counts <- FeatureMatrix(
  fragments = frags.NT,
  features = combined.peaks,
  cells = rownames(md.NT)
)

OE.counts.2 <- FeatureMatrix(
  fragments = frags.OE.2,
  features = combined.peaks,
  cells = rownames(md.OE.2)
)

NT.counts.2 <- FeatureMatrix(
  fragments = frags.NT.2,
  features = combined.peaks,
  cells = rownames(md.NT.2)
)

#Create the objects
# get gene annotations 
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA data
counts_OE <- Read10X_h5("/media/Coco/MOSAIC LIVER/Experiments/nuclei/OE_outs/filtered_feature_bc_matrix.h5")
OE <- CreateSeuratObject(
  counts = counts_OE$`Gene Expression`,
  assay = "RNA", project = "OE" 
)

# create ATAC assay and add it to the object
OE[["ATAC"]] <- CreateChromatinAssay(
  counts = OE.counts[, colnames(counts_OE$`Gene Expression`)],
  sep = c(":", "-"),
  fragments = frags.OE,
  annotation = annotation
)
OE

#call peaks using MACS2
DefaultAssay(OE) <- "ATAC"
peaks_OE <- CallPeaks(OE, macs2.path="/home/borrellc/miniconda3/envs/macs2-new/bin/macs2", outdir = tempdir("~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei"),
                      fragment.tempdir = tempdir("~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei"))

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks_OE <- keepStandardChromosomes(peaks_OE, pruning.mode = "coarse")
peaks_OE <- subsetByOverlaps(x = peaks_OE, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(OE),
  features = peaks_OE,
  cells = colnames(OE)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
OE[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = frags.OE,
  annotation = annotation
)

#now for NT
counts_NT <- Read10X_h5("/media/Coco/MOSAIC LIVER/Experiments/nuclei/NT_outs/filtered_feature_bc_matrix.h5")

NT <- CreateSeuratObject(
  counts = counts_NT$`Gene Expression`,
  assay = "RNA", project = "NT" 
)

# create ATAC assay and add it to the object
NT[["ATAC"]] <- CreateChromatinAssay(
  counts = NT.counts[, colnames(counts_NT$`Gene Expression`)],
  sep = c(":", "-"),
  fragments = frags.NT,
  annotation = annotation
)
NT


#call peaks using MACS2
DefaultAssay(NT) <- "ATAC"
peaks_NT <- CallPeaks(NT, macs2.path="/home/borrellc/miniconda3/envs/macs2-new/bin/macs2", outdir = tempdir("~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei"),
                      fragment.tempdir = tempdir("~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei"))

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks_NT <- keepStandardChromosomes(peaks_NT, pruning.mode = "coarse")
peaks_NT <- subsetByOverlaps(x = peaks_NT, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts_NT <- FeatureMatrix(
  fragments = Fragments(NT),
  features = peaks_NT,
  cells = colnames(NT)
)

# create a NTw assay using the MACS2 peak set and add it to the Seurat object
NT[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts_NT,
  fragments = frags.NT,
  annotation = annotation
)


# create a Seurat object containing the RNA data
counts_OE.2 <- Read10X_h5("/media/Coco/MOSAIC LIVER/Experiments/nuclei/OE2/outs/filtered_feature_bc_matrix.h5")
OE.2 <- CreateSeuratObject(
  counts = counts_OE.2$`Gene Expression`,
  assay = "RNA", project = "OE.2" 
)

# create ATAC assay and add it to the object
OE.2[["ATAC"]] <- CreateChromatinAssay(
  counts = OE.counts.2[, colnames(counts_OE.2$`Gene Expression`)],
  sep = c(":", "-"),
  fragments = frags.OE.2,
  annotation = annotation
)
OE.2

#call peaks using MACS2
DefaultAssay(OE.2) <- "ATAC"
peaks_OE.2 <- CallPeaks(OE.2, macs2.path="/home/borrellc/miniconda3/envs/macs2-new/bin/macs2", outdir = tempdir("~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei"),
                      fragment.tempdir = tempdir("~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei"))

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks_OE.2 <- keepStandardChromosomes(peaks_OE.2, pruning.mode = "coarse")
peaks_OE.2 <- subsetByOverlaps(x = peaks_OE.2, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts.OE2 <- FeatureMatrix(
  fragments = Fragments(OE.2),
  features = peaks_OE.2,
  cells = colnames(OE.2)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
OE.2[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts.OE2,
  fragments = frags.OE.2,
  annotation = annotation
)

#now for NT
counts_NT.2 <- Read10X_h5("/media/Coco/MOSAIC LIVER/Experiments/nuclei/NT2/outs/filtered_feature_bc_matrix.h5")

NT.2 <- CreateSeuratObject(
  counts = counts_NT.2$`Gene Expression`,
  assay = "RNA", project = "NT.2" 
)

# create ATAC assay and add it to the object
NT.2[["ATAC"]] <- CreateChromatinAssay(
  counts = NT.counts.2[, colnames(counts_NT.2$`Gene Expression`)],
  sep = c(":", "-"),
  fragments = frags.NT.2,
  annotation = annotation
)
NT.2


#call peaks using MACS2
DefaultAssay(NT.2) <- "ATAC"
peaks_NT.2 <- CallPeaks(NT.2, macs2.path="/home/borrellc/miniconda3/envs/macs2-new/bin/macs2", outdir = tempdir("~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei"),
                      fragment.tempdir = tempdir("~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei"))

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks_NT.2 <- keepStandardChromosomes(peaks_NT.2, pruning.mode = "coarse")
peaks_NT.2 <- subsetByOverlaps(x = peaks_NT.2, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts_NT.2 <- FeatureMatrix(
  fragments = Fragments(NT.2),
  features = peaks_NT.2,
  cells = colnames(NT.2)
)

# create a NTw assay using the MACS2 peak set and add it to the Seurat object
NT.2[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts_NT.2,
  fragments = frags.NT.2,
  annotation = annotation
)

######MERGE######
# add information to identify dataset of origin
OE$sample <- 'sgPlxnb2_1'
NT$sample <- 'sgNT_1'
OE.2$sample <- 'sgPlxnb2_2'
NT.2$sample <- 'sgNT_2'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(OE, c(NT, OE.2, NT.2), add.cell.ids = c("OE", "NT","OE.2", "NT.2" ))

combined[["ATAC"]]

#Quality control
DefaultAssay(combined) <- "ATAC"

combined <- NucleosomeSignal(combined)
combined <- TSSEnrichment(combined)

VlnPlot(
  object = combined,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
combined <- subset(
  x = combined,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 900 &
    nCount_RNA > 900 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
combined

saveRDS(combined, file = "combined.rds")
save.image("/media/Coco/MOSAIC LIVER/Experiments/nuclei/4samples.RData")
