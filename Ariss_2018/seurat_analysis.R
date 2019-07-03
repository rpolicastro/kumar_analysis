#!/usr/bin/env Rscript

library("tidyverse")
library("Seurat")

############################################
## scRNA-seq Analysis of Ariss et al., 2018
############################################

## Preparing counts files.
## ----------

counts <- list.files("counts", pattern="*DMS") %>% file.path("counts",.) %>%
	setNames(basename(.) %>% substr(., 1, nchar(.)-9)) %>%
	split(., names(.)) %>% lapply(., unname) %>%
	map(
		~ read.delim(., sep="\t", header=T, stringsAsFactors=F) %>%
		column_to_rownames("GENE") %>%
		as.sparse
	)

## Creating seurat objects.
## ----------

seurat.obj <- map(names(counts), ~ CreateSeuratObject(counts=counts[[.]], project=.))

## Basic cell quality control.
## ----------

## Plotting number of genes and read counts.

pdf("./plots/QC_violin-plot.pdf")
map(
	names(count.files),
	~ VlnPlot(seurat.obj[[.]], features=c("nCount_RNA", "nFeature_RNA"), ncol=2)
)
dev.off()

## Filtering out cells by outlying number of reads and read counts.

seurat.obj <- map(
	seurat.obj,
	~ subset(., subset =
		nFeature_RNA > 250 &
		nFeature_RNA < 1000 &
		nCount_RNA > 500 &
		nCount_RNA < 5000
	)
)

## Sample Normalization.
## ----------

seurat.obj <- map(seurat.obj, ~ SCTransform(.))

## Determining Ideal Principle Components
## ----------

## Running initial PCA.

seurat.obj <- map(seurat.obj, ~ RunPCA(., npcs=50))

## Plotting elbow plots to pick ideal PCA number.

pdf("./plots/elbow-plots.pdf")
map(
	names(count.files),
	~ ElbowPlot(seurat.obj[[.]], ndims=50)
)
dev.off()

# Elbow is around 10-15 as determined by visual inspection of elbow plot.
# Will use 20 to be safe.

## Integrating data.
## ----------

## Finding integration anchors.

anchors <- FindIntegrationAnchors(seurat.obj, dims=1:20)

## Using anchors to integrate data.

integrated.data <- IntegrateData(anchors, dims=1:20)

## Setting integrated data as default assay.

DefaultAssay(integrated.data) <- "integrated"

## Scale data.

integrated.data <- ScaleData(integrated.data)

## Dimensionality reduction.
## ----------

## Running PCA.

integrated.data <- RunPCA(integrated.data, npcs=50)

## Plotting elbow plot.

pdf("./plots/integrated-elbow-plot.pdf")
ElbowPlot(integrated.data, ndims=50)
dev.off()

# Elbow is around 20-25 as determined by visual inspection of elbow plot.
# Will use 30 to be safe.

## UMAP for dimensionality reduction.

integrated.data <- RunUMAP(integrated.data, reduction="pca", dims=1:30)

## Clustering.
## ----------

## Finding neighbors.

integrated.data <- FindNeighbors(integrated.data, dims=1:30)

## Finding clusters.

integrated.data <- FindClusters(integrated.data, resolution=c(seq(0.2, 1.4, 0.2)))

## Plotting dimplots to find ideal resolution.

pdf("./plots/dimplot.pdf")
for (n in seq(0.2, 1.4, 0.2)) {
        Idents(integrated.data) <- paste0("integrated_snn_res.",n)
        p <- DimPlot(integrated.data, reduction="umap", label=TRUE)
        print(p)
}
dev.off()

# Ideal resolution seems to be ~1.2.

Idents(integrated.data) <- "integrated_snn_res.1.2"

pdf("./plots/dimplot_ideal-resolution.pdf")
DimPlot(integrated.data, reduction="umap", label=TRUE) + scale_color_viridis_d()
dev.off()

## Marker analysis.
## ----------

markers <- FindAllMarkers(
	integrated.data,
	logfc.threshold=log(1),
	test.use="wilcox",
	min.pct=0.2,
	only.pos=FALSE,
	return.thresh=0.05
)

## Export data for cell browser.
## ----------


