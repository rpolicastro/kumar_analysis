#!/usr/bin/env Rscript

library("tidyverse")
library("Seurat")

############################################
## scRNA-seq Analysis of Ariss et al., 2018
############################################

## Preparing counts files.
## ----------

count.files <- list.files("counts", pattern="*DMS") %>% file.path("counts",.) %>%
	setNames(basename(.) %>% substr(., 1, nchar(.)-9)) %>%
	split(., names(.)) %>% lapply(., unname) %>%
	map(
		~ read.delim(., sep="\t", header=T, stringsAsFactors=F) %>%
		column_to_rownames("GENE") %>%
		as.sparse
	)

## Creating seurat objects.
## ----------

## Create object.

seurat.obj <- map(counts, ~CreateSeuratObject(counts=.,))

## Basic cell quality control.

# Plotting number of genes and read counts.
pdf("QC_violin-plot.pdf")
map(
	names(count.files),
	~ VlnPlot(seurat.obj[[.]], features=c("nCount_RNA", "nFeature_RNA"), ncol=2)
)
dev.off()

# Filtering out cells by outlying number of reads and read counts.
seurat.obj <- map(
	seurat.obj,
	~ subset(., subset =
		nFeature_RNA > 250 &
		nFeature_RNA < 1000 &
		nCount_RNA > 500 &
		nCount_RNA < 5000
	)
)
