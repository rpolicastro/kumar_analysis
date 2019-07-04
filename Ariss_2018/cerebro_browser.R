#!/usr/bin/env Rscript

library("Seurat")
library("cerebroPrepare")

#########################################################
## Preparing Ariss et al., 2018 Data for Cerebro Browser
#########################################################

## Load data from seurat.

integrated.data <- readRDS("integrated_data.RDS")
Idents(integrated.data) <- "integrated_snn_res.1.2"

## Get most expressed genes.

cerebro.obj <- getMostExpressedGenes(
	integrated.data,
	column_sample="orig.ident",
	column_cluster="integrated_snn_res.1.2"
)

## Get marker genes.

cerebro.obj <- getMarkerGenes(
	cerebro.obj,
	organism="D. melanogaster",
	column_sample="orig.ident",
	column_cluster="integrated_snn_res.1.2",
	only_pos=FALSE,
	min_pct=0.25,
	thresh_use=log(1.5),
	return_thresh=0.05,
	test_use="wilcox"
)

## Export initial cerebro object.

exportFromSeurat(
	cerebro.obj,
	file="cerebro/ariss_2018_browser.crb",
	experiment_name="Ariss_2018",
	organism="D. melanogaster",
	column_sample="orig.ident",
	column_cluster="integrated_snn_res.1.2",
	column_nUMI="nCount_RNA",
	column_nGene="nFeature_RNA"
)

## Load back up cerebro file to change expression values.

cerebro.file <- readRDS("
