#!/usr/bin/env Rscript

library("tidyverse")
library("Seurat")

############################################
## scRNA-seq Analysis of Ariss et al., 2018
############################################

## Preparing counts files.

counts <- list.files("counts", pattern="*DMS") %>%
	file.path("counts",.) %>%
	map(
		~ read.delim(., sep="\t", header=T, stringsAsFactors=F) %>%
		column_to_rownames("GENE") %>%
		as.matrix
	)
