#!/bin/bash

##################################
## Analysis of Ariss et al., 2018
##################################

## Unpacking Counts Files

cd counts
tar -xzvf GSE115476_counts.tar.gz
cd ..

## Preparing Singularity Container
## ----------

## download container if it doesn't exist

if [ ! -f ../singularity/kumar_analysis_latest.sif ]
then
	cd ../singularity
	singularity pull shub://rpolicastro/kumar_analysis
fi

cd Ariss_2018

## initialize container and conda environment

singularity shell \
-eCB .. \
-H $PWD \
../singularity/kumar_analysis_latest.sif

source activate dropseq-analysis
