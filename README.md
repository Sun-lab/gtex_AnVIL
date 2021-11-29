# gtex_AnVIL

A workflow to process GTEx data at AnVIL cloud space. 

This repository include the workflow and docker file to process raw RNA-seq data from NHGRI AnVIL (Genomic Data Science Analysis, Visualization, and Informatics Lab-space, https://anvil.terra.bio/). 

## Docker image 
The docker file can be easily used in other computation environment and it is built on R/bioconductor, with addition of asSeq R package together with an R script "get_TReC_ASReC.R" that extracts both total read count (TReC) and allele-specific read count (ASReC) data. 

The docker configuration file [Dockerfile](Docker/Dockerfile) and the command to build the docker [build_docker.sh](Docker/build_docker.sh) saved in the folder Docker.
