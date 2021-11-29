# gtex_AnVIL

A workflow to process GTEx bam files at AnVIL cloud space to obtain Total Read Count and Allele-specific read count per gene. 

This repository include the workflow and docker file to process raw RNA-seq data from NHGRI AnVIL (Genomic Data Science Analysis, Visualization, and Informatics Lab-space, https://anvil.terra.bio/). 

## Docker image 
The docker file can be easily used in other computation environment and it is built on R/bioconductor, with addition of asSeq R package together with an R script "get_TReC_ASReC.R" that extracts both total read count (TReC) and allele-specific read count (ASReC) data. 

The docker configuration file [Dockerfile](Docker/Dockerfile) and the command to build the docker image [build_docker.sh](Docker/build_docker.sh) are saved in the folder Docker. The content of the docker image is relatively simple. It inherits a docker instance of R bioconductor, installs a few R packages from bioconductor, and our R package asSeq. Then add an R script [get_TReC_ASReC.R](Docker/get_TReC_ASReC.R).

## Workflow 
The workflow [collect_TReC_ASReC.wdl](collect_TReC_ASReC.wdl) is written by Workflow Description Language (WDL), which specifies data processing workflows with a human-readable and writeable syntax. A WDL script chains several tasks to accomplish specific goals. Here our workflow is designed to run on google cloud environment through the AnVIL interface, and it has two tasks. One is to run the R script get_TReC_ASReC.R, and the other one is to move the output file to a specified location. 

The key part of the workflow is [get_TReC_ASReC.R](Docker/get_TReC_ASReC.R), which takes four arguments, and outputs total fragment count and two allele-specific counts (for two haplotypes) for each gene. the four arguments are 

1. the name of a bam file from which to collect count information
2. the name of a gene annotation file. For example, file [exon_by_genes_gencode.v26.GRCh38.rds](_prepare_gene_anno/exon_by_genes_gencode.v26.GRCh38.rds), which is produced by a short script [step0_build_TxDb_local.R](_prepare_gene_anno/step0_build_TxDb_local.R).
3. sam_name, e.g., something like "GTEX-1117F-0426-SM-5EGHI".
4. the name of a text file for the list of heterozygous SNPs, e.g., "GTEX-1117F.txt".
 

### get_TReC_ASReC.R

## Running this workflow in a local environment. 
Since this workflow 



