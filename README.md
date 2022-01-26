# gtex_AnVIL

A workflow to process GTEx bam files stored at NHGRI AnVIL (Genomic Data Science Analysis, Visualization, and Informatics Lab-space, https://anvil.terra.bio/) cloud space to obtain Total Read Count and Allele-specific read count per gene. 

This repository includes the workflow and docker configuration file . 

## Docker image 
The docker image is available as `sunway1999/bioconductor_trecase:0.1` from [DockerHub](https://hub.docker.com/). The image was built based on R/bioconductor (RELEASE_3_10), with addition of asSeq R package (asSeq_0.99.501.tar.gz) together with an R script [get_TReC_ASReC.R](Docker/get_TReC_ASReC.R) that extracts both total read/fragment count (TReC) and allele-specific read/fragment count (ASReC) data from bam files. 

The docker configuration file [Dockerfile](Docker/Dockerfile) and the command to build the docker image [build_docker.sh](Docker/build_docker.sh) are saved in the folder Docker. 

## Workflow 
The workflow [collect_TReC_ASReC.wdl](collect_TReC_ASReC.wdl) is written by Workflow Description Language (WDL), which specifies data processing workflows with a human-readable and writeable syntax. A WDL script chains several tasks to accomplish specific goals. Here our workflow is designed to run on google cloud environment through the AnVIL interface, and it has two tasks. One is to run the R script get_TReC_ASReC.R, and the other one is to move the output file to a specified location. 

The key part of the workflow is [get_TReC_ASReC.R](Docker/get_TReC_ASReC.R), which takes four arguments, and outputs total fragment count and two allele-specific counts (for two haplotypes) for each gene. The four arguments are 

1. The name of a bam file from which to collect count information.
2. The name of a gene annotation file. For example, file [exon_by_genes_gencode.v26.GRCh38.rds](_prepare_gene_anno/exon_by_genes_gencode.v26.GRCh38.rds), which is produced by a short script [step0_build_TxDb_local.R](_prepare_gene_anno/step0_build_TxDb_local.R).
3. Sample name, e.g., "GTEX-1117F-0426-SM-5EGHI".
4. The name of a text file for the list of heterozygous SNPs for one sample. The content of this file looks like the following, without header, with four columns: chromosome, position, and two alleles on haplotype 1 and 2, respectively. 
```
chr1	777	G	T
chr1	999	A	G
 ```
We collect the lists of heterozygous SNPs for each GTEx sample using phased genotype data from file: GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased, which is available from the `AnVIL_GTEx_V8_hg38` workspace in https://anvil.terra.bio/.


## Running this workflow in a other environments.
 
In a local environment, there is no need to move file and one can simply run the R script [get_TReC_ASReC.R](Docker/get_TReC_ASReC.R), for example, using the code like the following 
```
Rscript --vanilla get_TReC_ASReC.R sample1.bam exon_by_genes_gencode.v26.GRCh38.rds sample1 sample1_hetSNP.txt
```

The functions that we use from R bioconductor include `scanBamFlag` to filter bam files and `summarizeOverlaps` to obtain read/fragment counts. These functions should not vary in more recent R bioconductor releases. If it is desirable to run the analysis using the exact version of R environment we used, one can use a wdl script similar to [get_TReC_ASReC.R](Docker/get_TReC_ASReC.R), which uses the docker image: "sunway1999/bioconductor_trecase:0.1" from DockerHub. Outside AnVIL, an wdl script can be run by [Cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/), an execution engine supporting three types of platforms: local machine, a local cluster, or a cloud platform. For example, by a command  `java -jar ~/cromwell/cromwell-71.jar run my.wdl -i my.json` with input arguments specified by `my.json`.

We have inlcuded some sample data (bam files and SNP genotype data for two genes in 25 samples of 1000 Genome Project) in https://github.com/Sun-lab/asSeq_pipelines/tree/main/pipeline_1KGP/data, as well as a sample R code [test_scr.R](https://github.com/Sun-lab/asSeq_pipelines/blob/main/pipeline_1KGP/test_scr.R) to iterately all function ```get_TReC_ASReC.R``` to calculate counts from those bam files. 




