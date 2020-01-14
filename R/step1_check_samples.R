
library(data.table)

# ------------------------------------------------------------------------
# read in sample information
# ------------------------------------------------------------------------

sam_file = "../data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
sams = fread(sam_file)

dim(sams)
sams[1:2,]

table(sams$SMTS)
table(sams$SMTSD)
table(sams$SMAFRZE)

# ------------------------------------------------------------------------
# The samples with WGS data
# ------------------------------------------------------------------------

sams4WGS = sams$SAMPID[which(sams$SMAFRZE=="WGS")]
length(sams4WGS)
length(unique(sams4WGS))
sams4WGS[1:5]

sams4WGS = strsplit(sams4WGS, split="-")
table(sapply(sams4WGS, length))
sams4WGS = sapply(sams4WGS, function(v){paste(v[1:2], collapse='-')})
length(unique(sams4WGS))
sams4WGS[1:5]

# ------------------------------------------------------------------------
# Compare with the samples with het SNP data
# ------------------------------------------------------------------------

hetSNPs = list.files(path="~/research/data/_GTEx/v8/WGS_VCF/het_snps/byindnid_chr")
length(hetSNPs)
hetSNPs[1:5]

hetSNPs = gsub(".txt", "", hetSNPs, fixed=TRUE)
hetSNPs[1:5]

setequal(sams4WGS, hetSNPs)

# ------------------------------------------------------------------------
# The samples with RNA-seq data
# ------------------------------------------------------------------------

sams = sams[which(sams$SMAFRZE == "RNASEQ"),]
dim(sams)

sams4RNA = strsplit(sams$SAMPID, split="-")
sams4RNA = sapply(sams4RNA, function(v){paste(v[1:2], collapse='-')})

table(sams4RNA %in% sams4WGS)

sams = sams[which(sams4RNA %in% sams4WGS),]
dim(sams)

t1 = table(sams$SMTSD)
t1

# ------------------------------------------------------------------------
# compare with the summary table from GTEx portal
# ------------------------------------------------------------------------

t0 = read.csv("../data/GTEx_Portal.csv", as.is=TRUE)
dim(t0)
t0

table(names(t1) == t0$Tissue)
table(t1 - t0$X..RNASeq.and.Genotyped.samples)

# ------------------------------------------------------------------------
# check the tissues for one sample
# ------------------------------------------------------------------------

sams4RNA = strsplit(sams$SAMPID, split="-")
sams4RNA = sapply(sams4RNA, function(v){paste(v[1:2], collapse='-')})
sams$sams4RNA = sams4RNA

sams$SMTSD[which(sams4RNA == 'GTEX-1117F')]

# ------------------------------------------------------------------------
# compare with the sample list from GTEx files from google storage (gs)
# https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_GTEx_V8_hg38/data
# in the tab of TABLES, with 979 participants and 17382 samples. 
# ------------------------------------------------------------------------

sams_gs = fread("../data/sample.tsv")
dim(sams_gs)
sams_gs[1:2,]

table(sams$SAMPID %in% sams_gs$`entity:sample_id`)
sams_gs = sams_gs[match(sams$SAMPID, sams_gs$`entity:sample_id`),]
table(sams$SAMPID == sams_gs$`entity:sample_id`)

table(sams$SMTSD == sams_gs$tissue_site_detail)

# ------------------------------------------------------------------------
# generate list of bam files and bai files
# ------------------------------------------------------------------------

table(t1 > 300)
table(t1 > 250)
table(t1 > 200)
table(t1 > 100)

tissue2use = names(t1)[t1 > 300]
tissue_ids = gsub(" - ", "_", tissue2use, fixed = TRUE)
tissue_ids = gsub(" (", "_", tissue_ids, fixed = TRUE)
tissue_ids = gsub(")", "", tissue_ids,   fixed = TRUE)
tissue_ids = gsub(" ", "_", tissue_ids,  fixed = TRUE)

tissue_ids

for(k in 1:length(tissue2use)){
  tissue.k  = tissue2use[k]
  tissue.id = tissue_ids[k]
  w2kp      = which(sams_gs$tissue_site_detail==tissue.k)
  bam_files = sams_gs$bam_file[w2kp]
  bai_files = sams_gs$bam_index[w2kp]
  
  geno_dir = "gs://fc-secure-06f8c803-9dbf-48b5-a7d1-676b4867e9fb/byindnid_chr/"
  geno_files = paste0(geno_dir, sams_gs$participant[w2kp], ".txt")
  
  cat(bam_files, file=paste0("../list/", tissue.id, "_bams.list"), sep="\n")
  cat(bai_files, file=paste0("../list/", tissue.id, "_bais.list"), sep="\n")
  cat(geno_files, file=paste0("../list/", tissue.id, "_hetSNP.list"), sep="\n")
}

sessionInfo()
q(save="no")





