# Prep files for gtex_anvil

rm(list=ls()) # clear workspace
getwd() # assuming working directory is ./gtex_AnVIL/R
setwd("../")
repo_dir 	= getwd(); repo_dir
data_dir	= file.path(repo_dir,"data")
tissue		= "BrainFrontalCortex"

# Local functions/Libraries
library(data.table)
smart_merge = function(x,y,mess=NULL,...){
	if( !is.null(mess) ){
		intersect_vars = paste(intersect(names(x),names(y)),collapse=", ")
		cat(paste0("Merging dataframes on variables = { ",intersect_vars," }\n"))
	}
	
	merge(x,y,by=intersect(names(x),names(y)),...)
}
smart_table = function(...){
	table(...,useNA='ifany')
}

# Combine data for one tissue
anno_fn 	= file.path(data_dir,"GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
pheno_fn 	= file.path(data_dir,"GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
anno 			= data.table::fread(anno_fn,header = TRUE,sep = "\t",data.table = FALSE)
pheno 		= data.table::fread(pheno_fn,header = TRUE,sep = "\t",data.table = FALSE)

dim(pheno); pheno[1:2,]
names(pheno)[names(pheno) == "AGE"] = "AGE_RANGE"
names(pheno)[names(pheno) == "SEX"] = "SEX_NUM"
dim(anno); anno[1:2,]
sort(names(anno))

smart_table(anno$SMTS) # tissue type
smart_table(anno$SMTSD[which(anno$SMTS == "Brain")]) # more specific sub-tissue type
smart_table(anno$SMAFRZE)

## Combine anno vars with sex/age_range
anno$SUBJID = sapply(anno$SAMPID,function(xx) 
	paste(strsplit(xx,"-")[[1]][1:2],collapse="-"),
	USE.NAMES = FALSE)
dat = smart_merge(pheno,anno)
dim(dat); dat[1:2,]

## Subset sub-tissue
if( tissue == "BrainFrontalCortex" ){
	dat = dat[which(dat$SMTSD == "Brain - Frontal Cortex (BA9)"
		& dat$SMAFRZE == "RNASEQ"),]
} else {
	stop("That tissue type hasn't been coded for yet!")
}
dim(dat); length(unique(dat$SUBJID)) 
dat[1:2,]

## Get google cloud link for bam/bai
samp_fn = file.path(data_dir,"sample.tsv")
samp = data.table::fread(samp_fn,header = TRUE,sep = "\t",data.table = FALSE)
dim(samp); samp[1:2,]
names(samp)[names(samp) == "entity:sample_id"] = "SAMPID"
names(samp)[names(samp) == "participant"] = "SUBJID"
# smart_table(samp$tissue_site_detail)
samp = samp[which(samp$SAMPID %in% dat$SAMPID),]
dim(samp); samp[1:2,]
dat = smart_merge(dat,samp)
dim(dat); dat[1:2,]

## Append other columns from participant.tsv
part_fn = file.path(data_dir,"participant.tsv")
part = data.table::fread(part_fn,header = TRUE,sep = "\t",data.table = FALSE)
dim(part); part[1:2,]
names(part)[names(part) == "entity:participant_id"] = "SUBJID"
part = part[which(part$SUBJID %in% dat$SUBJID),]
dat = smart_merge(dat,part)
dim(dat); dat[1:2,]

## Subset subjects with both rnaseq and genotype
sort(names(dat))
smart_table(dat$SEX_NUM,dat$sex)
smart_table(dat[,c("has_genotype","has_rnaseq")])
dat = dat[which(dat$has_genotype & dat$has_rnaseq),]
dim(dat)



q("no")

##