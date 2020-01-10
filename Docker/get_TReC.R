
library(GenomicFeatures)
library(GenomicAlignments)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
  fname = args[1]
}

genes = readRDS("exon_by_genes_gencode.v26.GRCh38.rds")

# filenames = list.files(path="~/research/data/_GTEx/v8/RNAseq", 
#                        pattern=".bam$", full.names=TRUE)
# filenames
# fname   = filenames[1]

bamfile = BamFileList(fname, yieldSize=1000000)
bamfile

date()
se = summarizeOverlaps(features=genes, reads=bamfile, mode="Union",
             singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
date()

as1 = as.data.frame(assay(se))
dim(as1)
head(as1)

write.table(as1, file = sprintf("%s_TReC.txt", fname), 
  quote = FALSE, sep = "\t", eol = "\n")

sessionInfo()
q(save="no")
