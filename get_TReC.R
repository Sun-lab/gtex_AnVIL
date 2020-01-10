
library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")

genes = readRDS("exon_by_genes_gencode.v26.GRCh38.rds")

filenames = list.files(path="~/research/data/_GTEx/v8/RNAseq", 
                       pattern=".bam$", full.names=TRUE)
filenames

fname   = filenames[1]
bamfile = BamFileList(fname, yieldSize=1000000)
bamfile

date()
se = summarizeOverlaps(features=genes, reads=bamfile, mode="Union",
             singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
date()

se
colData(se)
rowRanges(se)
str(metadata(rowRanges(se)))

as1 = assay(se)
as1 = as.data.frame(as1)
dim(as1)
head(as1)

write.table(as1, file = sprintf("%s_TReC.txt", fname), append = FALSE,
  quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,
  col.names = TRUE)

sessionInfo()
q(save="no")
