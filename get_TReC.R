
library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")

genes = readRDS("exon_by_genes_gencode.v26.GRCh38.rds")

filenames = list.files(path="~/research/data/_GTEx/v8/RNAseq", 
                       pattern=".bam$", full.names=TRUE)
filenames

bamfiles  = BamFileList(filenames, yieldSize=1000000)
bamfiles

date()
se = summarizeOverlaps(features=genes, reads=bamfiles, mode="Union",
             singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
date()

se
colData(se)
rowRanges(se)
str(metadata(rowRanges(se)))

as1 = assay(se)

head(as1)
cts = rbind(cts, colSums(as1))

setwd(pipeline_dir)
write.table(as1, file = sprintf("%s/%s_gene_level_counts_filterIt_total.txt", proj_output, sam), append = FALSE,
  quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,
  col.names = TRUE)

row.names(cts) = c("nMappedReads", "nFragMappedToGene")
write.table(cts, file = sprintf("%s/%s_TReC_liberal_filterIt_total.txt", proj_output, sam), append = FALSE,
quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,
col.names = TRUE)

sessionInfo()
q(save="no")
