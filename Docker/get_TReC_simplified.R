
library(GenomicAlignments); library(GenomicFeatures); library(Rsamtools); genes = readRDS('${gene_anno_file_name}'); bamfile = BamFileList('${input_bam_name}', yieldSize=1000000); se = summarizeOverlaps(features=genes, reads=bamfile, singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE); as1 = as.data.frame(assay(se)); write.table(as1, file = paste0('${input_bam_name}', '_TReC.txt'));
