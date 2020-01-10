workflow collect_TReC {
    Int disk_size
    String memory_size
    Int preemptible

    call TReC {
        input:
        disk_size = disk_size, 
        memory_size = memory_size, 
        preemptible = preemptible
    }
    
    output {
        File trec_file = TReC.output_trec
    }
}

task TReC {
    File input_bam
    File bamIndex
    File gene_anno
    String input_bam_name
    String gene_anno_file_name

    Int disk_size
    String memory_size
    Int preemptible

    command {
        # Rscript --vanilla /get_TReC.R ${input_bam}
        Rscript -e "library(GenomicAlignments); library(GenomicFeatures); library(Rsamtools); genes = readRDS('${gene_anno_file_name}'); bamfile = BamFileList('${input_bam_name}', yieldSize=1000000); se = summarizeOverlaps(features=genes, reads=bamfile, singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE); as1 = as.data.frame(assay(se)); write.table(as1, file = paste0('${input_bam_name}', '_TReC.txt'));"
    }
    
    output {
        File output_trec = "${input_bam_name}.trec.txt"
    }
    
    runtime {
        continueOnReturnCode: false
        docker: "bioconductor/bioconductor_trecase:0.1"
        memory: memory_size
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible 
    }

}
