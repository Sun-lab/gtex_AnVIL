
workflow collect_TReC_ASReC {
    Int disk_size
    Int n_CPU
    String memory_size
    String tissue
    Int preemptible

    Array[File] bam_files
    Array[File] bam_indices
    Array[File] het_snp_files
    File gene_anno

    Array[Int] indices = range(length(bam_files))
    
    scatter (i in indices) {

        call TReC_ASReC {
            input:
            disk_size   = disk_size, 
            n_CPU       = n_CPU, 
            memory_size = memory_size, 
            preemptible = preemptible,
            file_bam    = bam_files[i],
            file_bai    = bam_indices[i],
            file_hetSNP = het_snp_files[i],
            gene_anno   = gene_anno
        }
    }
    
    call move_file {
        input:
        disk_size   = disk_size, 
        memory_size = memory_size, 
        preemptible = preemptible,
        cts_files   = TReC_ASReC.output_cts,
        tissue      = tissue
    }
    
    output {
        Array[File] final_cts_files = move_file.trec_asrec_files
    }
}

task TReC_ASReC {
    Int disk_size
    Int n_CPU
    String memory_size
    Int preemptible

    File file_bam
    File file_bai
    File file_hetSNP
    File gene_anno

    String sam_name = sub(basename(file_hetSNP), "\\.txt", "")
    
    command {
        Rscript --vanilla /get_TReC_ASReC.R ${file_bam} ${gene_anno} ${sam_name} ${file_hetSNP} 
    }
    
    output {
        File output_cts = "${sam_name}.trecase.txt"
    }
    
    runtime {
        continueOnReturnCode: false
        docker: "sunway1999/bioconductor_trecase:0.1"
        memory: memory_size
        cpu: n_CPU
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible 
    }

    meta {
        author: "Wei Sun"
    }

}

task move_file {
    Int disk_size
    String memory_size
    Int preemptible
    String tissue

    Array[File] cts_files
    
    command {
        mkdir ${tissue}
        for file1 in ${sep=' ' cts_files}; do
            mv $file1 ${tissue}/${tissue}_${file1}
        done
    }
    
    output {
        Array[File] trec_asrec_files = glob("${tissue}/*.*")
    }
    
    runtime {
        continueOnReturnCode: false
        docker: "sunway1999/bioconductor_trecase:0.1"
        memory: memory_size
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible 
    }

    meta {
        author: "Wei Sun"
    }

}
