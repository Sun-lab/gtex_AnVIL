workflow collect_TReC {
    String memory
    String disks
    Int preemptible


    call TReC {memory = memory, disks = disks, preemptible = preemptible}
}

task TReC {
    File input_bam
    String input_bam
    File bamIndex
    File gene_anno

    String memory
    String disks
    Int preemptible

    command {
        Rscript --vanilla /get_TReC.R ${input_bam}
    }
    
    output {
        File trec_file = "${input_bam}.trec.txt"
    }
    
    runtime {
        continueOnReturnCode: false
        docker: "bioconductor/bioconductor_docker:latest"
        memory: memory
        disks: disks
        preemptible: preemptible 
    }

}
