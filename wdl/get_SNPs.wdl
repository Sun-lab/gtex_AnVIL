# modified version of github repo file under broadinstitute/gtex-pipeline/genotype/participant_vcfs.wdl
# We want to extract all SNPs for eQTL mapping and subset hetero SNPs for ASReC

task get_SNPs {

    File vcf_file
    String SUBJID
		
		# vcf_file, the large GTEx WGS phased vcf.gz file
		# SUBJID, the SUBJID variable from dat_<tissue>.tsv

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        date +"[%b %d %H:%M:%S] Generating participant VCF (SNPs only)"
        # select SNPs, filter out missing sites
        bcftools view --no-update -s ${SUBJID} -v snps ${vcf_file} | bcftools view --no-update -e 'GT=".|."' -Oz -o ${SUBJID}.snps.vcf.gz
        tabix ${SUBJID}.snps.vcf.gz

        date +"[%b %d %H:%M:%S] Subsetting het sites for ASE"
        bcftools view --no-update -i 'GT="het"' -Oz -o ${SUBJID}.snps.het.vcf.gz ${SUBJID}.snps.vcf.gz
        tabix ${SUBJID}.snps.het.vcf.gz

        date +"[%b %d %H:%M:%S] Done"
    }

    output {
        File snps_vcf = "${SUBJID}.snps.vcf.gz"
        File snps_vcf_index = "${SUBJID}.snps.vcf.gz.tbi"
        File snps_het_vcf = "${SUBJID}.snps.het.vcf.gz"
        File snps_het_vcf_index = "${SUBJID}.snps.het.vcf.gz.tbi"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Paul Little"
    }
}


workflow get_SNPs_workflow {
    call get_SNPs
}