# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_10

# Install required Bioconductor package
RUN R -e 'BiocManager::install("GenomicAlignments")'
RUN R -e 'BiocManager::install("GenomicFeatures")'
RUN R -e 'BiocManager::install("zlibbioc")'
RUN R -e 'BiocManager::install("Rsamtools")'

COPY asSeq_0.99.501.tar.gz .

RUN R -e 'install.packages("asSeq_0.99.501.tar.gz", repos = NULL, type = "source")'

ADD ./get_TReC_ASReC.R get_TReC_ASReC.R

