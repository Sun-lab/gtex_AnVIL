# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

# Install required Bioconductor package
RUN R -e 'BiocManager::install("GenomicAlignments")'
RUN R -e 'BiocManager::install("GenomicFeatures")'
RUN R -e 'BiocManager::install("Rsamtools")'

ADD ./get_TReC.R get_TReC.R

