# Environment for the paired tumour/normal RNA-seq analysis.
# Base: Bioconductor 3.23 (R 4.6.1), pinned by the multi-architecture index
# digest, so the same file builds on x86_64 and on arm64.
FROM bioconductor/bioconductor_docker:RELEASE_3_23@sha256:821dbf9ac119eac41f177531c7ca8fc7084c99c04eb218a1aa7f52cb63bad5d9

# Packages used by the scripts in R/. Installed as root at build time.
RUN Rscript -e 'BiocManager::install(c("DESeq2", "recount3", "edgeR", "apeglm"), ask = FALSE, update = FALSE, Ncpus = 2)'

# Fail the build if any package is missing, and print the versions installed
# so the build log records them.
RUN Rscript -e 'pk <- c("DESeq2", "recount3", "edgeR", "apeglm"); stopifnot(all(pk %in% rownames(installed.packages()))); for (p in pk) cat(p, as.character(packageVersion(p)), "\n")'

WORKDIR /home/rstudio/project
