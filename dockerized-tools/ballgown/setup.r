#!/usr/bin/env Rscript

# load utils package
library("utils")

# install packages from CRAN
install.packages(c("dplyr", "RColorBrewer", "RCurl"))

# install packages from bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite(ask=FALSE)
biocLite(c("ballgown","genefilter"), ask=FALSE)
update.packages(ask=FALSE)
