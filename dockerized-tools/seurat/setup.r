#!/usr/bin/env Rscript

# load utils package
library("utils")

# install packages from CRAN
install.packages(c("devtools", "dplyr", "RColorBrewer", "RCurl", "Cairo"))

# install Seurat
library("devtools")
install_github("satijalab/seurat")