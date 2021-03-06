#!/usr/bin/env Rscript
# pass the pheno data csv file as the only argument 
args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
# assume no output directory argument was given to rnaseq_pipeline.sh
#  pheno_data_file <- paste0(getwd(), "/chrX_data/geuvadis_phenodata.csv")
#} else {
#  pheno_data_file <- args[1]
#}

pheno_data_file <- args[1]

library(ballgown)
library(genefilter)
library(dplyr)
library(Cairo)

## Read phenotype sample data
pheno_data <- read.csv(pheno_data_file)

## Read in expression data
bg_chrX <- ballgown(dataDir = "/data/read-coverage-tables", samplePattern="ERR", pData=pheno_data)

## Filter low abundance genes
bg_chrX_filt <- subset(bg_chrX, "rowVars(texpr(bg_chrX)) > 1", genomesubset=TRUE)

## DE by transcript
results_transcripts <-  stattest(bg_chrX_filt, feature='transcript', covariate='sex', 
         adjustvars=c('population'), getFC=TRUE, meas='FPKM')

## DE by gene
results_genes <-  stattest(bg_chrX_filt, feature='gene', covariate='sex', 
         adjustvars=c('population'), getFC=TRUE, meas='FPKM')

## Add gene name
results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),
          geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)

## Sort results from smallest p-value
results_transcripts <- arrange(results_transcripts, pval)
results_genes <-  arrange(results_genes, pval)

## Write results to CSV
write.csv(results_transcripts, "transcripts_results.csv", row.names=FALSE)
write.csv(results_genes, "genes_results.csv", row.names=FALSE)

## Filter for genes with q-val <0.05
subset(results_transcripts, results_transcripts$qval <=0.05)
subset(results_genes, results_genes$qval <=0.05)

## Plotting setup
#tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
#palette(tropical)

## Plotting gene abundance distribution
fpkm <- texpr(bg_chrX, meas='FPKM')
fpkm <- log2(fpkm +1)
CairoPDF("gene_abundance_distribution.pdf", width=5, height=5)
par(mar=c(4,4,1,1)+0.1)
boxplot(fpkm, col=c("#4DAF4A","#984EA3"), las=2,ylab='log2(FPKM+1)')
dev.off()

## Plot individual transcripts
#ballgown::transcriptNames(bg_chrX)[12]
#plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),
#     main=paste(ballgown::geneNames(bg_chrX)[12], ' : ',ballgown::transcriptNames(bg_chrX)[12]),
#     pch=19, xlab="Sex", ylab='log2(FPKM+1)')
#points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))

## Plot gene of transcript 1729
#plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX,
#                main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

## Plot average expression
#plotMeans(ballgown::geneIDs(bg_chrX)[203], bg_chrX_filt, groupvar="sex", legend=FALSE)

