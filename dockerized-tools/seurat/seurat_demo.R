# Seurat demo script

library(Seurat)
library(Matrix)
library(Cairo)

pbmc33k.data <- Read10X("../input")
pbmc33k  <- new("seurat", raw.data = pbmc33k.data)

#setup setting do.scale and do.center to F - this means that we will NOT scale genes by default (to speed things up)
pbmc33k <- Setup(pbmc33k, min.cells = 3, min.genes = 200, project = "10X_PBMC", do.scale = F, do.center = F, names.field = 2,
                 names.delim = "\\-")
mito.genes <- grep("^MT-", rownames(pbmc33k@data), value = T)
percent.mito <- colSums(expm1(pbmc33k@data[mito.genes, ])) / colSums(expm1(pbmc33k@data))

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
pbmc33k <- AddMetaData(pbmc33k, percent.mito, "percent.mito")

#Vln Plot
CairoPDF("VlnPlot.pdf", width=8, height=8)
par(mar=c(4,4,1,1)+0.1)
VlnPlot(pbmc33k, c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

#We filter out cells that have unique gene counts over 2,500 and under 500, and > 5% mitochondrial percentage
pbmc33k <- SubsetData(pbmc33k, subset.name = "nGene", accept.high = 2500)
pbmc33k <- SubsetData(pbmc33k, subset.name = "percent.mito", accept.high = 0.05)
pbmc33k <- SubsetData(pbmc33k, subset.name = "nGene", accept.low = 500)

#mean-variability plot
CairoPDF("MeanVarPlot.pdf", width=5, height=5)
par(mar=c(4,4,1,1)+0.1)
MeanVarPlot(pbmc33k, x.low.cutoff = 0, y.cutoff = 0.8)
dev.off()

