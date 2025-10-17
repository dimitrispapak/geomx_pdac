suppressPackageStartupMessages({
	library(Biobase)
	library(BisqueRNA)
	library(Seurat)
	library(knitr)
	library(NanoStringNCTools
	library(GeomxTools)
})
# Deconvolution of CD68 segments unsing BisqueRNA, utilizes approx 80G of RAM
# momac_verse from paper https://www.cell.com/immunity/article/S1074-7613(21)00293-4/fulltext
# https://nextcloud.gustaveroussy.fr/s/LZNjJmsFfwQzSW9

set.seed(42)
geomx_dataset = 'target_demoData.RData'
momac_verse_file = "/mnt/backup/Workspace/phd/geomx/data/2021_MoMac_VERSE.RDS"
load(geomx_dataset)

cd68_subset <- target_demoData[,pData(target_demoData)$segment == 'CD68' & pData(target_demoData)$pancreas_colon == 'pancreas']
bulk = assayDataElement(cd68_subset, elt = "q_norm" )
metadata = pData(cd68_subset)

bulk.matrix = as.matrix(bulk)
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)

print(1)
momac_seurat = readRDS(momac_verse_file)
print(ncol(momac_seurat))
num_cells <- 30000  
subset_cells <- sample(colnames(momac_seurat), num_cells)
momac_seurat <- subset(momac_seurat, cells = subset_cells)
print(ncol(momac_seurat))
expr_matrix <- GetAssayData(momac_seurat, assay = "RNA", slot = "data")

print(2)
cell_metadata <- momac_seurat@meta.data
print(colnames(cell_metadata))
print(3)
sample.ids = as.character(rownames(cell_metadata))
print(4)
individual.labels = as.character(cell_metadata$`Patient No`)
print(5)
cell.type.labels = as.character(cell_metadata$Clusters)
print(6)
sc.pheno <- data.frame(check.names=F, check.rows=F,
		       stringsAsFactors=F,
		       row.names=sample.ids,
		       SubjectName=individual.labels,
		       cellType=cell.type.labels)
print(7)
sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
		      row.names=c("SubjectName", "cellType"))

print(8)
sc.pdata <- new("AnnotatedDataFrame",
		data=sc.pheno,
		varMetadata=sc.meta)
print(9)
sc.eset <- Biobase::ExpressionSet(assayData=as.matrix(expr_matrix), phenoData=sc.pdata)
print(10)
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=F)
print(11)
ref.based.estimates <- t(res$bulk.props)
print(12)
write.csv(ref.based.estimates,file = 'cd68_momac_deconvoluted.csv')
print(13)
knitr::kable(ref.based.estimates, digits=2)
