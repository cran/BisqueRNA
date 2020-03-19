## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "##",
  highlight = TRUE,
  prompt = FALSE,
  results = "markup"
)

## ----load, echo=T, results='hide', message=F, warning=F-----------------------
library(Biobase)
library(BisqueRNA)

## ----eval=FALSE---------------------------------------------------------------
#  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)

## ----eval=FALSE---------------------------------------------------------------
#  sample.ids <- colnames(sc.counts.matrix)
#  # individual.ids and cell.types should be in the same order as in sample.ids
#  sc.pheno <- data.frame(check.names=F, check.rows=F,
#                         stringsAsFactors=F,
#                         row.names=sample.ids,
#                         SubjectName=individual.labels,
#                         cellType=cell.type.labels)
#  sc.meta <- data.frame(labelDescription=c("SubjectName",
#                                           "cellType"),
#                        row.names=c("SubjectName",
#                                    "cellType"))
#  sc.pdata <- new("AnnotatedDataFrame",
#                  data=sc.pheno,
#                  varMetadata=sc.meta)
#  sc.eset <- Biobase::ExpressionSet(assayData=sc.counts.matrix,
#                                    phenoData=sc.pdata)

## ----eval=FALSE---------------------------------------------------------------
#  sc.eset <- BisqueRNA::SeuratToExpressionSet(seurat.obj, delimiter="-", position=2, version="v3")

## ----simulate, echo=FALSE-----------------------------------------------------
set.seed(42)
cell.types <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial Cells")
avg.props <- c(.5, .2, .2, .07, .03)

expr.data <- BisqueRNA::SimulateData(n.ind=2, n.genes=10, n.cells=10, cell.types=cell.types, avg.props=avg.props)
sc.eset <- expr.data$sc.eset
bulk.eset <- expr.data$bulk.eset

## ----example_input------------------------------------------------------------
sampleNames(sc.eset)
sc.eset$SubjectName
sc.eset$cellType
sampleNames(bulk.eset)

## ----example_input_2----------------------------------------------------------
cell.types <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial Cells")
avg.props <- c(.5, .2, .2, .07, .03)
sim.data <- SimulateData(n.ind=10, n.genes=100, n.cells=500, cell.types=cell.types, avg.props=avg.props)
sc.eset <- sim.data$sc.eset[,sim.data$sc.eset$SubjectName %in% as.character(6:10)]
bulk.eset <- sim.data$bulk.eset
true.props <- sim.data$props
markers <- sim.data$markers

## ----cleanup, echo=FALSE------------------------------------------------------
rm(sim.data)

## ----reference_based----------------------------------------------------------
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=TRUE)

## ----ref_results--------------------------------------------------------------
ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)

## ----ref_results_2------------------------------------------------------------
r <- cor(as.vector(ref.based.estimates), 
         as.vector(true.props[row.names(ref.based.estimates),colnames(ref.based.estimates)]))
knitr::knit_print(sprintf("R: %f", r))

## ----marker_example, echo=FALSE-----------------------------------------------
marker.data.frame <- data.frame(gene=paste("Gene", 1:6),
                                cluster=c("Neurons", "Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial Cells"),
                                avg_logFC=c(0.82, 0.59, 0.68, 0.66, 0.71, 0.62))
knitr::kable(marker.data.frame)

## ----marker_based-------------------------------------------------------------
res <- BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers, weighted=F)

## ----marker_results-----------------------------------------------------------
marker.based.estimates <- res$bulk.props
knitr::kable(marker.based.estimates, digits = 2)

## ----marker_comparison--------------------------------------------------------
scaled.true.props <- t(scale(t(true.props)))[rownames(marker.based.estimates),]
r <- cor(as.vector(marker.based.estimates),
         as.vector(scaled.true.props))
knitr::knit_print(sprintf("R: %f", r))

