JS_generate_sce_from_df <- function(df) {
  #turns a df into a single cell object. Add Mitochondria & ERCC spikein metadata
  sce<-SingleCellExperiment(assays = list(counts=df))
  is.mito <- grepl("^MT-", rownames(sce))
  is.spike <- grepl("^ERCC", rownames(sce))
  is.spike <- grepl("^ERCC-", rownames(sce))
  is.mito <- grepl("^MT-", rownames(sce))
  sce <- perCellQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))
  return(sce)
}
