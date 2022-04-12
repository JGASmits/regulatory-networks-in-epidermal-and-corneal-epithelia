JS_sce_cell_filter_nmad <- function(sce, nmads) {
  libsize.drop <- isOutlier(as.numeric(sce$total_counts), nmads=nmads, type="both", log=TRUE)
  feature.drop <- isOutlier(as.numeric(sce$total_features_by_counts), nmads=nmads, type="both", log=TRUE)
  MT.drop <- isOutlier(sce$pct_counts_Mt, nmads=nmads, type="higher")
  spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=nmads, type="higher")
  sce <- sce[,!(libsize.drop | feature.drop | spike.drop | MT.drop)]
  return(sce)
}

