plot_dist_matrix <- function(normalized_deseq_int, output_file, annotation1 = NA, annotation2 = NA) {
  #normalized_deseq_int = either rlog or vst normalied deseq2 object
  #output_file = location in which the output pdf is generated
  #annotation1 = optional annotation name, should be present as metadata to the normalized_deseq_int object
  #annotation2 = optional second annotation name, should be present as metadata to the normalized_deseq_int object
  library('ComplexHeatmap')
  library("RColorBrewer")
  sampleDists <- dist(t(assay(normalized_deseq_int)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9,"Blues")) )(1000)
  
  if (isNA(annotation1)) {
    myplot <- 
      Heatmap(sampleDistMatrix,
              row_dend_side = "right",
              show_column_dend = FALSE,
              col = colors)}
  else {
    if (isNA(annotation2)){
      ha = rowAnnotation(annotation1 = normalized_deseq_int[[annotation1]])}
    else {ha = rowAnnotation(annotation1 = normalized_deseq_int[[annotation1]],
                             annotation2 = normalized_deseq_int[[annotation2]])}
    myplot <- 
      Heatmap(sampleDistMatrix,
              row_dend_side = "right",
              show_column_dend = FALSE,
              col = colors,
              right_annotation = ha )}
  
  pdf(output_file,width=6,height=6,paper='special') 
  print(myplot)
  dev.off()
}