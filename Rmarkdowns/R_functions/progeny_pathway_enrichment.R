Progeny_enrichment <- function (object, annotation_column_dds = 'genotype', control_annotation ='CTR', comp_annotation = 'all', output_pdf = '~/PROGENy_barplot_analysis.pdf', output_matrix = '~/progeny_data.csv') {
  #object = rlog or vst normalized deseq2 object
  #annotation_column_dds = the name of the metadata column in the deseq2 object in which you specify which are control samples
  #control_annotation = the annotation in the #annotation_column_dds that are the control samples
  #comp_annotation = the annotation in the #annotation_column_dds that will be compared with the control samples, if not specified use all other samples
  
  library(progeny)
  library(dplyr)
  library(data.table)
  
  
  expression_matrix <- assay(object)
  pathways <- progeny(expression_matrix, scale = T)
  controls <- object[[annotation_column_dds]] == control_annotation
  if (comp_annotation != 'all') {
    print('notall')
    comp_columns <- object[[annotation_column_dds]] == comp_annotation
    select_rows <- (controls | comp_columns)
    controls <- controls[select_rows]
    pathways <- pathways[select_rows,]
  }
  
  ctl_mean = apply(pathways[controls,], 2, mean)
  ctl_sd = apply(pathways[controls,], 2, sd)
  pathways <- t(apply(pathways, 1, function(x) x - ctl_mean))
  pathways <- apply(pathways, 1, function(x) x / ctl_sd)
  
  #checking pvallue of difference compared to the control samples
  result = apply(pathways, 1, function(x) {
    broom::tidy(lm(x ~ !controls)) %>%
      filter(term == "!controlsTRUE") %>%
      dplyr::select(-term)
  })
  pathwayspval <- mutate(bind_rows(result), pathway=names(result))
  pathways <- as.data.frame(pathways)
  pathways$pathway <- row.names(pathways)
  
  pathways <- merge(pathways,pathwayspval[c('p.value','pathway')])
  pathways <- as.data.frame(pathways)
  
  controls <- append(controls, FALSE)
  controls <- prepend(controls, FALSE)
  
  
  variable_scores <- pathways[!controls]
  variable_scores <- as.data.frame(variable_scores)
  variable_scores_long <- melt(setDT(variable_scores), id.vars = c("pathway","p.value"), variable.name = "sample")
  pathway_plot <- ggplot(variable_scores_long, aes(x= reorder(pathway,value, FUN = median), y=value, fill = p.value)) + geom_boxplot() +    scale_fill_gradient2(
    low="orange",high="grey",mid="grey",trans="log",breaks=c(0,0.01,0.1,1,10,100,1000))+ ylab(paste('0 = the cell type: ', control_annotation, sep = '')) + xlab(paste('mean pathway intensity in the cell type: ',comp_annotation, sep = '' ))
  
  pdf(output_pdf ,width=8,height=8,paper='special') 
  print(pathway_plot)
  dev.off()
  
  write.csv(as.data.frame(pathways), file= output_matrix)
}