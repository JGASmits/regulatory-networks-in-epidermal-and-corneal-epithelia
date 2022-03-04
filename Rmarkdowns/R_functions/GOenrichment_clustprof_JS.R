GOenrichment_clustprof_JS <- function(gene_list,gene_background_list, output_file, qval_cutoff = 0.05,gene_ID = 'SYMBOL', GO_ontology="BP", simplify_value = 0.7 ) {
  #run clusterprofiler GO-term enrichment on the genelist
  #simplify's the GO term results 
  #return a list with 1. a GO_term object, 2 a df of the GO_term object
  #write GO_term output table to a csv file in output_file
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyverse)
  if (gene_background_list == "NA"){
    GO_obj <- enrichGO(gene = gene_list,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = gene_ID,
                     ont           = GO_ontology,
                     pAdjustMethod = "BH",
                     qvalueCutoff  = qval_cutoff)
  }else{
  GO_obj <- enrichGO(gene = gene_list,
                     universe = gene_background_list,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = gene_ID,
                     ont           = GO_ontology,
                     pAdjustMethod = "BH",
                     qvalueCutoff  = qval_cutoff)}

  if(nrow(as.data.frame(GO_obj))==0){
    return('no sig go terms')}else{
      GO_obj_s <- clusterProfiler::simplify(GO_obj, cutoff = simplify_value)
      GO_df <- as.data.frame(GO_obj_s)
      str_split(GO_df$GeneRatio, '/')
      GO_df <- GO_df %>% separate(GeneRatio, c("genes_found", "GO_size"), "/")
      GO_df$genes_found <- as.numeric(GO_df$genes_found)
      GO_df <- GO_df[order(GO_df$genes_found, decreasing = T),]
      write.table(as.data.frame(GO_df), file= output_file, sep = ',', row.names = T)
      return(list(GO_obj_s,GO_df))
      }
}
