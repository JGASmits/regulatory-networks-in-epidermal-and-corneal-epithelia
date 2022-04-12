GOenrichment_plot_counts_JS <- function(GO_object, dds, res, output_file, topn_GOterms, topn_genes, dds_grouping_coldata){
  #generates a dotplot fo a GO_object and visualize the counts 
  # give a deseq object, a results table (of which the most sig genes are picked) 
  # generating topn_genes contributing to each of the topn_GOTerms
  #finally, add some metadata to group
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyverse)
  
  GO_df <- as.data.frame(GO_object)
  str_split(GO_df$GeneRatio, '/')
  GO_df <- GO_df %>%
    separate(GeneRatio, c("genes_found", "GO_size"), "/")
  GO_df$genes_found <- as.numeric(GO_df$genes_found)
  GO_df <- GO_df[order(GO_df$genes_found, decreasing = T),]
  rownames(GO_df) <- NULL
  pdf(output_file ,width=10,height=5,paper='special')
  print(dotplot(GO_object,showCategory = topn_GOterms))#print overview dotplot
  
  if(dim(GO_df)[1] > topn_GOterms){max_examples <- topn_GOterms}else{max_examples <- dim(GO_df)[1]}
  
  for (i in 1:max_examples){
    count_plot <- NULL
    genes_enriched <- str_split(GO_df[i,]$geneID,'/')[[1]]
    genes_enriched_res <- res[genes_enriched,]
    top <- rownames(head(genes_enriched_res[order(genes_enriched_res$padj),],topn_genes))  #take top X genes associated to the GO_term
    tcounts <- t(log2((counts(dds[top, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
      merge(colData(dds), ., by="row.names") %>%
      gather(gene, expression, (ncol(.)-length(top)+1):ncol(.))
    
    count_plot <- ggplot(tcounts, aes_string(dds_grouping_coldata, 'expression', fill = dds_grouping_coldata)) + 
      geom_boxplot() + 
      facet_wrap(~gene, scales="free_y") + 
      labs(x=GO_df[i,]$Description, 
           y="Expression (log normalized counts)", 
           title=paste0("Top GO-term genes ",GO_df[i,]$Description))
    
    print(count_plot)}
  dev.off()
}

