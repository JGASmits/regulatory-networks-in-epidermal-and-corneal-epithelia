JS_generate_pseudo_bulk <- function(seurat_object,types) {
  #generates pseudo_bulk counttable for each of the ID's from a seurat file give a 
  #vector containing ID's you want pseudobulk from. Natively it looks into the metadata column
  #called: "cell_type", if this is not present change the function to the correct column.

  subset_seurat <- function(seurat_object,type){
    #grab cells from a specific 'cell_type' metadata column and sum their counts to a pseudocount
    j_cells <- as.vector(seurat_object@meta.data$cell_type == type)#select cells to keep 
    small_seurat <- seurat_object[,j_cells]
    mat <- as.data.frame(small_seurat@assays$RNA@counts)#maka a df of the counttable from the seurat object
    row.names(mat) <- row.names(small_seurat@assays$RNA@counts)
    mat <- rowSums(mat)#summ all reads into pseudobulk
    return(mat)
  }
  
  counts <- lapply(types, function(x)  subset_seurat(seurat_object, x))
  counts <- as.data.frame(counts)
  colnames(counts) <- types
  return(counts)
}