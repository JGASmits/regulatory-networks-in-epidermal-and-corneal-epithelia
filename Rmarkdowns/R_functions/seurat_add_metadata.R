JS_add_metadata <- function(seurat_object, metadata_file) {
  #Opens a metadata file, and adds each column as a metadata object to the seurat object.
  #first column of the metadata column has to correspond to the samplenames from the celseq2 pipeline
  #seurat_object <- seur_obj
  metadata <- read.csv2(metadata_file, sep = ',')
  library('stringr')
  library('Seurat')
  #first add a column with the file the cell originated from (cutting the barcode from it and storing this seperate)
  #remove columns with no associated barcode
  seurat_object <- seurat_object[,grep("*.BC.*",colnames(seurat_object)),]
  print(paste('adding filename metadata'))
  barcode <- str_split_fixed(colnames(seurat_object), ".BC.", 2)[,2]
  file <- str_split_fixed(colnames(seurat_object), ".BC.", 2)[,1]
  seurat_object <- AddMetaData(seurat_object, barcode, col.name = 'barcode')
  seurat_object <- AddMetaData(seurat_object, file, col.name = 'filename')
  
  metadata_list <- colnames(metadata)[-1]
  #map the metadata column based on the filename mapped before.
  for (datatype in metadata_list){
      print(paste('adding metadata column ',datatype))
      tmp_metadata <- plyr::mapvalues(x = seurat_object$filename, from = as.vector(metadata[,1]), 
      to = as.vector(metadata[,which(colnames(metadata)==datatype)]))
      seurat_object <- AddMetaData(seurat_object, tmp_metadata, col.name = datatype)
      }
  return(seurat_object)
}
