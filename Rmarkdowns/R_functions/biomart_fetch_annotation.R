JS_fetch_biomart_annotation <- function(gene_list, location_biomart_annotation_file){
  #takes a file location & list of genes(hgnc names!) as input, fetches the go_id and go name from the biomart server
  #gene_list<- row.names(counts)
  #location_biomart_annotation_file <- '/R_analysis/Biomart_hcnc2GO.csv2'
  #x_cores <- 40
  if (!file.exists(location_biomart_annotation_file)){
    print(paste('Biomart object ',location_biomart_annotation_file,'does not yet exists, fetching data'))
    print('this might take a while & crash depending on the ammount of trafic on the biomart server, it sucks but I cant fix that')
    library("biomaRt")
    library("dplyr")
    library("stringr")
    library('foreach')
    library('doParallel')
    options(future.globals.maxSize = 4000 * 1024^2)
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl",host = "https://asia.ensembl.org/", ensemblRedirect = FALSE)
    G_list <- vector("list")
    genes <- gene_list
    G_list <- NULL
    G_list <- getBM(attributes=c("go_id","hgnc_symbol"),filters=c('hgnc_symbol'), values=genes, mart=ensembl)
    saveRDS(G_list, location_biomart_annotation_file)
    return(G_list)
    } else {print(paste('Biomart object ',location_biomart_annotation_file,' already present,loaded data from there'))
    return(readRDS(location_biomart_annotation_file))
      }
}


