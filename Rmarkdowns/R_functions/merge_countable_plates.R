JS_merge_countable_plates <- function(location,location_QC_plate) {
  #merges the output of the celseq2 pipeline into a single UMI countmatrix
  #outputs a QC file on a location specified at the second variable
  library('forcats')
  
  file_list <- list.files(location)#make a file list from all files in the /expr/ folder of the homedir
  scData <- list()
  
  myplots <- vector("list")
  low <- "#313695"
  mid="#FFFFBF"
  high="#A50026" 
  RdYlBu.orig <- colorRampPalette(c(low,mid,high))(91)
  palette <- colorRampPalette(rev(brewer.pal(n = 11,name = "RdYlBu")))(10)
  ERCC.palette <- colorRampPalette(rev(brewer.pal(n = 6,name = "RdYlBu")), bias = 2)(10)
  
  for (file in file_list) {
    file_path <- paste(location,file,'/expr.csv',sep="")
    tmp_data <- read.csv(file_path,sep=",")
    
    plate <- NULL
    plate <- data.frame(Col = 1:24, Row = unlist(lapply(rep(LETTERS[1:16]), rep, 24)))
    plate$Row <- fct_rev(plate$Row)
    #plate$Row <- factor(plate$Row, levels=rev(levels(plate$Row)))
    set.seed(1)
    tmp_data.ercc <- data.matrix(tmp_data[grepl("^ERCC-", tmp_data[,1]),-1])
    
    percent_good_wells<- paste(">1500 unique reads :", round(length(which(colSums(data.matrix(tmp_data))[2:length(tmp_data)]>1500))/length(tmp_data)*100),"%")
    percentage_good_ERCC <- paste(">40 ERCCs :",round(length(which(colSums(tmp_data.ercc)>40))/length(tmp_data)*100),"%")
    
    #Plot the ammount of counts per well
    plate[1:384,'val'] <- 0
    plate[1:(length(tmp_data)-1),'val'] <- log10(colSums(data.matrix(tmp_data)))[2:length(tmp_data)]
    myplots[[paste(file,'_counts')]] <- ggplot(plate, aes(Col, Row , fill=val)) + 
      geom_raster() + scale_fill_gradientn(colours=palette)+
      labs(title = paste('counts of', file),
           subtitle = percent_good_wells)
    
    #plot the ammounts of ERCCs sequenced per well
    plate[1:384,'val'] <- 0
    plate[1:(length(tmp_data)-1),'val'] <- colSums(tmp_data.ercc)
    myplots[[paste(file,'_ERCC')]] <- ggplot(plate, aes(Col, Row , fill=val)) + 
      geom_raster() + scale_fill_gradientn(colours=ERCC.palette)+
      labs(title = paste('ERCCs of', file),
           subtitle = percentage_good_ERCC)
    
    cell_type <- gsub("_","-",file)
    #colnames(tmp_data) <- paste(cell_type,sprintf("%1$d", 0:(length(tmp_data)-1)),colnames(tmp_data),sep=':')
    colnames(tmp_data) <- paste(cell_type,colnames(tmp_data),sep=':')
    scData[[length(scData)+1]] <- tmp_data
    print(paste("done for file:",file))
  }
  
  pdf(location_QC_plate ,width=10,height=6,paper='special')
  print(myplots)
  dev.off()
  
  scData_df  <-  as.data.frame(scData)
  scData_df <-  na.omit(scData_df)
  scData_df <- scData_df %>% distinct(scData_df[,1], .keep_all = TRUE)
  row.names(scData_df)  <-scData_df[,1]
  scData_df[,1] <- NULL
  scData_df <- scData_df[1:(length(scData_df)-1)]
  scData_df <- data.matrix(scData_df)
  return(scData_df)
}
