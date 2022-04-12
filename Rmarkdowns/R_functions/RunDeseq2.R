args = commandArgs(trailingOnly=TRUE)


if (length(args)!=3) {
  stop("A countmatrix csv file (first arg), a sample csv file (second arg), design (third arg) and output file location (fourth arg) should be supplied", call.=FALSE)
} else if (length(args)==3) {
  print('using the files supplied')
  print((paste("supplied count_matrix: ", args[1])))
  print((paste("supplied sample_file: ", args[2])))
  print((paste("output_file_loc: ", args[3])))
}


counts <- read.csv2(args[1], header = TRUE, sep = '\t')
row.names(counts) <- counts$loc
regions <- counts$loc
counts$loc <- NULL

counts  <- mapply(counts , FUN=as.numeric)

sample_data<- read.csv2(args[2], header = TRUE, sep = ',')
row.names(sample_data) <- sample_data$samples
sample_data$samples <- NULL


sample_data$condition <- factor(sample_data$condition)

suppressMessages(
library("DESeq2"))


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_data,
                              design =  ~ condition)

dds <- DESeq(dds)
res <- results(dds)

peaks <- as.data.frame(res)
peaks$region <- regions 
write.csv(peaks, file =args[3])
