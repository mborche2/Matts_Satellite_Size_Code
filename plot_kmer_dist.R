library(ggplot2)
setwd("/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/")

for (i in 1:21){
  filenam <- paste("plots/kmer_frequency_asp/kmer_frequency_",toString(i),"_array.tsv",sep = "")
  array_specifics <- read.table(filenam,header = FALSE)
  freq_table <- table(array_specifics[,2])
  freq_df <- as.data.frame(freq_table)
  freq_df[,2] <- log(freq_df[,2])
  filenam2 <- paste("kmer_frequency_",toString(i),"_array_standard_axes.png",sep="")
  png(filenam2)
  plot(freq_df,xlab="Number of Occurences of Kmer in Array",ylab="Log (ln) Frequency of Kmer Occurence Number",axes=FALSE)
  axis(side=1, at=seq(0,4000,by=25))
  axis(side=2, at=seq(0, 10, by=1))
  box()
  dev.off()
}
