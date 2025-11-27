library(data.table)
library(dplyr)
library(readr)
library(limma)
library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No input folder provided", call. = FALSE)
} else {
  INPUT_DIR <- args[1]
  cat("Input folder:", INPUT_DIR, "\n")
}

# make results folder
OUT_DIR<- "results/fragment_stats"
system(paste("mkdir -p", OUT_DIR))

#stats
### import files
FILES<- list.files(path =INPUT_DIR, pattern ="_length_freq.tsv" )%>%unlist()
#SAMPLES<- substr(files , 1 , nchar(files)-16)
#mean, median, SD, N, min, max
summary<- data.frame( sample=character(0),mean=numeric(0), median=numeric(0),
                      sd=numeric(0),N=integer(0), min=integer(0), max=integer(0))
for (file in FILES){
  ###read length frequency file
  tsv<- read_tsv(file=paste0(INPUT_DIR,"/",file))
  sample<-  substr(file , 1 , nchar(file)-16)
  Mean=sum(tsv$lengths*tsv$Count)/sum(tsv$Count)
  frags<- c()
  for ( i in 1:nrow(tsv))
    frags<- c(frags, rep(tsv$lengths[i], tsv$Count[i]))
  Median<- median(frags)
  SD<- sd(frags, na.rm = TRUE)
  tmp<- data.frame(sample=sample, mean=Mean, median=Median, sd=SD,
                   N=sum(tsv$Count), min=min(tsv$lengths),
                   max=max(tsv$lengths))
  summary<- rbind(summary, tmp)
  
  #histogram
  df<- data.frame(length=frags)
  pdf(file=paste0(OUT_DIR,"/hist_", sample,".pdf"), width = 5, height = 5)
  print(
    ggplot(data=df,aes(x=length))+
      geom_histogram(binwidth=1)
  )
  dev.off()

}

#summary of all data set
write_tsv(summary, file = paste0(OUT_DIR,"/summary.tsv"))


