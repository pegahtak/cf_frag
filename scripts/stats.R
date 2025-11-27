suppressMessages( library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(limma))
suppressMessages(library(ggplot2))


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No input folder provided \n arguments: input, output, sample", call. = FALSE)
} else {
  INPUT<- args[1]
  OUTDIR<- args[2]
  sample<- args[3]
}


#stats
### import files

  tsv<- read_tsv(file= INPUT)
  Mean=sum(tsv$lengths*tsv$Count)/sum(tsv$Count)
  frags<- c()
  for ( i in 1:nrow(tsv))
    frags<- c(frags, rep(tsv$lengths[i], tsv$Count[i]))
  Median<- median(frags)
  SD<- sd(frags, na.rm = TRUE)
  Min=min(tsv$lengths)
  Max=max(tsv$lengths)
  N=sum(tsv$Count)
  summary<- data.frame(sample=sample, mean=Mean, median=Median, sd=SD,
                   N=N, min=Min,
                   max=Max)

  
  text_box <- sprintf(
    "Mean: %.2f\nMedian: %.2f\nSD: %.2f\nMin: %d\nMax: %d\nTotal Fragments: %d",
    Mean, Median, SD, Min, Max, N)
  
  #histogram
  df<- data.frame(length=frags)
  pdf(file.path(OUTDIR,"hist.pdf"), width = 5, height = 5)
  print(
    ggplot(data=df,aes(x=length))+
      geom_histogram(binwidth=1)+
      annotate(
        "text",
        x = Inf, y = Inf,             
        label = text_box,
        hjust = 1.1, vjust = 1.1,
        family = "mono",
        size = 3.5,
        color = "black",
        lineheight = 1.1,
        fontface = "plain"
      )+
      geom_vline(xintercept = c(167, 334, 501), linetype = "dashed", color = "red")+
      labs( x = "Fragment length (bp)",
            y = "Frequency",
            title = "cfDNA fragment length distribution",
            subtitle = "expected mono/di/tri-nucleosome")
  )
  dev.off()
  


#summary of all data set
write_tsv(summary, file.path(OUTDIR, paste0(sample,"_summary.tsv")))


