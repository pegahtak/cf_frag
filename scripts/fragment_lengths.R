library(data.table)
library(dplyr)
library(readr)
# parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Both input and output are required", call. = FALSE)
} else {
  input_file <- args[1]
  output_file <- args[2]
  
  # cat("Input file:", input_file, "\n")
  # cat("Output file:", output_file, "\n")
}


lengths<- fread(file=input_file)%>%unlist()
# lengths<- fread(input = "case1_lengths.txt")%>%unlist()
t <- table(factor(lengths, levels = min(lengths):max(lengths))) %>%
  as.data.frame()
colnames(t) <- c("lengths", "Count")
t$Freq<- t$Count/sum(t$Count)

write_tsv(t, file = output_file)
