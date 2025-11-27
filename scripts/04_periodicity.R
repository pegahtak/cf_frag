library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(stats)
library(tseries)
library(zoo)

rm(list = ls())
# parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("INPUT_DIR required", call. = FALSE)
} else {
  INPUT_DIR <- args[1]
  cat("Input dir:", INPUT_DIR, "\n")

}

# make results folder
OUT_DIR<- "results/fragment_stats"
# system(paste("mkdir -p", OUT_DIR))

#stats
### import files
FILES<- list.files(path =INPUT_DIR, pattern ="_length_freq.tsv" )%>%unlist()

for (file in FILES){
  ###read length frequency file
  tsv<- read_tsv(file=paste0(INPUT_DIR,"/",file))
  sample<-  substr(file , 1 , nchar(file)-16)
  df<- tsv
  
  all_lengths <- min(tsv$lengths):max(tsv$lengths)
  counts <- rep(0, length(all_lengths))
  counts[match(df$lengths, all_lengths)] <- df$Count

  signal <- counts - mean(counts)  # remove DC offset
  N <- length(signal)

  fft_vals <- fft(signal)                # complex FFT
  amplitude <- Mod(fft_vals)[1:(N/2)]    # magnitude spectrum
  freq <- (0:(N/2 - 1)) / N              # normalized frequency (cycles per bp)

  period <- 1 / freq
  power <- amplitude^2

  period_bp <- 1 / freq
  power_df <- data.frame(freq = freq, period_bp = period_bp, amplitude = amplitude)

  # --- Find strongest periodic component (excluding DC, low-freq noise) ---
  power_df <- power_df[is.finite(power_df$period_bp) & power_df$period_bp <= 50, ]  # look up to 50 bp
  peak <- power_df[which.max(power_df$amplitude), ]

  cat(sprintf("Dominant periodicity: %.2f bp\n", peak$period_bp))

  write.table(peak$period, file = paste0(OUT_DIR,"/",sample,"_period.txt"),
              col.names = F, row.names = F, quote = F, sep = "\t")

  # --- 7. Plot FFT amplitude vs. period ---
  pdf(file=paste0(OUT_DIR,"/",sample,"_FFT_plot.pdf"), width = 8, height = 5)
print(
      ggplot(df_fft, aes(x = period, y = amplitude)) +
      geom_line(color = "steelblue", linewidth = 1) +
      geom_vline(xintercept = peak$period, linetype = "dashed", color = "red") +
      labs(title = sprintf("fragment-length distribution (peak ≈ %.1f bp)",
                           peak$period_bp),
      x = "Period (bp per cycle)",
      y = "Amplitude"
         ) +
        theme_minimal(base_size = 14)

    )
  dev.off()
  
  #shendure FFT
  df<- tsv
  counts <- df$Count
  # Daniell window of width 3 with half weights at ends
  daniell_weights <- c(0.5, 1, 0.5)
  daniell_weights <- daniell_weights / sum(daniell_weights)
  
  smoothed <- rollapply(counts, width = 3,
                        FUN = function(x) sum(x * daniell_weights[seq_along(x)]),  # weighted mean
                        align = "center", partial = TRUE)
  x <- df$lengths
  y <- smoothed 
  y_centered <- y - mean(y, na.rm = TRUE)
  fit <- lm(y_centered ~ x)          # fit linear trend
  detrended <- resid(fit)  
  
  # check smoothing
  pdf(file = paste0(OUT_DIR,"/",sample,"_smoothed.pdf"), width = 8, height = 5)
  plot(x, y_centered, type = "l", col = "grey70",
       xlab = "Fragment length (bp)", ylab = "Signal",
       main = "After mean subtraction and linear detrending")
  lines(x, detrended, col = "steelblue", lwd = 2)
  abline(h = 0, lty = "dashed", col = "black")
  dev.off()
  
  signal <- detrended
  # define filter frequencies
  freqs <- 1 / seq(5, 100, 4)
  filtered <- signal
  
 #  #recursive filtering
 #  for (f in freqs) {
 #    k <- round(1 / f)
 #    k <- max(2, min(k, length(filtered)))
 #    coeff <- rep(1 / k, k)
 #    
 #    # apply a one-sided moving average manually
 #    newf <- rep(NA_real_, length(filtered))
 #    for (i in seq_along(filtered)) {
 #      lo <- max(1, i - k + 1)
 #      newf[i] <- mean(filtered[lo:i])
 #    }
 #    
 #    filtered <- newf
 # }
 #  
 #  
 #  # correct for the last 24 elements
 #  filtered <- as.numeric(filtered)
 #  n <- length(filtered)
 #  if (n > 48) {
 #    filtered[(n-23):n] <- filtered[(n-47):(n-24)]
 #  }
 #  # fix for NAs
 #  filtered[is.na(filtered)] <- 0
 # pdf(file = paste0(OUT_DIR,"/",sample,"_recursive_filtered.pdf"), width = 8 , height = 5)
 #  plot(x, signal, type = "l", col = "grey80",
 #       xlab = "Fragment length (bp)", ylab = "Signal",
 #       main = "Recursive filtered cfDNA trajectory")
 #  lines(x, filtered, col = "steelblue", lwd = 2)
 #  abline(h = 0, lty = "dashed", col = "black")
 # dev.off()  

 # histogram plot:
  Mean=sum(tsv$lengths*tsv$Count)/sum(tsv$Count)
  frags<- c()
  for ( i in 1:nrow(tsv))
    frags<- c(frags, rep(tsv$lengths[i], tsv$Count[i]))
  Median<- median(frags)
  SD<- sd(frags, na.rm = TRUE)
  N=sum(tsv$Count)
  Min=min(tsv$lengths)
  Max=max(tsv$lengths)
  
  text_box <- sprintf(
    "Mean: %.2f\nMedian: %.2f\nSD: %.2f\nMin: %d\nMax: %d\nTotal Fragments: %d",
    Mean, Median, SD, Min, Max, N)
  
  df<- data.frame(length=frags)
  pdf(file=paste0(OUT_DIR,"/histogram_", sample,".pdf"), width = 8, height = 5)
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
      )
  )
  dev.off()
}



