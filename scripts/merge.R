library(ggplot2)
rm (list = ls())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No input folder provided \n arguments:sample_file , OUTDIR", call. = FALSE)
} else {
  sample_file<- args[1]
  OUTDIR<- args[2]

}


samples <- read.table(
  sample_file,
  header = FALSE,
  col.names = c("sample", "group","path")
)

df_list <- lapply(seq_len(nrow(samples)), function(i) {
  s <- samples$sample[i]
  g <- samples$group[i]
  f <- file.path("results", s, "smoothed.txt")
  
  dat <- read.table(f, header = TRUE)
  dat$sample <- s
  dat$group  <- g
  dat
})

all_freq <- do.call(rbind, df_list)

pdf(file.path(OUTDIR, "histogram_all.pdf" ))
ggplot(all_freq, aes(x = length, y = smooth,
                     group = sample, color = group)) +
  geom_line(alpha = 0.7) +
  labs(x = "Fragment length (bp)",
       y = "Smoothed frequency",
       color = "Group") +
  theme_minimal()

dev.off()
