library(data.table)
rm (list = ls())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No input folder provided \n arguments: input, output, sample", call. = FALSE)
} else {
  frag_file<- args[1]
  OUTDIR<- args[2]
  sample<- args[3]
}



# Read 
frag_raw <- fread(frag_file, header = FALSE, col.names = "length")
# sanity check
frag_raw <- frag_raw[length > 0 & length <= 1000]  
summary(frag_raw$length)

# Define bin range 
min_len <- min(frag_raw$length)
max_len <- max(frag_raw$length)

# Create a factor over all possible lengths to get 1-bp bins
len_levels <- min_len:max_len
tab <- table(factor(frag_raw$length, levels = len_levels))

frag <- data.frame(length = len_levels,count  = as.numeric(tab))
frag$freq <- frag$count / sum(frag$count)

## Subset fragments in 100–180 bp range
sub <- frag[frag$length >= 100 & frag$length <= 180, ]
# Use frequency as the signal and center it
x <- sub$freq
x <- x - mean(x)

# FFT
fft_res <- fft(x)
N <- length(x)

# Frequency in cycles per bp (index step = 1 bp)
freqs <- (0:(N - 1)) / N

# Power spectrum
power <- Mod(fft_res)^2

# Use only positive frequencies up to Nyquist
idx <- 1:floor(N / 2)
freqs_pos <- freqs[idx]
power_pos <- power[idx]

# Convert frequency to period in bp
period_bp <- 1 / freqs_pos

# Focus on periods between 5 and 15 bp
period_window <- period_bp >= 5 & period_bp <= 15

cand <- data.frame(
  period_bp = period_bp[period_window],
  power     = power_pos[period_window]
)

# Inspect strongest periodicity in this range
# cand<-cand[order(-cand$power), ][1:min(10, nrow(cand)), ]
cand<- cand[order(cand$power, decreasing = T),]
write.table(cand,file.path(OUTDIR, "period_power.txt"))
#plot periodogram 
maxPower<- cand$period_bp[1]
pdf(file.path(OUTDIR, "periodogram.pdf"))
plot(
  period_bp, power_pos,
  type = "l",
  xlim = c(5, 30),
  xlab = "Period (bp)",
  ylab = "Power",
  main = "Periodogram of cfDNA (100–170 bp window)")
abline(v = maxPower, lty = 2, col="red")

dev.off()


## ---------de-trending----#
L   <- sub$length
# Remove slow envelope with centered moving average (k_trend = 21 bp)
k_trend <- 21
trend <- stats::filter(x, rep(1 / k_trend, k_trend), sides = 2)
x_detr <- x - trend
x_detr[is.na(x_detr)] <- 0  # handle edges

pdf(file.path(OUTDIR, "residual_distribution.pdf"))
plot(L, x, type = "l", xlab = "Length (bp)", ylab = "Freq",
     main = "100–170 bp: raw, envelope, residual")
lines(L, trend, col = "red")
lines(L, x_detr, col = "blue")
legend("topleft",
       legend = c("raw", "envelope (trend)", "residual (detrended)"),
       col = c("black", "red", "blue"), lty = 1, bty = "n")
dev.off()

# Apply Hann window to reduce spectral leakage
N <- length(x_detr)
w <- 0.5 - 0.5 * cos(2 * pi * (0:(N - 1)) / (N - 1))
xw <- x_detr * w

# FFT
fft_res <- fft(xw)

# Frequency axis (cycles per bp)
freqs <- (0:(N - 1)) / N

# Power spectrum
power <- Mod(fft_res)^2

# Use positive frequencies up to Nyquist (skip k = 0)
idx <- 1:floor(N / 2)
freqs_pos <- freqs[idx]
power_pos <- power[idx]

# Convert to period in bp
period_bp <- 1 / freqs_pos

# Focus on realistic nucleosomal oscillation range (8–14 bp)
mask <- period_bp >= 8 & period_bp <= 14

best_idx    <- which.max(power_pos[mask])
best_period <- period_bp[mask][best_idx]
best_power  <- power_pos[mask][best_idx]

cat("Estimated dominant 10-bp-like period in 100–170 bp window:",
    round(best_period, 2), "bp\n")

# Plot periodogram for inspection
pdf(file.path(OUTDIR, "periodogram_detrend.pdf"))
plot(
  period_bp, power_pos,
  type = "l",
  xlim = c(5, 30),
  xlab = "Period (bp)",
  ylab = "Power",
  main = "Periodogram (detrended, windowed, 100–170 bp)",
  sub= paste0("best period at ", round(best_period, 2))
)
abline(v = best_period, lty = 3, col="red")
dev.off()
##------------------------------------------------
## 3. Detect nucleosome peaks on full distribution
##    (mono/di/tri ladder)
##------------------------------------------------

# Smooth frequencies with Daniell smoothing (0.5, 1, 0.5)
daniell_weights <- c(0.5, 1, 0.5)
daniell_weights <- daniell_weights / sum(daniell_weights)

smoothed_freq <- stats::filter(frag$freq, daniell_weights, sides = 2)
y <- as.numeric(smoothed_freq)
smooth_df<- data.frame(length=frag$length, smooth=y)
write.table(smooth_df, file.path(OUTDIR, "smoothed.txt"),  quote = F, row.names = F)

# Local maxima: y[i-1] < y[i] > y[i+1]
dy1 <- diff(y)
peak_idx <- which(dy1[-1] < 0 & dy1[-length(dy1)] > 0) + 1

peak_lengths <- frag$length[peak_idx]
peak_freqs   <- y[peak_idx]

peaks_table <- data.frame(
  length = peak_lengths,
  freq   = peak_freqs
)

# Look at major peaks in 120–400 bp (mono + di region)
peaks_120_400 <- peaks_table[
  peaks_table$length >= 120 & peaks_table$length <= 400, 
]

peaks_120_400 <- peaks_120_400[order(-peaks_120_400$freq), ]

cat("\nTop nucleosome-scale peaks (120–400 bp):\n")
print(head(peaks_120_400, 10))

# Plot full distribution with detected peaks
pdf(file.path(OUTDIR, "dristribution_nucleosome.pdf"))
plot(
  frag$length, frag$freq,
  type = "h",
  xlab = "Fragment length (bp)",
  ylab = "Frequency",
  main = "cfDNA length distribution with detected peaks"
)
lines(frag$length, y, col = "blue", lwd = 2)
points(peak_lengths, peak_freqs, pch = 19, col = "darkgreen", lwd = 1)
abline(v = c(167, 334, 501), lty = 2, col = "red", lwd = 4)
legend("topleft",
       legend = c("expected nucleosome peaks", "smoothed frequency", "called peaks"),
       col = c("red", "blue", "darkgreen"), lty = 1, bty = "n")
dev.off()
