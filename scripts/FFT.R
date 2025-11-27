library(ggplot2)
library(zoo)
library(data.table)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("INPUT_DIR required", call. = FALSE)
} else {
  INPUT_DIR <- args[1]
  cat("Input dir:", INPUT_DIR, "\n")
  
}

# make results folder
OUT_DIR<- "results/fragment_stats"

FILES<- list.files(path =INPUT_DIR, pattern ="_length_freq.tsv" )%>%unlist()

for (file in FILES){
  
  tsv<- read_tsv(file=paste0(INPUT_DIR,"/",file))
  sample<-  substr(file , 1 , nchar(file)-16)
  df<- subset(tsv,lengths >= 30 & lengths <= 180 )
  
  sub <- df %>% filter(lengths >= 30, lengths <= 180)
  
  sub <- sub %>%mutate(smooth = rollmean(Freq, k = 100, fill = NA),
           residual = Freq - smooth)
  period_file<-  paste0(OUT_DIR,"/",sample,"_period.txt")
  period<- read.table(file =period_file)%>%unlist()%>% as.vector()%>%as.numeric()%>%round(., digits = 3)
  ten_bp_lines <- data.frame(x = seq(30, 180, by = period))
  
  p1 <- ggplot(sub, aes(x = lengths, y = Freq)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_vline(data = ten_bp_lines, aes(xintercept = x),
               linetype = "dotted", color = "red", alpha = 0.5) +
    labs(
      title = "cfDNA fragment-length distribution (30–180 bp)",
      subtitle = "Red dotted lines: 10 bp periodic spacing",
      x = "Fragment length (bp)",
      y = "Normalized count"
    ) +
    theme_minimal(base_size = 14)
  
  # --- Plot 2: residuals showing oscillation directly ---
  p2 <- ggplot(sub, aes(x = length, y = residual)) +
    geom_line(color = "darkorange", linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Residual oscillation after trend removal",
      x = "Fragment length (bp)",
      y = "Residual (normalized count – smoothed trend)"
    ) +
    theme_minimal(base_size = 14)
  
  
}
# --- restrict to 30–180 bp ---
sub <- df %>% filter(length >= 30, length <= 180)

# --- smooth overall trend with moving average (window = 5 bp) ---
sub <- sub %>%
  mutate(smooth = rollmean(norm_count, k = 5, fill = NA),
         residual = norm_count - smooth)

# --- create 10 bp tick marks ---
ten_bp_lines <- data.frame(x = seq(30, 180, by = 10))

# --- Plot 1: zoomed histogram with 10 bp markers ---
p1 <- ggplot(sub, aes(x = length, y = norm_count)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_vline(data = ten_bp_lines, aes(xintercept = x),
             linetype = "dotted", color = "red", alpha = 0.5) +
  labs(
    title = "cfDNA fragment-length distribution (30–180 bp)",
    subtitle = "Red dotted lines: 10 bp periodic spacing",
    x = "Fragment length (bp)",
    y = "Normalized count"
  ) +
  theme_minimal(base_size = 14)

# --- Plot 2: residuals showing oscillation directly ---
p2 <- ggplot(sub, aes(x = length, y = residual)) +
  geom_line(color = "darkorange", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Residual oscillation after trend removal",
    x = "Fragment length (bp)",
    y = "Residual (normalized count – smoothed trend)"
  ) +
  theme_minimal(base_size = 14)

# --- save plots ---
ggsave("results/plots/fragment_length_zoom_10bp.png", p1, width = 7, height = 4, dpi = 300)
ggsave("results/plots/fragment_length_residual_10bp.png", p2, width = 7, height = 4, dpi = 300)
