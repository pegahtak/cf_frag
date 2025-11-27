suppressMessages(library(data.table))
suppressMessages(library(ggplot2))


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No input folder provided", call. = FALSE)
} else {
  SAMPLE <- args[1]
  OUTPUT_DIR <- args[2]
  cat("plotting norm coverage for:", SAMPLE, "\n")
}

# OUTPUT_DIR<- "results/fragment_stats"

dt <- fread(file=
              file.path(OUTPUT_DIR, paste0(SAMPLE,"_regions_normalized_cov.bed")),
            header = TRUE)

# Quick check
colnames(dt)<- c("chr","start","end","gene_id","score","strand","cov_region", "cov_up", "cov_down", "bg_cov", "cov_norm")
#head(dt)

# 1) Histogram of normalized coverage
pdf(file = file.path(OUTPUT_DIR, "regions_hist.pdf"))
print(
ggplot(dt, aes(x = cov_norm)) +
  geom_histogram(bins = 50) +
  xlab("Normalized coverage (region / flanks)") +
  ylab("Number of regions") +
  ggtitle("Distribution of normalized coverage")
)
dev.off()
# 2) Scatter: region vs background coverage
pdf(file = file.path(OUTPUT_DIR,"regions_vs_backgrouns_hist.pdf"))
print(
ggplot(dt, aes(x = bg_cov, y = cov_region)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("Background coverage (mean of flanks)") +
  ylab("Region coverage") +
  ggtitle("Region vs local background coverage")
)
dev.off()
# regions ( promoter, enhancer, etc)
# dt[, type := tstrsplit(name, ":", fixed = TRUE)[[1]]]
# 
# ggplot(dt, aes(x = type, y = cov_norm)) +
#   geom_boxplot(outlier.size = 0.5) +
#   xlab("Region type") +
#   ylab("Normalized coverage") +
#   ggtitle("Normalized coverage by region type") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
