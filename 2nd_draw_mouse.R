# Setup
work_dir <- "D:/Bishe/1ng/project/wmx/20260330"
setwd(work_dir)

# Configuration
AVERAGE_READ_LEN <- 150  
date_prefix <- "20260330"

library("ggplot2")
library("tidyr")
library("dplyr")

# Functions & Theme
cv <- function(x) { sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) }

gini_index <- function(x) {
  if (all(x == 0)) return(0)
  x <- sort(x)
  n <- length(x)
  G <- sum((2 * (1:n) - n - 1) * x)
  G / (n * sum(x))
}

calculate_metric_autosomes <- function(data_list, metric_func, metric_name, 
                                       autosomes = paste0("chr", 1:19)) {
  result_list <- list()
  for (sample_name in names(data_list)) {
    df <- data_list[[sample_name]]
    autosome_depths <- df[df$V1 %in% autosomes, ]$V3
    val <- metric_func(autosome_depths)
    res <- data.frame(sample = sample_name, chromosome = "autosomes")
    res[[metric_name]] <- val
    result_list[[length(result_list) + 1]] <- res
  }
  return(do.call(rbind, result_list))
}

common_theme <- theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5) 
  )

# Coverage Analysis
cat(">>> Processing Coverage...\n")
cov_file <- paste0(date_prefix, ".coverage.count")

if(file.exists(cov_file)){
  coverage_df <- read.table(cov_file, header = F, sep = "", fill = TRUE)
  sample_ids <- sort(unique(coverage_df$V1)) 
  coverage_df$V1 <- factor(coverage_df$V1, levels = sample_ids)
  
  w_base <- max(14, length(sample_ids) * 0.6) 
  w_group <- max(14, length(sample_ids) * 0.8)
  
  coverage_df$percentage <- coverage_df$V2 / coverage_df$V3
  coverage_df$label <- paste0(round(coverage_df$percentage * 100, 1), "%")
  
  p_coverage <- ggplot(data = coverage_df, mapping = aes(x = V1, y = percentage, fill = V1, label = label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    labs(x = "Sample", y = "Percentage", title = "Genome Coverage") + 
    common_theme + theme(legend.position = "none")
  
  print(p_coverage)
  ggsave(plot = p_coverage, filename = paste0(date_prefix, ".coverage.pdf"), width = w_base, height = 8)
  ggsave(plot = p_coverage, filename = paste0(date_prefix, ".coverage.png"), width = w_base, height = 8)
} else {
  stop("Coverage file not found!")
}

# Relative Depth Analysis (CV & Gini)
cat(">>> Processing Relative Depth (CV & Gini)...\n")
file_paths <- paste0(sample_ids, ".pe.F904", ".s.bam.100000.relative.depth")

if(all(file.exists(file_paths))) {
  relative_depth_list <- setNames(lapply(file_paths, function(fp) read.table(fp, header = FALSE, sep = "\t")), sample_ids)
  
  # CV Calculation
  cv_df <- calculate_metric_autosomes(relative_depth_list, cv, "cv")
  cv_df$platform <- "MDA"
  cv_df$sample <- factor(cv_df$sample, levels = sample_ids) 
  
  p_cv <- ggplot(data = cv_df, aes(x = sample, y = cv, fill = platform)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(mapping = aes(label = round(cv, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    labs(x = "", y = "Coefficient of Variation (CV)", title = "CV of Autosomes Relative Depth") +
    scale_fill_brewer(palette = "YlGnBu") + common_theme + theme(legend.position = "none")
  
  ggsave(plot = p_cv, filename = paste0(date_prefix, ".autosomes_cv.pdf"), width = w_base, height = 6)
  ggsave(plot = p_cv, filename = paste0(date_prefix, ".autosomes_cv.png"), width = w_base, height = 6)
  
  # Gini Calculation
  gini_df <- calculate_metric_autosomes(relative_depth_list, gini_index, "gini")
  gini_df$platform <- "MDA"
  gini_df$sample <- factor(gini_df$sample, levels = sample_ids) 
  
  p_gini <- ggplot(data = gini_df, aes(x = sample, y = gini, fill = platform)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(mapping = aes(label = round(gini, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    labs(x = "", y = "Gini Coefficient", title = "Gini Index of Autosomes Relative Depth") +
    scale_fill_brewer(palette = "YlGnBu") + common_theme + theme(legend.position = "none")
  
  ggsave(plot = p_gini, filename = paste0(date_prefix, ".autosomes_gini.pdf"), width = w_base, height = 6)
  ggsave(plot = p_gini, filename = paste0(date_prefix, ".autosomes_gini.png"), width = w_base, height = 6)
}

# Chimeric Ratios Analysis
cat(">>> Processing Chimeric Ratios...\n")
chimera_file <- paste0(date_prefix, ".Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.count")

if(file.exists(chimera_file)) {
  chimera_df <- read.table(chimera_file, header = F, sep = "\t")
  chimera_df$V1 <- factor(chimera_df$V1, levels = sample_ids) 
  
  chimera_df$total_fragments <- chimera_df$V5 / 2 
  chimera_df$ratio <- chimera_df$V4 / chimera_df$total_fragments
  chimera_df$label <- paste0(round(chimera_df$ratio * 100, 1), "%")
  
  # All Types (Including Normal)
  p_chimera <- ggplot(data = chimera_df, mapping = aes(x = V1, y = ratio, fill = V3, label = label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
    labs(x = "Sample", y = "Percentage", title = "Chimeric Reads Percentage (All Types)") + common_theme
  
  ggsave(plot = p_chimera, filename = paste0(date_prefix, ".chimera_all.pdf"), width = w_group, height = 8)
  ggsave(plot = p_chimera, filename = paste0(date_prefix, ".chimera_all.png"), width = w_group, height = 8)
  
  # Excluding Normal for plot
  chimera_no_normal_df <- chimera_df %>% filter(V3 != "Normal")
  p_chimera_noNormal <- ggplot(data = chimera_no_normal_df, mapping = aes(x = V1, y = ratio, fill = V3, label = label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
    labs(x = "Sample", y = "Percentage", title = "Chimeric Reads Percentage (Excl. Normal)") + common_theme
  
  ggsave(plot = p_chimera_noNormal, filename = paste0(date_prefix, ".chimera_noNormal.pdf"), width = w_group, height = 8)
  ggsave(plot = p_chimera_noNormal, filename = paste0(date_prefix, ".chimera_noNormal.png"), width = w_group, height = 8)
}

# Base Counts & Normalized Breakpoints
cat(">>> Processing Base Counts & Breakpoints...\n")
base_count_file <- paste0(date_prefix, ".total_base.count")

# Total Base Count
if(file.exists(base_count_file)) {
  base_count_raw <- read.table(base_count_file, header = F, sep = "\t")
  base_count_agg <- base_count_raw %>% 
    group_by(V1) %>% summarise(sum_V3 = sum(V3)) %>% 
    mutate(sum_V3_sci = formatC(sum_V3, format = "e", digits = 2)) 
  base_count_agg$V1 <- factor(base_count_agg$V1, levels = sample_ids) 
  
  p_basecount <- ggplot(data = base_count_agg, mapping = aes(x = V1, y = sum_V3, fill = V1, label = sum_V3_sci)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.75)  +
    geom_text(vjust = -0.5, size = 3) +
    labs(x = "Sample", y = "Base Count (Raw)", title = "Sequencing Base Count") +
    common_theme + theme(legend.position = "none")
  
  ggsave(plot = p_basecount, filename = paste0(date_prefix, ".basecount.pdf"), width = w_base, height = 8)
  ggsave(plot = p_basecount, filename = paste0(date_prefix, ".basecount.png"), width = w_base, height = 8)
}

# Breakpoint Normalization
if(exists("chimera_df")) {
  chimera_norm_df <- chimera_df %>%
    mutate(mapped_bases = V5 * AVERAGE_READ_LEN, breakpoint_per_10kb = (V4 / mapped_bases) * 10000) %>%
    filter(V3 != "Normal")
  
  p_breakpoints <- ggplot(data = chimera_norm_df, mapping = aes(x = V1, y = breakpoint_per_10kb, fill = V3, label = round(breakpoint_per_10kb, 2))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
    labs(x = "Sample", y = "Breakpoint per 10kb (Mapped)", title = "Normalized Chimeric Breakpoints") +
    common_theme
  
  print(p_breakpoints)
  ggsave(plot = p_breakpoints, filename = paste0(date_prefix, ".breakpoint_per_10kb.pdf"), width = w_group, height = 8)
  ggsave(plot = p_breakpoints, filename = paste0(date_prefix, ".breakpoint_per_10kb.png"), width = w_group, height = 8)
}

# Split Mapping Analysis
cat(">>> Processing Split Mapping Rates...\n")
split_file <- paste0(date_prefix, ".splitmapping.count")

if(file.exists(split_file)) {
  split_df <- read.table(split_file, header = F, sep = "\t")
  colnames(split_df) <- c("Sample", "SA_Count", "Total_Count")
  split_df$Sample <- factor(split_df$Sample, levels = sample_ids) 
  
  split_df$Rate <- split_df$SA_Count / split_df$Total_Count
  
  split_df$Label <- sprintf("%.2f%%", split_df$Rate * 100)
  
  p_split <- ggplot(data = split_df, mapping = aes(x = Sample, y = Rate, fill = Sample, label = Label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
    geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    labs(x = "Sample", y = "Split Mapping Rate", title = "Split-Read Chimeric Rate") +
    common_theme + theme(legend.position = "none") 
  
  print(p_split)
  ggsave(plot = p_split, filename = paste0(date_prefix, ".splitmapping_rate.pdf"), width = w_base, height = 6)
  ggsave(plot = p_split, filename = paste0(date_prefix, ".splitmapping_rate.png"), width = w_base, height = 6)
}

cat(">>> Analysis Complete.\n")
