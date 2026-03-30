library(tidyverse)
library(RColorBrewer)
#All txt in one file!!!!Input format:(*.report.txt)
read_kraken_report <- function(file_path) {
  cols <- c("Percentage", "Clade_reads", "Taxon_reads", "Rank", "TaxID", "Name")
  sample_name <- str_replace(basename(file_path), "\\.report$|\\.txt$", "")
  
  read_tsv(file_path, col_names = cols, show_col_types = FALSE) %>%
    filter(Rank %in% c("S", "U")) %>%
    mutate(
      Name = str_trim(Name),
      Sample = sample_name
    ) %>%
    select(Sample, Name, Percentage, Rank)
}

#Set input PATH here
#Set output PATH at the end of the code
report_files <- list.files(path = "D:/Bishe/1ng/project/xyz/20260324_2", pattern = "*.report", full.names = TRUE)
df_all <- map_dfr(report_files, read_kraken_report)

df_unclass <- df_all %>% filter(Rank == "U") %>% select(-Rank)
df_species <- df_all %>% filter(Rank == "S")

top5_species <- df_species %>%
  group_by(Name) %>%
  summarise(Mean_Pct = mean(Percentage), .groups = "drop") %>%
  arrange(desc(Mean_Pct)) %>%
  slice_head(n = 5) %>%
  pull(Name)

df_plot <- df_species %>%
  mutate(Species = if_else(Name %in% top5_species, Name, "Other")) %>%
  group_by(Sample, Species) %>%
  summarise(Percentage = sum(Percentage), .groups = "drop") %>%
  bind_rows(df_unclass %>% rename(Species = Name))


level_order <- c(top5_species, "Other", "unclassified")
df_plot$Species <- factor(df_plot$Species, levels = level_order)


my_colors <- c(brewer.pal(5, "Set2"), "#E0E0E0", "#555555")


p <- ggplot(df_plot, aes(x = Sample, y = Percentage, fill = Species)) +

  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", linewidth = 0.2, width = 0.7) +
  

  geom_text(
    aes(label = sprintf("%.1f%%", Percentage)), 
    position = position_dodge(width = 0.8), 
    vjust = -0.5, 
    size = 2.5,
    color = "black"
  ) +
  
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(title = "Relative Abundance", x = "Sample ID", y = "Percentage (%)", fill = "Category") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

print(p)

#Set output PATH at the end of the code
ggsave("D:/Bishe/1ng/project/xyz/20260324_2/kraken_abundance_grouped.pdf", plot = p, width = 20, height = 6, dpi = 300)
