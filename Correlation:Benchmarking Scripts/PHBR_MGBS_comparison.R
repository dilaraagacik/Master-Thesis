# PHBR and MGBS comparison

# Load required libraries
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ConsensusTME)
library(edgeR)
library(scales)
library(pheatmap)
library(broom)
library(tidyr)
library(ggpubr)

# Load data
phbr_all_tcga <- readRDS("~/Desktop/phbr/phbr_merged_all_samples.rds")
immune_cells <- readRDS("~/Desktop/ConsensusTME_results_geneCounts.rds")
meta_data <- readRDS("~/Desktop/patient_info_from_uuid_SERVER_ALL_tcga.rds")
MGBS <- readRDS("~/Desktop/MGBS_TCGA.rds")

# Correct barcode_to_cancer mapping
barcode_to_cancer <- setNames(meta_data$cancer_type, meta_data$submitter_id)

phbr_joined <- phbr_all_tcga %>%
  mutate(cancer_type = barcode_to_cancer[sample_id]) %>%
  group_by(sample_id, cancer_type) %>%
  summarise(phbr_mhc1 = median(phbr_mhc1, na.rm = TRUE),
            phbr_mhc2 = median(phbr_mhc2, na.rm = TRUE),
            .groups = "drop")

# Clean MGBS patient ID
MGBS$patient <- gsub("\\.[0-9]+$", "", MGBS$patient)

# Build MGBS1 dataframe
mgbs1_df <- MGBS %>%
  dplyr::rename(sample_id = patient) %>%
  mutate(cancer_type = barcode_to_cancer[sample_id]) %>%
  select(sample_id, cancer_type, MGBS1)

combined_scores_mhci <- phbr_joined %>%
  select(sample_id, cancer_type, phbr_mhc1) %>%
  inner_join(mgbs1_df, by = c("sample_id", "cancer_type")) %>%
  mutate(MGBS_score = MGBS1, PHBR_score = phbr_mhc1,
         MGBS_type = "MGBS-I", PHBR_type = "phbr_mhc1") %>%
  select(sample_id, cancer_type, MGBS_score, PHBR_score, MGBS_type, PHBR_type)

# Build MGBS2 dataframe
mgbs2_df <- MGBS %>%
  dplyr::rename(sample_id = patient) %>%
  mutate(cancer_type = barcode_to_cancer[sample_id]) %>%
  select(sample_id, cancer_type, MGBS2)

combined_scores_mhcii <- phbr_joined %>%
  select(sample_id, cancer_type, phbr_mhc2) %>%
  inner_join(mgbs2_df, by = c("sample_id", "cancer_type")) %>%
  mutate(MGBS_score = MGBS2, PHBR_score = phbr_mhc2,
         MGBS_type = "MGBS-II", PHBR_type = "phbr_mhc2") %>%
  select(sample_id, cancer_type, MGBS_score, PHBR_score, MGBS_type, PHBR_type)

# Combine both
combined_scores_separate <- bind_rows(combined_scores_mhci, combined_scores_mhcii)

# Correlation calculation across ALL cancer types
cor_results_all <- combined_scores_separate %>%
  group_by(cancer_type, MGBS_type, PHBR_type) %>%
  summarise(r = cor(PHBR_score, MGBS_score, method = "pearson"),
            p = cor.test(PHBR_score, MGBS_score, method = "pearson")$p.value,
            .groups = "drop") %>%
  filter(!is.na(r), !is.na(p)) %>%
  mutate(adj_p = p.adjust(p, method = "fdr"))

# Filter for CHOL and prepare plot data
plot_data <- combined_scores_separate %>%
  filter(cancer_type %in% c("CHOL")) %>%
  mutate(color_group = case_when(
    MGBS_type == "MGBS-I" ~ "MHC-I",
    MGBS_type == "MGBS-II" ~ "MHC-II"
  ))

manual_colors <- c("MHC-I" = "#1f78b4", "MHC-II" = "#e31a1c")

# Prepare correlation labels with proper y position
cor_results_plot <- cor_results_all %>%
  filter(cancer_type %in% c("CHOL")) %>%
  mutate(
    label = paste0("r = ", round(r, 2), "\nFDR = ", signif(adj_p, 2),"\np-value = ", signif(p, 2)),
    x = 0.5,
    y = case_when(
      MGBS_type == "MGBS-I" ~ 0.028,
      MGBS_type == "MGBS-II" ~ 0.28
    )
  )
# MHC-I plot
p1 <- ggplot(plot_data %>% filter(MGBS_type == "MGBS-I"),
             aes(x = PHBR_score, y = MGBS_score, color = color_group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 1.1) +
  facet_grid(MGBS_type ~ cancer_type) +
  geom_text(
    data = cor_results_plot %>% filter(MGBS_type == "MGBS-I"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 5,
    color = "black",
    hjust = 0,
    vjust = 1
  ) +
  scale_y_continuous(limits = c(0, 0.03)) +
  scale_color_manual(values = manual_colors) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",   # REMOVE LEGEND
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
  ) +
  labs(
    title = "MHC-I",
    x = "PHBR Score",
    y = "MGBS Score"
  )

# MHC-II plot
p2 <- ggplot(plot_data %>% filter(MGBS_type == "MGBS-II"),
             aes(x = PHBR_score, y = MGBS_score, color = color_group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 1.1) +
  facet_grid(MGBS_type ~ cancer_type) +
  geom_text(
    data = cor_results_plot %>% filter(MGBS_type == "MGBS-II"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 5,
    color = "black",
    hjust = 0,
    vjust = 1
  ) +
  scale_y_continuous(limits = c(0, 0.30)) +
  scale_color_manual(values = manual_colors) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",   # REMOVE LEGEND
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
  ) +
  labs(
    title = "MHC-II",
    x = "PHBR Score",
    y = "MGBS Score"
  )

# Combine side by side with main title
ggarrange(p1, p2, ncol = 2, widths = c(1, 1)) %>%
  annotate_figure(
    top = text_grob("PHBR and MGBS Correlation in CHOL Cancer Type", face = "bold", size = 20)
  )
ggsave("PHBR_MGBS_CHOL_correlation.png", width = 16, height = 8, dpi = 300)





