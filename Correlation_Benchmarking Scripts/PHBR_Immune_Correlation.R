# Load libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(dplyr)
library(broom)



phbr_all_tcga = readRDS("~/Desktop/phbr/phbr_merged_all_samples.rds")
immune_cells = readRDS("~/Desktop/ConsensusTME_results_geneCounts.rds")
meta_data = readRDS("~/Desktop/patient_info_from_uuid_SERVER_ALL_tcga.rds")
barcode_to_cancer <- setNames(meta_data$cancer_type, meta_data$submitter_id)


immune_df <- bind_rows(
  lapply(names(immune_cells), function(cancer) {
    mat <- immune_cells[[cancer]]  # rows = immune cell types, columns = samples
    df <- as.data.frame(t(mat))    # transpose to have samples as rows
    df$sample_id <- rownames(df)
    df$cancer_type <- cancer
    df
  }),
  .id = "cancer_index"
)


# Join PHBR scores with immune scores
phbr_joined <- phbr_all_tcga %>%
  mutate(cancer_type = barcode_to_cancer[sample_id]) %>%
  group_by(sample_id, cancer_type) %>%
  summarise(phbr_mhc1 = median(phbr_mhc1, na.rm = TRUE),
            phbr_mhc2 = median(phbr_mhc2, na.rm = TRUE),
            .groups = "drop")


immune_long <- immune_df %>%
  pivot_longer(cols = -c(sample_id, cancer_type, cancer_index), names_to = "cell_type", values_to = "score")

# Merge immune and PHBR data
plot_data <- immune_long %>%
  inner_join(phbr_joined, by = c("sample_id", "cancer_type"))


calc_correlation <- function(df, phbr_column) {
  df %>%
    group_by(cancer_type, cell_type) %>%
    summarise(
      cor = if (sum(is.finite(score) & is.finite(.data[[phbr_column]])) >= 2) {
        cor(score, .data[[phbr_column]], method = "pearson", use = "complete.obs")
      } else { NA },
      pval = if (sum(is.finite(score) & is.finite(.data[[phbr_column]])) >= 3) {
        cor.test(score, .data[[phbr_column]], method = "pearson")$p.value
      } else { NA },
      .groups = "drop"
    )
}


# Calculate correlations
cor_mhc1 <- calc_correlation(plot_data, "phbr_mhc1")
cor_mhc2 <- calc_correlation(plot_data, "phbr_mhc2")
cor_combined <- bind_rows(cor_mhc1, cor_mhc2)
cor_mhc1$score_type <- "MHC I"
cor_mhc2$score_type <- "MHC II"
cor_combined <- bind_rows(cor_mhc1, cor_mhc2)
cor_combined$adj.P <- p.adjust(cor_combined$pval, method = "fdr")
cor_combined$logAdjP <- -log10(cor_combined$adj.P + 1e-300)  # Avoid -Inf
cor_combined$size_cut <- cut(
  cor_combined$logAdjP,
  breaks = c(-Inf, 1, 2, 3, Inf),
  labels = c("≤1", "1–2", "2–3", ">3")
)

# Exclude rows with NA adjusted p-values
cor_combined <- cor_combined %>%
  mutate(score_type = ifelse(is.na(score_type), "unknown", score_type)) %>%
  filter(!is.na(adj.P)) %>%
  mutate(size_cut = cut(
    adj.P,
    breaks = c(-Inf, 1e-3, 1e-2, 5e-2, Inf),
    labels = c("≤0.001", "0.001–0.01", "0.01–0.05", ">0.05")
  ))
# First compute logP for dot size (as you did for MGBS)
cor_combined <- cor_combined %>%
  filter(!is.na(adj.P)) %>%
  mutate(logP = -log10(adj.P))  # avoid -Inf

# Optional: consistent facet labels
facet_labels <- c("MHC I" = "MHC-I", "MHC II" = "MHC-II")

# Define dot size range consistent with MGBS style
dot_size_range <- c(1, 12)

# Define color palette consistent with MGBS style
color_palette <- c("#053061", "#2166ac", "#d1e5f0", "white", "#fddbc7", "#b2182b", "#67001f")

# Now the plot
ggplot(cor_combined, aes(x = cancer_type, y = cell_type, fill = cor, size = adj.P)) +
  geom_point(shape = 21, color = "grey30", stroke = 0.2, alpha = 0.92) +
  facet_wrap(~ score_type, ncol = 2, scales = "free_x",
             labeller = as_labeller(facet_labels)) +
  scale_fill_gradientn(
    colors = color_palette,
    limits = c(-0.4, 0.4),
    oob = scales::squish,
    name = "Pearson r"
  ) +
  scale_size_continuous(
    range = dot_size_range,
    breaks = c(0.2, 0.3, 0.4, 0.7, 0.8, 0.9, 1),
    labels = c("0.2", "0.3", "0.4", "0.7", "0.8", "0.9", "1"),
    trans = "reverse",  # Important: smaller FDR → larger dot!
    name = "FDR"
  ) +
  scale_x_discrete(guide = guide_axis(angle = 40)) +
  scale_y_discrete(expand = expansion(mult = c(0.07, 0.13))) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(size = 26, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    panel.spacing = unit(3, "lines"),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
    legend.key.height = unit(1.5, "cm"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    strip.background = element_rect(fill = "grey96", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) + guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  )+
  labs(
    title = "PHBR–Immune Cells Correlation Across Cancer Types",
    x = "Cancer Type",
    y = "Immune Cell Type"
  )

# Save matching MGBS
ggsave("DotPlot_PHBR_Immune_CancerTypes_adjP_updated.png", width = 30, height = 17, dpi = 300)
ggsave("DotPlot_PHBR_Immune_CancerTypes_adjP_updated.pdf", width = 23, height = 11)

# One panel → prepare the label coordinates for "left" placement:

cor_labels <- cor_combined %>%
  filter(cell_type == "Immune_Score") %>%
  select(score_type, cancer_type, cor, pval, adj.P) %>%
  filter(cancer_type == "CHOL") %>%  # only CHOL
  mutate(PHBR_type = recode(score_type, "MHC I" = "MHC-I", "MHC II" = "MHC-II")) %>%
  mutate(label = paste0(
    "r = ", round(cor, 2), "\nFDR = ", signif(adj.P, 2), "\np = ", signif(pval, 2)
  )) %>%
  mutate(x = -Inf, y = Inf) %>%   # force top-left placement
  select(PHBR_type, label, x, y)

# Create merged PHBR dataframe for plotting
merged_phbr_df <- plot_data %>%
  filter(cell_type == "Immune_Score") %>%  # only Immune_Score
  pivot_longer(
    cols = c(phbr_mhc1, phbr_mhc2),
    names_to = "PHBR_type",
    values_to = "PHBR"
  ) %>%
  mutate(
    PHBR_type = recode(PHBR_type, "phbr_mhc1" = "MHC-I", "phbr_mhc2" = "MHC-II"),
    immune_score = score
  ) %>%
  select(sample_id, cancer_type, PHBR_type, PHBR, immune_score)


# Prepare plot data
plot_df <- merged_phbr_df %>%
  filter(cancer_type == "CHOL")  # only CHOL

# Plot
ggplot(plot_df, aes(x = PHBR, y = immune_score, color = PHBR_type)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 1.2) +
  facet_wrap(~ PHBR_type, ncol = 2, scales = "free_x") +
  geom_text(
    data = cor_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,    # left-top corner
    size = 5,
    color = "black",
    fontface = "bold"
  ) +
  scale_color_manual(values = c("MHC-I" = "#1f78b4", "MHC-II" = "#e31a1c")) +
  labs(
    title = "PHBR and Immune Score Correlation in CHOL",
    x = "PHBR Score",
    y = "Immune Score",
    color = "PHBR Type"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none",
    panel.spacing = unit(1.2, "lines")
  )

# Save
ggsave("ImmuneScore_vs_PHBR_Correlation_CHOL_pretty_left_label.png", width = 12, height = 6, dpi = 300)
