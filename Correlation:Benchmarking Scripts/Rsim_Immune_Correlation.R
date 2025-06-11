# Load necessary libraries
library(tidyverse)     # For data manipulation and visualization
library(ggplot2)       # For plotting
library(reshape2)      # For reshaping data
library(scales)        # For scaling in plots
library(dplyr)         # For data manipulation
library(cowplot)       # For combining ggplots
library(readr)         # For reading CSV files
library(stats)         # For statistical functions
library(tidyr)         # For data tidying
library(purrr)         # For functional programming
library(broom)         # For tidying model outputs

# Load data
rsim_df <- read.csv("~/Desktop/Robs_Rexp_Rratio_summary.csv")    # Rsim obs/exp ratio data
immune_cells = readRDS("~/Desktop/ConsensusTME_results_geneCounts.rds")    # ConsensusTME immune cell scores
meta_data = readRDS("~/Desktop/patient_info_from_uuid_SERVER_ALL_tcga.rds") # TCGA sample metadata
barcode_to_cancer <- setNames(meta_data$cancer_type, meta_data$submitter_id)  # Create lookup for cancer types

# Convert immune_cells list to long-format data frame
immune_df <- bind_rows(
  lapply(names(immune_cells), function(cancer) {
    mat <- immune_cells[[cancer]]
    df <- as.data.frame(t(mat))
    df$sample_id <- rownames(df)
    df$cancer_type <- cancer
    df
  }),
  .id = "cancer_index"
)

# Reshape immune data to long format
immune_long <- immune_df %>%
  pivot_longer(cols = -c(sample_id, cancer_type, cancer_index),
               names_to = "cell_type", values_to = "score")

# Define correlation function (cor + p-value)
calc_correlation <- function(df, variable) {
  df %>%
    group_by(cancer_type, cell_type) %>%
    summarise(
      cor = if (sum(complete.cases(.data[[variable]], score)) >= 2) {
        cor(.data[[variable]], score, use = "complete.obs")
      } else { NA_real_ },
      
      pval = if (sum(complete.cases(.data[[variable]], score)) >= 2) {
        cor.test(.data[[variable]], score, method = "pearson")$p.value
      } else { NA_real_ },
      
      .groups = "drop"
    )
}

# Merge immune and Rsim data
plot_data_rsim <- immune_long %>%
  inner_join(rsim_df, by = c("sample_id", "cancer_type"))

# Calculate correlations for MHC-I Rratio only
cor_rsim1 <- calc_correlation(plot_data_rsim, "Rratio_MHC1") %>%
  mutate(score_type = "MHC I")

# Now apply FDR correction only on MHC-I
cor_rsim1 <- cor_rsim1 %>%
  mutate(adj.P = p.adjust(pval, method = "fdr")) %>%
  filter(!is.na(adj.P)) %>%
  mutate(size_cut = cut(
    adj.P,
    breaks = c(-Inf, 1e-3, 1e-2, 5e-2, Inf),
    labels = c("≤0.001", "0.001–0.01", "0.01–0.05", ">0.05")
  )) %>%
  mutate(logP = -log10(adj.P))

# This is your final data to plot
cor_combined_rsim <- cor_rsim1

# Define facet label
facet_labels <- c("MHC I" = "MHC-I")
color_palette <- c("#053061", "#2166ac", "#d1e5f0", "white", "#fddbc7", "#b2182b", "#67001f")

# Plot
p_rsim <- ggplot(cor_combined_rsim, aes(x = cancer_type, y = cell_type, fill = cor, size = adj.P)) +
  geom_point(shape = 21, color = "grey30", stroke = 0.2, alpha = 0.92) +
  facet_wrap(~score_type, ncol = 1, scales = "free_x", labeller = as_labeller(facet_labels)) +
  scale_fill_gradientn(
    colors = color_palette,
    limits = c(-0.4, 0.4),
    oob = scales::squish,
    name = "Pearson r"
  ) + guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  )+scale_size_continuous(
    range = c(2,15),
    breaks = c(0.25, 0.3, 0.4, 0.5, 0.7, 0.9,1),
    labels = c("0.25", "0.3", "0.4", "0.5", "0.7", "0.9","1"),
    trans = "reverse",
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
  ) +
  labs(
    title = "Immune Cells Correlation with Rsim",
    x = "Cancer Type", y = "Immune Cell Type"
  )

# Save
ggsave("DotPlot_Rsim_Immune_CancerTypes_MHC1_only_correctFDR.png", plot = p_rsim, width = 20, height = 17, dpi = 300)
ggsave("DotPlot_Rsim_Immune_CancerTypes_MHC1_only_correctFDR.pdf", plot = p_rsim, width = 15, height = 11)

# Select immune score of interest
immune_score_name <- "Immune_Score"

# Filter for ImmuneScore in cor_combined_rsim → get r and adj.P
cor_labels <- cor_combined_rsim %>%
  filter(cell_type == immune_score_name) %>%
  mutate(
    label = paste0("r = ", round(cor, 2), "\nFDR = ", signif(adj.P, 2), "\np-value = ", signif(pval, 2))
  ) %>%
  select(cancer_type, score_type, label)

# Prepare plot data for CHOL MHC-I
plot_chol_mhci <- plot_data_rsim %>%
  filter(cancer_type == "CHOL", cell_type == immune_score_name) %>%
  mutate(Group = "CHOL")

# Prepare plot data for KIRP MHC-I
plot_kirp_mhci <- plot_data_rsim %>%
  filter(cancer_type == "KIRP", cell_type == immune_score_name) %>%
  mutate(Group = "KIRP")

# Prepare plot data for STAD MHC-I
plot_stad_mhci <- plot_data_rsim %>%
  filter(cancer_type == "STAD", cell_type == immune_score_name) %>%
  mutate(Group = "STAD")

# Combine all for side by side plot
plot_combined <- bind_rows(plot_chol_mhci, plot_kirp_mhci, plot_stad_mhci)

# Merge labels (for each Group)
cor_labels <- cor_labels %>%
  mutate(Group = case_when(
    score_type == "MHC I" & cancer_type == "CHOL" ~ "CHOL",
    score_type == "MHC I" & cancer_type == "KIRP" ~ "KIRP",
    score_type == "MHC I" & cancer_type == "STAD" ~ "STAD"
  )) %>%
  filter(!is.na(Group))

# Final scatter plot
ggplot(plot_combined, aes(x = Rratio_MHC1, y = score, color = Group)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 1.3) +
  facet_wrap(~Group, ncol = 3) +
  geom_text(data = cor_labels,
            aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = -0.1, vjust = 1.2, size = 6, fontface = "bold") +
  scale_color_manual(values = c(
    "CHOL" = "#1f78b4",
    "KIRP" = "#33a02c",
    "STAD" = "#e31a1c"
  )) +
  theme_bw(base_size = 16) +
  labs(
    title = "Immune Score vs Rsim Correlation (CHOL, KIRP, STAD) — MHC-I",
    x = "Rsim Ratio (MHC-I)",
    y = "Immune Score",
    color = ""
  ) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# Save
ggsave("Rsim_Immune_Score_Correlation_CHOL_KIRP_STAD_MHCI_no_MHCI_in_labels.png", width = 15, height = 8, dpi = 300)


# -----------------------
# Rsim boxplot preparation — MHC-I only
# -----------------------

# Reshape Rsim data to long format
rsim_long <- rsim_df %>%
  select(sample_id, cancer_type, Rratio_MHC1) %>%   # Only Rratio_MHC1
  pivot_longer(
    cols = starts_with("Rratio_MHC"),
    names_to = "score_type",
    values_to = "rsim"
  ) %>%
  mutate(score_type = recode(score_type, Rratio_MHC1 = "MHC-I")) %>%
  filter(!is.na(rsim), rsim != 0)

# Outlier filtering (IQR-based cutoff)
rsim_long <- rsim_long %>%
  group_by(score_type, cancer_type) %>%
  filter(rsim <= quantile(rsim, 0.75, na.rm = TRUE) + 1.0 * IQR(rsim, na.rm = TRUE)) %>%
  ungroup()

# Order cancer types by median Rsim
rsim_long <- rsim_long %>%
  group_by(score_type, cancer_type) %>%
  summarise(median_rsim = median(rsim, na.rm = TRUE), .groups = "drop") %>%
  group_by(score_type) %>%
  mutate(cancer_type_ordered = fct_reorder(cancer_type, median_rsim)) %>%
  right_join(rsim_long, by = c("score_type", "cancer_type")) %>%
  ungroup()

# Wilcoxon test vs. 1 + add significance stars
pvals_annotate <- rsim_long %>%
  group_by(score_type, cancer_type) %>%
  summarise(
    pval = tryCatch(suppressWarnings(wilcox.test(rsim, mu = 1)$p.value), error = function(e) NA),
    end_x = quantile(rsim, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    adj_p = p.adjust(pval, method = "fdr"),
    stars = case_when(
      adj_p < 0.001 ~ "***",
      adj_p < 0.01 ~ "**",
      adj_p < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(stars)) %>%
  left_join(
    rsim_long %>% select(score_type, cancer_type, cancer_type_ordered) %>% distinct(),
    by = c("score_type", "cancer_type")
  )

# Plot Rsim boxplot with significance stars — MHC-I only
main_plot = ggplot(rsim_long, aes(y = cancer_type_ordered, x = rsim, fill = cancer_type)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, width = 0.6, color = "black", coef = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray30") +
  geom_text(data = pvals_annotate,
            aes(x = end_x + 0.01, y = cancer_type_ordered, label = stars),
            inherit.aes = FALSE,
            size = 4.5, hjust = 0) +
  labs(
    x = "Rsim (obs/exp)",
    y = NULL,
    title = "Observed/Expected Mutation Ratio (Rsim)\nin MHC-I Region"
  )+
  facet_wrap(~ score_type, ncol = 1, scales = "fixed") +   # ncol = 1 → only MHC-I
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 2.5))

# Final plot with legend note
final_plot <- ggdraw(main_plot) +
  draw_text("* FDR < 0.05\n** FDR < 0.01\n*** FDR < 0.001",
            x = 0.90, y = 0.10,
            hjust = 1, size = 12,
            fontface = "italic")

# Save Rsim boxplot — MHC-I only
ggsave("Rratio_MHC1_boxplot.png", plot = final_plot, width = 5.5, height = 10, dpi = 300)

