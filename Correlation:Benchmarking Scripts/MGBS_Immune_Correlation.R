#MGBS and Immune Cells Correlation Across Cancer Types

library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(grid)




# Load data
cancer_results <- readRDS("~/Desktop/ConsensusTME_results_geneCounts.rds")
MGBS <- readRDS("~/Desktop/MGBS_TCGA.rds")
meta_data = readRDS("~/Desktop/patient_info_from_uuid_SERVER_ALL_tcga.rds") # TCGA sample metadata
barcode_to_cancer <- setNames(meta_data$cancer_type, meta_data$submitter_id)  # Create lookup for cancer types


MGBS$patient <- gsub("\\.[0-9]+$", "", MGBS$patient)


# Initialize result container
cor_long <- data.frame()

# Loop over cancer types
for (cancer in names(cancer_results)) {
  
  immune_mat <- cancer_results[[cancer]]
  patients <- colnames(immune_mat)
  
  # Match patients to MGBS dataframe
  idx <- match(patients, MGBS$patient)
  has_mgbs <- !is.na(idx)
  
  cat("Cancer:", cancer, "- Patients with MGBS:", sum(has_mgbs), "\n")
  
  # For each MGBS column (MGBS1, MGBS2, MGBSd)
  for (mgbs_col in c("MGBS1", "MGBS2", "MGBSd")) {
    
    if (sum(has_mgbs) > 2) {
      for (cell in rownames(immune_mat)) {
        
        v_immune <- as.numeric(immune_mat[cell, has_mgbs])
        v_mgbs <- as.numeric(MGBS[[mgbs_col]][idx[has_mgbs]])
        
        # Debug print
        cat("  ", mgbs_col, "-", cell, "- Immune len:", length(v_immune), 
            " MGBS len:", length(v_mgbs), 
            " sd immune:", sd(v_immune), 
            " sd mgbs:", sd(v_mgbs), "\n")
        
        if (length(v_immune) > 2 && sd(v_immune) > 0 && sd(v_mgbs) > 0) {
          test <- cor.test(v_immune, v_mgbs, method = "pearson")
          
          cor_long <- rbind(cor_long, data.frame(
            Cancer = cancer,
            Cell = cell,
            ScoreType = mgbs_col,
            Cor = test$estimate,
            P = test$p.value
          ))
        }
      }
    } else {
      cat("  Skipped", mgbs_col, "- Not enough patients with MGBS\n")
    }
  }
}



# Adjust p-values
cor_long$adj.P <- p.adjust(cor_long$P, method = "fdr")

cor_long <- cor_long %>% filter(ScoreType != "MGBSd")


# Make sure p_value or P is numeric — use the correct column name!
cor_long$logP <- -log10(cor_long$adj.P)

# Clip Pearson r scale to -0.4 to 0.4
color_palette <- c("#053061", "#2166ac", "#d1e5f0", "white", "#fddbc7", "#b2182b", "#67001f")

# Optional: You can define color_breaks if you want
color_breaks <- c(-0.4, -0.2, 0, 0.2, 0.4)

# Compute min and max logP
logP_range <- range(cor_long$logP, na.rm = TRUE)

# Define reasonable breaks within that range
logP_breaks <- pretty(logP_range, n = 5)

# Define dot size range — adjust if needed
dot_size_range <- c(0, 15)
facet_labels <- c("MGBS1" = "MGBS-I", "MGBS2" = "MGBS-II")
# Plot
p <- ggplot(cor_long, aes(x = Cancer, y = Cell, fill = Cor, size = logP)) +
  geom_point(shape = 21, color = "grey30", stroke = 0.2, alpha = 0.92) +
  facet_wrap(~ScoreType, ncol = 2, scales = "free_x",
             labeller = as_labeller(facet_labels))+
  scale_fill_gradientn(
    colors = color_palette,
    limits = c(-0.4, 0.4),
    oob = scales::squish,
    name = "Pearson r"
  ) +
  scale_size_continuous(
    range = dot_size_range,
    breaks = logP_breaks,
    labels = scales::number_format(accuracy = 0.1)(10^(-logP_breaks)),
    name = "FDR"
  )+
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
    panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),  # no deprecated warning
    legend.key.height = unit(1.5, "cm"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    strip.background = element_rect(fill = "grey96", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Immune Cells Correlation Across Cancer Types",
    x = "Cancer Type", y = "Immune Cell Type"
  )

# Save
ggsave("DotPlot_Immune_CancerTypes_adjP_updated.png", plot = p, width = 30, height = 17, dpi = 300)
ggsave("DotPlot_Immune_CancerTypes_adjP_updated.pdf", plot = p, width = 23, height = 11)



# Extract all immune signatures
immune_sigs <- rownames(cancer_results[[1]])

# Create summary matrix
summary_mat <- matrix(NA, nrow = length(immune_sigs), ncol = length(cancer_results))
rownames(summary_mat) <- immune_sigs
colnames(summary_mat) <- names(cancer_results)

for (cancer in names(cancer_results)) {
  mat <- cancer_results[[cancer]]
  common_sigs <- intersect(rownames(mat), immune_sigs)
  summary_mat[common_sigs, cancer] <- apply(mat[common_sigs, , drop=FALSE], 1, median, na.rm=TRUE)
}

# Scale the matrix row-wise (z-score)
summary_mat_scaled <- t(scale(t(summary_mat)))
rownames(summary_mat_scaled) <- gsub("_", " ", rownames(summary_mat_scaled))
# Step 1: Get patient counts per cancer type
patient_counts <- sapply(cancer_results, ncol)  # number of samples per cancer

# Step 2: Rename columns in summary matrix
colnames(summary_mat_scaled) <- paste0(
  colnames(summary_mat_scaled), " (n=", patient_counts[colnames(summary_mat_scaled)], ")"
)

# Create the heatmap
pheatmap(
  summary_mat_scaled,
  color = colorRampPalette(c("#053061", "#2166ac", "#4393c3", "#92c5de", 
                             "#d1e5f0", "#fddbc7", "#f4a582", "#d6604d", 
                             "#b2182b", "#67001f"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Immune Cell Landscape Across TCGA Cancers",
  width = 10,
  height = 8,
  filename = "ImmuneDeconv_Heatmap_FromRDS.png"
)

immune_df <- bind_rows(
  lapply(names(cancer_results), function(cancer) {
    mat <- cancer_results[[cancer]]
    df <- as.data.frame(t(mat))
    df$sample_id <- rownames(df)
    df$cancer_type <- cancer
    df
  }),
  .id = "cancer_index"
)

immune_score_df <- immune_df %>%
  select(-cancer_index) %>%
  rowwise() %>%
  mutate(immune_score = Immune_Score) %>%
  ungroup() %>%
  select(sample_id, cancer_type, immune_score) %>%
  distinct()


# Prepare cor_labels from cor_long (Immune_Score only)
cor_labels <- cor_long %>%
  filter(Cell == "Immune_Score") %>%
  filter(Cancer %in% c("CHOL", "KICH", "UCS", "MESO")) %>%
  mutate(Type = recode(ScoreType, MGBS1 = "MGBS-I", MGBS2 = "MGBS-II")) %>%
  mutate(label = paste0(
    "r = ", round(Cor, 2), "\nFDR = ", signif(adj.P, 2), "\np = ", signif(P, 2)
  )) %>%
  select(Cancer, Type, label)
# make column name match the facet!

# Add custom x/y positions for label to avoid Inf problems
cor_labels <- cor_labels %>%
  mutate(x = 0.4,  # slight margin from right
         y = case_when(
           Type == "MGBS-I" ~ 0.03,  # adjust this if needed based on your Y scale for MGBS-I
           Type == "MGBS-II" ~ 0.25  # adjust this if needed based on your Y scale for MGBS-II
         ))

cor_labels <- cor_labels %>%
  rename(cancer_type = Cancer)


# Prepare plot data
plot_df_long <- MGBS %>%
  mutate(cancer_type = barcode_to_cancer[patient]) %>%
  select(patient, cancer_type, MGBS1, MGBS2) %>%
  inner_join(immune_score_df, by = c("patient" = "sample_id", "cancer_type")) %>%
  pivot_longer(cols = c(MGBS1, MGBS2), names_to = "Type", values_to = "Score") %>%
  mutate(Type = recode(Type, MGBS1 = "MGBS-I", MGBS2 = "MGBS-II")) %>%
  filter(cancer_type %in% c("CHOL", "KICH", "UCS", "MESO"))

# Split plot data
plot_df_MGBS1 <- plot_df_long %>% filter(Type == "MGBS-I")
plot_df_MGBS2 <- plot_df_long %>% filter(Type == "MGBS-II")

# Color palette
color_palette <- c(
  "CHOL" = "#1f78b4",
  "KICH" = "#33a02c",
  "UCS"  = "#e31a1c",
  "MESO" = "#ff7f00"
)

# Plot MGBS-I
p1 <- ggplot(plot_df_MGBS1, aes(x = immune_score, y = Score, color = cancer_type)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 1.2) +
  facet_wrap(~ cancer_type, scales = "free_y") +
  geom_text(
    data = cor_labels %>% filter(Type == "MGBS-I"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 5,
    color = "black",
    fontface = "bold"
  ) +
  scale_color_manual(values = color_palette) +
  labs(title = "MGBS-I", x = "Immune Score", y = "MGBS Score") +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing = unit(1.2, "lines")
  )

# Plot MGBS-II
p2 <- ggplot(plot_df_MGBS2, aes(x = immune_score, y = Score, color = cancer_type)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 1.2) +
  facet_wrap(~ cancer_type, scales = "free_y") +
  geom_text(
    data = cor_labels %>% filter(Type == "MGBS-II"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 5,
    color = "black",
    fontface = "bold"
  ) +
  scale_color_manual(values = color_palette) +
  labs(title = "MGBS-II", x = "Immune Score", y = "MGBS Score") +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing = unit(1.2, "lines")
  )

# Combine p1 and p2
library(ggpubr)

combined_plot <- ggarrange(p1, p2, ncol = 2, widths = c(1, 1), common.legend = TRUE, legend = "top") %>%
  annotate_figure(
    top = text_grob("Immune Score vs MGBS-I and MGBS-II Correlation in CHOL, KICH, UCS, MESO", face = "bold", size = 20)
  )

# Save
ggsave("ImmuneScore_vs_MGBS_Correlation_CHOL_KICH_UCS_MESO_PRETTY_correct_labels_FINAL.png", combined_plot, width = 16, height = 8, dpi = 300)
