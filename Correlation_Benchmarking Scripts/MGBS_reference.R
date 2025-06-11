#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(patchwork)

# ----------------------------
# Load immune study RData
# IMPORTANT: load() returns object names — so we capture the name
immune_loaded <- load("~/Desktop/CCGGlab-mhc_immunotherapy-ee121bc/data/MHC_immunotherapy.RData")
cat("Loaded immune object(s):", immune_loaded, "\n")

# We assume the object is called ICB_study (as you used later)
# If not, check: print(ls()) or print(immune_loaded)

# ----------------------------
# Load your computed MGBS scores
my_results <- readRDS("~/Desktop/MGBS_scores.rds")

# ----------------------------
# Load ICB study dataframe (comes from MHC_immunotherapy.RData)
# Here you used:
icb_study <- ICB_study  # assuming ICB_study is in your RData

# ----------------------------
# Match patients — inner join
merged_df <- inner_join(
  my_results %>% select(patient, MGBS1_my = MGBS1, MGBS2_my = MGBS2),
  icb_study %>% select(patient, MGBS1_icb = MGBS1, MGBS2_icb = MGBS2),
  by = "patient"
)
cat("Merged", nrow(merged_df), "patients\n")

# ----------------------------
# Prepare for plotting — long format
# We pivot to long format, then wide again to have columns: icb, my
plot_df <- merged_df %>%
  pivot_longer(cols = c(MGBS1_my, MGBS2_my, MGBS1_icb, MGBS2_icb),
               names_to = c("score", "source"),
               names_pattern = "(MGBS[12])_(.*)",
               values_to = "value") %>%
  pivot_wider(names_from = source, values_from = value) %>%
  rename(
    my = my,
    icb = icb
  )

# ----------------------------
# Plot function for one score
plot_one <- function(score_name) {
  df <- plot_df %>% filter(score == score_name)
  
  # Correlation and p-value
  r_val <- cor(df$my, df$icb, use = "complete.obs")
  p_val <- cor.test(df$my, df$icb)$p.value
  
  # Create plot
  ggplot(df, aes(x = icb, y = my)) +
    geom_point(color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue", linewidth = 0.8) +
    labs(
      title = score_name,
      x = paste0("Claeys, A (2024) ", score_name),
      y = paste0("Observed ", score_name)
    ) +
    annotate("text",
             x = min(df$icb, na.rm = TRUE),
             y = max(df$my, na.rm = TRUE),
             hjust = 0, vjust = 1,
             label = sprintf("r = %.2f, p = %s",
                             r_val,
                             ifelse(p_val < 0.001, "<0.001", formatC(p_val, format = "f", digits = 3))),
             size = 5) +
    theme_bw(base_size = 14)
}

# ----------------------------
# Plot MGBS1 and MGBS2 side by side
p1 <- plot_one("MGBS1")
p2 <- plot_one("MGBS2")

# Combine plots
final_plot <- p1 + p2 + plot_annotation(
  caption = "Adapted from: Claeys & Van den Eynden, Commun Med 4, 184 (2024)"
)

# ----------------------------
# Save figure
ggsave("MGBS1_2vsClaeys.png", final_plot, width = 10, height = 5, dpi = 300)

# Show figure
print(final_plot)
