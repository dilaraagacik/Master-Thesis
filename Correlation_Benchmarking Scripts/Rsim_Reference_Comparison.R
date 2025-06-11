# ---- 1. Libraries ----
library(dplyr)
library(ggplot2)
library(readxl)
library(readr)

# ---- 2. Load Data ----
rsim <- read_csv("~/Desktop/Robs_Rexp_Rratio_summaryy.csv")
sheet_names <- excel_sheets("~/Desktop/Rsim/rsim.xlsx")
all_sheets <- lapply(sheet_names, function(sheet) {
  read_excel("~/Desktop/Rsim/rsim.xlsx", sheet = sheet, skip = 1)
})
names(all_sheets) <- sheet_names
sheet2 <- all_sheets[[2]]
names(sheet2)[1] <- "sample_id"
sheet2$sample_id <- as.character(sheet2$sample_id)
rsim$sample_id <- as.character(rsim$sample_id)

# ---- 3. Merge Results ----
comparison <- right_join(sheet2, rsim, by = "sample_id")
comparison$jimmy_ratio <- as.numeric(comparison$`dNHLA/dNnonHLA...10`)

# ---- 4. Cancer Type Annotation ----
patient_info <- readRDS("~/Desktop/patient_info_from_uuid_SERVER_ALL_tcga.rds")
patient_info <- patient_info %>%
  rename(sample_id = submitter_id)

comparison <- comparison %>%
  left_join(patient_info, by = "sample_id")

# Perform your joins as before
comparison <- right_join(sheet2, rsim, by = "sample_id")
comparison$jimmy_ratio <- as.numeric(comparison$`dNHLA/dNnonHLA...10`)

# Join with patient_info
comparison <- left_join(comparison, patient_info, by = "sample_id")



# Remove NAs (if any)
plot_df <- comparison %>%
  filter(!is.na(jimmy_ratio), !is.na(Rratio_MHC1), !is.na(cancer_type.y))
plot_df <- plot_df %>%
  filter(is.finite(jimmy_ratio), is.finite(Rratio_MHC1))

plot_df$cancer_type.y <- factor(plot_df$cancer_type.y)

# 1. Balanced sampling
plot_df_clean <- plot_df %>%
  filter(is.finite(jimmy_ratio), is.finite(Rratio_MHC1)) %>%
  filter(jimmy_ratio > 0, Rratio_MHC1 > 0)

patients_per_cancer <- 10
set.seed(42)
plot_df_sampled <- plot_df_clean %>%
  group_by(cancer_type.y) %>%
  sample_n(size = min(n(), patients_per_cancer)) %>%
  ungroup()

# 2. Recalculate axis limits
lims_sampled <- range(c(plot_df_sampled$jimmy_ratio, plot_df_sampled$Rratio_MHC1), na.rm = TRUE)
lims_sampled <- lims_sampled + c(-0.05, 0.05) * diff(lims_sampled)

# 3. Correlation
cor_test_sampled <- cor.test(plot_df_sampled$jimmy_ratio, plot_df_sampled$Rratio_MHC1, use = "complete.obs")
r_val_sampled <- cor_test_sampled$estimate
p_val_sampled <- cor_test_sampled$p.value

cor_label_sampled <- bquote(italic(r) == .(format(r_val_sampled, digits = 2)) ~ ", " ~ italic(P) == .(ifelse(p_val_sampled < 0.05, "<0.05", format.pval(p_val_sampled, digits = 2, eps = 0.01, scientific = FALSE))))

# 4. Plot
ggplot(plot_df_sampled, aes(x = jimmy_ratio, y = Rratio_MHC1, fill = cancer_type.y)) +
  geom_point(shape = 21, color = "black", size = 4, alpha = 0.92, stroke = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1.3) +
  scale_x_continuous(limits = lims_sampled, name = "Reference Rsim Score (Van den Eynden et al. 2019)") +
  scale_y_continuous(limits = lims_sampled, name = "Observed Rsim Score") +
  scale_fill_viridis_d(option = "turbo", name = "Cancer Type") +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray92"),
    aspect.ratio = 1,
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Reference vs Observed Rsim Score by Cancer Type (Sampled Patients)",
    caption = "Data: Nature Genetics (2019), https://doi.org/10.1038/s41588-019-0532-6.\nSampled max 10 patients per cancer type."
  ) +
  annotate("text",
           x = min(lims_sampled) + 0.04 * diff(lims_sampled),
           y = max(lims_sampled) - 0.08 * diff(lims_sampled),
           label = as.expression(cor_label_sampled),
           hjust = 0, vjust = 1, size = 6.2, color = "#222222")

# 5. Save
ggsave("rratio_vs_reference_cancer_Rsim_colored_thesis_sampled_patients.png", width = 15, height = 10, dpi = 300)



# 1. Balanced sampling and TMB calculation
plot_df_clean <- plot_df %>%
  filter(is.finite(jimmy_ratio), is.finite(Rratio_MHC1)) %>%
  filter(jimmy_ratio > 0, Rratio_MHC1 > 0)

patients_per_cancer <- 10
set.seed(42)
plot_df_sampled_tmb <- plot_df_clean %>%
  group_by(cancer_type.y) %>%
  sample_n(size = min(n(), patients_per_cancer)) %>%
  ungroup() %>%
  mutate(TMB = n_nonbinding_MHC1 + n_binding_MHC1) # Calculate TMB

# 2. Recalculate axis limits
lims_sampled <- range(c(plot_df_sampled_tmb$jimmy_ratio, plot_df_sampled_tmb$Rratio_MHC1), na.rm = TRUE)
lims_sampled <- lims_sampled + c(-0.05, 0.05) * diff(lims_sampled)

# 3. Correlation
cor_test_sampled <- cor.test(plot_df_sampled_tmb$jimmy_ratio, plot_df_sampled_tmb$Rratio_MHC1, use = "complete.obs")
r_val_sampled <- cor_test_sampled$estimate
p_val_sampled <- cor_test_sampled$p.value

cor_label_sampled <- bquote(italic(r) == .(format(r_val_sampled, digits = 2)) ~ ", " ~ italic(P) == .(ifelse(p_val_sampled < 0.05, "<0.05", format.pval(p_val_sampled, digits = 2, eps = 0.01, scientific = FALSE))))

ggplot(plot_df_sampled_tmb, aes(x = jimmy_ratio, y = Rratio_MHC1, size = log10(TMB)))+
  geom_point(shape = 21, color = "black", fill = "#4682B4", alpha = 0.8, stroke = 0.7) + # SteelBlue fill
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1.3) +
  scale_x_continuous(limits = lims_sampled, name = "Reference Rsim Score (Van den Eynden et al. 2019)") +
  scale_y_continuous(limits = lims_sampled, name = "Observed Rsim Score") +
  scale_size_continuous(
    name = "Tumor Mutation\nBurden",
    range = c(1, 14),
    breaks = log10(c(3,5,10, 30, 60, 150, 300)),
    labels = c("3","5","10", "30", "60", "150", "300+")
  )+coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray92"),
    aspect.ratio = 1,
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Reference vs Observed Rsim Score",
    subtitle = "Point size corresponds to Tumor Mutation Burden (TMB)",
    caption = "Data: Nature Genetics (2019), https://doi.org/10.1038/s41588-019-0532-6.\nSampled max 10 patients per cancer type."
  ) +
  annotate("text",
           x = min(lims_sampled) + 0.04 * diff(lims_sampled),
           y = max(lims_sampled) - 0.08 * diff(lims_sampled),
           label = cor_label_sampled,  # CORRECT
           hjust = 0, vjust = 1, size = 6.2, color = "#222222")


# 5. Save the new plot
ggsave("rratio_vs_reference_TMB_Rsim_sampled_patients.png", width = 15, height = 10, dpi = 300)




