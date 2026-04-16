# =============================================================================
# Myopathy Diagnostic Journey: Main Analysis
# =============================================================================
# Publication: "The diagnostic journey of patients being investigated for
#              myopathy in a tertiary centre in England"
#              Journal of Neurology (2025) 272:35
#              https://doi.org/10.1007/s00415-024-12737-y
#
# Description: This script reproduces the primary statistical analyses and
#              visualisations reported in the paper, covering patient
#              demographics, investigation utilisation, diagnostic timing,
#              and test combination performance.
#
# Input:  data/myopathy_main.csv   — main patient-level dataset
#         data/myopathy_intervals.csv — derived interval dataset (optional;
#                                       can be built from main dataset)
#
# Note:   Data are available upon request from the corresponding author.
#         Column naming conventions expected in the input CSV are documented
#         in the Data Dictionary section below.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries
# -----------------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ComplexUpset)
library(caret)
library(dunn.test)

# -----------------------------------------------------------------------------
# 1. Data Dictionary
# -----------------------------------------------------------------------------
# Expected columns in myopathy_main.csv:
#
#   study_id              — anonymised patient identifier
#   age                   — age at first hospital visit (years)
#   gender                — "Male" / "Female"
#   specialist            — referring specialty category
#   final_diagnosis       — "Myopathy" / "Non-myopathy"
#   myopathy_subtype      — specific myopathy subtype (if applicable)
#   myopathy_mimic        — alternative diagnosis (if non-myopathy)
#   ck                    — CK result (0 = not done, 1 = non-myopathy, 2/3 = myopathy)
#   myositis_abs          — myositis antibody result (same coding)
#   neuroaxis_mri         — neuroaxis MRI result (0 = not done, 1 = performed)
#   muscle_mri            — muscle MRI result (same coding as ck)
#   emg                   — EMG result (same coding as ck)
#   muscle_biopsy         — muscle biopsy result (same coding as ck)
#   genetics              — genetics result (same coding as ck)
#   interval_blood_wks    — weeks from referral to blood tests
#   interval_emg_wks      — weeks from referral to EMG
#   interval_neuroaxis_wks— weeks from referral to neuroaxis MRI
#   interval_muscle_mri_wks — weeks from referral to muscle MRI
#   interval_biopsy_wks   — weeks from referral to biopsy
#   interval_genetics_wks — weeks from referral to genetics
#   interval_dx_referral_wks — weeks from referral to diagnosis
#   interval_dx_symptoms_wks — weeks from symptom onset to diagnosis
#   symptom_duration_wks  — weeks from symptom onset to first hospital visit

# -----------------------------------------------------------------------------
# 2. Load Data
# -----------------------------------------------------------------------------
# Replace the path below with your local file path or working directory path.
data_raw <- read.csv("data/myopathy_main.csv", stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 3. Cohort Overview — Demographics
# -----------------------------------------------------------------------------

# 3.1 Cohort size by final diagnosis
data_raw %>% count(final_diagnosis)

# 3.2 Age: whole cohort and by diagnosis group
summary(data_raw$age)
summary(data_raw$age[data_raw$final_diagnosis == "Myopathy"])
summary(data_raw$age[data_raw$final_diagnosis == "Non-myopathy"])

# 3.3 Age comparison — Wilcoxon rank-sum test (Bonferroni-corrected)
wilcox_age_p <- wilcox.test(age ~ final_diagnosis, data = data_raw)$p.value
p.adjust(wilcox_age_p, method = "bonferroni")

# 3.4 Gender distribution and chi-squared test
gender_by_dx <- pivot_wider(
  data_raw %>% count(gender, final_diagnosis),
  names_from  = final_diagnosis,
  values_from = n
) %>% as.data.frame()

gender_by_dx
data_raw %>% count(gender)

chisq_gender <- chisq.test(gender_by_dx[, c("Myopathy", "Non-myopathy")])
p.adjust(chisq_gender$p.value, method = "bonferroni")

# 3.5 Referring specialty distribution and chi-squared test
specialist_by_dx <- pivot_wider(
  data_raw %>% count(specialist, final_diagnosis),
  names_from  = final_diagnosis,
  values_from = n
) %>% as.data.frame()

specialist_by_dx
data_raw %>% count(specialist)

chisq_specialist <- chisq.test(specialist_by_dx[, c("Myopathy", "Non-myopathy")])
p.adjust(chisq_specialist$p.value, method = "bonferroni")

# -----------------------------------------------------------------------------
# 4. Diagnosis Breakdown
# -----------------------------------------------------------------------------

# 4.1 Myopathy subtypes (Table 2)
data_raw %>%
  filter(final_diagnosis == "Myopathy") %>%
  group_by(myopathy_subtype) %>%
  summarise(
    count      = n(),
    percentage = count / nrow(filter(data_raw, final_diagnosis == "Myopathy")) * 100
  ) %>%
  arrange(desc(count))

# 4.2 Non-myopathy diagnoses (Table 2)
data_raw %>%
  filter(final_diagnosis == "Non-myopathy") %>%
  group_by(myopathy_mimic) %>%
  summarise(
    count      = n(),
    percentage = count / nrow(filter(data_raw, final_diagnosis == "Non-myopathy")) * 100
  ) %>%
  arrange(desc(count))

# -----------------------------------------------------------------------------
# 5. Investigation Utilisation
# -----------------------------------------------------------------------------

# 5.1 Build binary investigation dataset
#     Any result code ∈ values_positive indicates the test was performed.
investigation_cols <- c("ck", "myositis_abs", "neuroaxis_mri",
                        "muscle_mri", "emg", "muscle_biopsy", "genetics")

values_positive <- list(
  ck           = c(1, 2),
  myositis_abs = c(1, 2),
  neuroaxis_mri = c(1),
  muscle_mri   = c(1, 2),
  emg          = c(1, 2),
  muscle_biopsy = c(1, 2),
  genetics     = c(1, 2)
)

investigation_df <- data_raw %>%
  select(study_id, all_of(investigation_cols), final_diagnosis,
         specialist, myopathy_subtype, myopathy_mimic) %>%
  mutate(across(
    all_of(investigation_cols),
    ~ ifelse(. %in% values_positive[[cur_column()]], 1L, 0L)
  ))

# 5.2 Utilisation counts by diagnosis group
count_utilisation <- function(col) {
  list(
    myopathy     = sum(investigation_df[investigation_df$final_diagnosis == "Myopathy",     col] == 1),
    non_myopathy = sum(investigation_df[investigation_df$final_diagnosis == "Non-myopathy", col] == 1)
  )
}
lapply(setNames(investigation_cols, investigation_cols), count_utilisation)

# 5.3 Upset plot — test combinations (Figure 1b/c)
# Myopathy group
ComplexUpset::upset(
  filter(investigation_df, final_diagnosis == "Myopathy"),
  investigation_cols,
  n_intersections = 10,
  width_ratio     = 0.4,
  height_ratio    = 1,
  stripes         = "white",
  base_annotations = list(
    "Intersection size" = intersection_size(counts = FALSE)
  ),
  name   = " ",
  matrix = intersection_matrix(
    geom    = geom_point(size = 5),
    segment = geom_segment(linewidth = 1)
  ),
  set_sizes = upset_set_size(
    position = "right",
    geom     = geom_bar(width = 0.85)
  ) + ylab("Investigation set size"),
  themes = upset_modify_themes(list(
    "intersections_matrix" = theme(text = element_text(size = 12)),
    "overall_sizes"        = theme(text = element_text(size = 12)),
    "Intersection size"    = theme(text = element_text(size = 12))
  ))
)

# Non-myopathy group
ComplexUpset::upset(
  filter(investigation_df, final_diagnosis == "Non-myopathy"),
  investigation_cols,
  n_intersections = 10,
  width_ratio     = 0.4,
  height_ratio    = 1,
  stripes         = "white",
  base_annotations = list(
    "Intersection size" = intersection_size(counts = FALSE)
  ),
  name   = " ",
  matrix = intersection_matrix(
    geom    = geom_point(size = 5),
    segment = geom_segment(linewidth = 1)
  ),
  set_sizes = upset_set_size(
    position = "right",
    geom     = geom_bar(width = 0.85)
  ) + ylab("Investigation set size"),
  themes = upset_modify_themes(list(
    "intersections_matrix" = theme(text = element_text(size = 12)),
    "overall_sizes"        = theme(text = element_text(size = 12)),
    "Intersection size"    = theme(text = element_text(size = 12))
  ))
)

# -----------------------------------------------------------------------------
# 6. Number of Tests per Patient (Figure 1d)
# -----------------------------------------------------------------------------
myopathy_investigations     <- filter(investigation_df, final_diagnosis == "Myopathy")     %>% select(all_of(investigation_cols))
non_myopathy_investigations <- filter(investigation_df, final_diagnosis == "Non-myopathy") %>% select(all_of(investigation_cols))

n_tests_myopathy     <- rowSums(myopathy_investigations)
n_tests_non_myopathy <- rowSums(non_myopathy_investigations)

# Descriptive statistics
cat(sprintf("Myopathy:     mean = %.2f, SD = %.2f\n", mean(n_tests_myopathy),     sd(n_tests_myopathy)))
cat(sprintf("Non-myopathy: mean = %.2f, SD = %.2f\n", mean(n_tests_non_myopathy), sd(n_tests_non_myopathy)))

# Welch t-test
t.test(n_tests_myopathy, n_tests_non_myopathy, alternative = "two.sided", var.equal = FALSE)

# Boxplot
df_n_tests <- data.frame(
  n_tests       = c(n_tests_myopathy, n_tests_non_myopathy),
  final_diagnosis = c(
    rep("Myopathy",     length(n_tests_myopathy)),
    rep("Non-myopathy", length(n_tests_non_myopathy))
  )
)

ggplot(df_n_tests, aes(x = final_diagnosis, y = n_tests, fill = final_diagnosis)) +
  geom_boxplot(outlier.shape = NA, lwd = 1.5) +
  labs(x = "", y = "Number of tests") +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 7, 1)) +
  scale_fill_manual(values = c("Myopathy" = "#FD6467", "Non-myopathy" = "#7294D4")) +
  theme_minimal() +
  theme(
    text                 = element_text(size = 14),
    axis.line.y.left     = element_line(linewidth = 1),
    axis.line.x.bottom   = element_line(linewidth = 1),
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.background     = element_blank(),
    axis.line            = element_line(colour = "black")
  ) +
  geom_bracket(
    xmin       = 1, xmax = 2,
    y.position = 7.5,
    label      = "****",
    tip.length = c(0.03, 0.49),
    inherit.aes = FALSE,
    size        = 1.5,
    label.size  = 6
  ) +
  guides(fill = "none")

# -----------------------------------------------------------------------------
# 7. Investigation Timing (Figure 2)
# -----------------------------------------------------------------------------

# 7.1 Build intervals dataset
#     Neuroaxis and muscle MRI intervals are set to NA when the test was not done.
interval_neuroaxis_mri <- ifelse(data_raw$neuroaxis_mri == 0, NA, data_raw$interval_neuroaxis_wks)
interval_muscle_mri    <- ifelse(data_raw$muscle_mri    == 0, NA, data_raw$interval_muscle_mri_wks)

intervals_df <- data.frame(
  study_id              = data_raw$study_id,
  interval_blood        = data_raw$interval_blood_wks,
  interval_emg          = as.numeric(data_raw$interval_emg_wks),
  interval_neuroaxis    = interval_neuroaxis_mri,
  interval_muscle_mri   = interval_muscle_mri,
  interval_biopsy       = data_raw$interval_biopsy_wks,
  interval_genetics     = as.numeric(data_raw$interval_genetics_wks),
  interval_dx_referral  = as.numeric(data_raw$interval_dx_referral_wks),
  interval_dx_symptoms  = as.numeric(data_raw$interval_dx_symptoms_wks),
  symptom_duration      = as.numeric(data_raw$symptom_duration_wks),
  final_diagnosis       = data_raw$final_diagnosis,
  myopathy_subtype      = data_raw$myopathy_subtype,
  myopathy_mimic        = data_raw$myopathy_mimic,
  specialist            = data_raw$specialist
)

# 7.2 Wilcoxon rank-sum tests for investigation timing by diagnosis
investigation_interval_cols <- c(
  "interval_blood", "interval_emg", "interval_neuroaxis",
  "interval_muscle_mri", "interval_biopsy", "interval_genetics"
)
lapply(investigation_interval_cols, function(col) {
  wilcox.test(as.formula(paste(col, "~ final_diagnosis")), data = intervals_df)
})

# 7.3 Investigation timing boxplot (Figure 2)
intervals_long <- pivot_longer(
  intervals_df,
  cols      = all_of(investigation_interval_cols),
  names_to  = "investigation",
  values_to = "interval_weeks"
)

ggplot(intervals_long,
  aes(
    x    = factor(investigation,
                  levels = c("interval_biopsy", "interval_muscle_mri",
                             "interval_neuroaxis", "interval_emg",
                             "interval_blood")),
    y    = interval_weeks,
    fill = final_diagnosis
  )) +
  geom_boxplot(outlier.shape = NA, lwd = 1,
               position = position_dodge(width = 0.85)) +
  labs(x = "", y = "Duration since referral (weeks)", fill = "Final conclusion") +
  scale_x_discrete(labels = c(
    "interval_blood"      = "Blood",
    "interval_emg"        = "EMG",
    "interval_neuroaxis"  = "MRI Neuroaxis",
    "interval_muscle_mri" = "MRI Muscle",
    "interval_biopsy"     = "Biopsy"
  )) +
  scale_y_continuous(limits = c(0, 104), breaks = seq(0, 104, 12)) +
  scale_fill_manual(values = c("Myopathy" = "#FD6467", "Non-myopathy" = "#7294D4")) +
  theme_minimal() +
  theme(
    text               = element_text(size = 14),
    axis.line.y.left   = element_line(linewidth = 1.5),
    axis.line.x.bottom = element_line(linewidth = 1.5),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    legend.position    = "top",
    legend.direction   = "horizontal"
  ) +
  coord_flip() +
  geom_bracket(
    xmin        = c(4.8125, 3.8125, 2.8125, 1.8125, 0.8125),
    xmax        = c(5.1875, 4.1875, 3.1875, 2.1875, 1.1875),
    y.position  = rep(100, 5),
    label       = c("n.s.", "**", "n.s.", "n.s.", "*"),
    tip.length  = 0.0005,
    inherit.aes = FALSE,
    size        = 1.5,
    label.size  = 6,
    coord.flip  = TRUE
  )

# -----------------------------------------------------------------------------
# 8. Time to Diagnosis (Figure 3)
# -----------------------------------------------------------------------------

# 8.1 Symptom onset to first hospital visit
wilcox.test(symptom_duration ~ final_diagnosis, data = intervals_df)

# 8.2 Symptom onset to diagnosis
wilcox.test(interval_dx_symptoms ~ final_diagnosis, data = intervals_df)
summary(intervals_df$interval_dx_symptoms[intervals_df$final_diagnosis == "Myopathy"])
summary(intervals_df$interval_dx_symptoms[intervals_df$final_diagnosis == "Non-myopathy"])

# 8.3 First hospital visit to diagnosis
wilcox.test(interval_dx_referral ~ final_diagnosis, data = intervals_df)
summary(intervals_df$interval_dx_referral[intervals_df$final_diagnosis == "Myopathy"])
summary(intervals_df$interval_dx_referral[intervals_df$final_diagnosis == "Non-myopathy"])

# 8.4 Kruskal-Wallis and Dunn post-hoc for myopathy subtypes — time to diagnosis
kruskal.test(interval_dx_referral ~ myopathy_subtype,
             data = filter(intervals_df, final_diagnosis == "Myopathy"))

dunn.test(
  x     = filter(intervals_df, final_diagnosis == "Myopathy")$interval_dx_referral,
  g     = filter(intervals_df, final_diagnosis == "Myopathy")$myopathy_subtype,
  method = "bonferroni"
)

# 8.5 Time to diagnosis by referring specialist
# By referral interval
lapply(unique(intervals_df$specialist), function(sp) {
  cat("Specialist:", sp, "\n")
  wilcox.test(
    interval_dx_referral ~ final_diagnosis,
    data = filter(intervals_df, specialist == sp)
  )
})

# By symptom onset interval
lapply(unique(intervals_df$specialist), function(sp) {
  cat("Specialist:", sp, "\n")
  wilcox.test(
    interval_dx_symptoms ~ final_diagnosis,
    data = filter(intervals_df, specialist == sp)
  )
})

# 8.6 Boxplot — time from referral to diagnosis by specialist (Supplementary Figure)
ggplot(intervals_df,
  aes(x = specialist, y = interval_dx_referral, fill = final_diagnosis)) +
  geom_boxplot(width = 0.95, outlier.shape = NA) +
  labs(
    x    = "",
    y    = str_wrap("Time from referral to diagnosis (weeks)", width = 35),
    fill = "Final conclusion"
  ) +
  scale_y_continuous(limits = c(0, 384), breaks = seq(0, 384, 48)) +
  scale_fill_manual(values = c("Myopathy" = "#FD6467", "Non-myopathy" = "#7294D4")) +
  theme_minimal() +
  theme(
    text               = element_text(size = 14),
    axis.line.y.left   = element_line(linewidth = 1.5),
    axis.line.x.bottom = element_line(linewidth = 1.5),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    legend.position    = "bottom"
  ) +
  geom_bracket(
    xmin        = c(0.7875, 1.7875, 2.7875, 3.7875),
    xmax        = c(1.2125, 2.2125, 3.2125, 4.2125),
    y.position  = rep(370, 4),
    label       = c("n.s.", "n.s.", "*", "n.s."),
    tip.length  = 0.005,
    inherit.aes = FALSE,
    size        = 1.5,
    label.size  = 6
  ) +
  guides(fill = "none")

# 8.7 Boxplot — time from symptom onset to diagnosis by specialist (Supplementary Figure)
ggplot(intervals_df,
  aes(x = specialist, y = interval_dx_symptoms, fill = final_diagnosis)) +
  geom_boxplot(width = 0.95, outlier.shape = NA) +
  labs(
    x    = "",
    y    = str_wrap("Time from symptom onset to diagnosis (weeks)", width = 35),
    fill = "Final conclusion"
  ) +
  scale_y_continuous(limits = c(0, 595), breaks = seq(0, 595, 72)) +
  scale_fill_manual(values = c("Myopathy" = "#FD6467", "Non-myopathy" = "#7294D4")) +
  theme_minimal() +
  theme(
    text               = element_text(size = 14),
    axis.line.y.left   = element_line(linewidth = 1.5),
    axis.line.x.bottom = element_line(linewidth = 1.5),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black"),
    legend.position    = "bottom"
  ) +
  geom_bracket(
    xmin        = c(0.7875, 1.7875, 2.7875, 3.7875),
    xmax        = c(1.2125, 2.2125, 3.2125, 4.2125),
    y.position  = rep(576, 4),
    label       = c("n.s.", "*", "**", "*"),
    tip.length  = 0.0003,
    inherit.aes = FALSE,
    size        = 1.5,
    label.size  = 6
  ) +
  guides(fill = "none")

