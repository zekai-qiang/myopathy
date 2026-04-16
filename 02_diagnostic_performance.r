# =============================================================================
# Myopathy Diagnostic Journey: Diagnostic Performance Analysis
# =============================================================================
# Publication: "The diagnostic journey of patients being investigated for
#              myopathy in a tertiary centre in England"
#              Journal of Neurology (2025) 272:35
#              https://doi.org/10.1007/s00415-024-12737-y
#
# Description: Computes sensitivity, specificity, PPV and NPV for each
#              diagnostic investigation (Table 3), and evaluates the
#              performance of test combinations using an OR operator
#              (Table 4 / Figure 1d analogue).
#
# Input:  data/myopathy_main.csv   — main patient-level dataset
#         data/diagnosis_results.csv — per-test binary correctness dataset
#                                      (columns: study_id, final_diagnosis,
#                                       ck, myositis_abs, emg, muscle_mri,
#                                       biopsy, genetics)
#
# Note:   Data are available upon request from the corresponding author.
#         See Data Dictionary in 01_myopathy_main_analysis.R for column
#         specifications.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries
# -----------------------------------------------------------------------------
library(tidyverse)
library(caret)
library(gtools)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
# Replace paths with your local file paths.
data_raw       <- read.csv("data/myopathy_main.csv",       stringsAsFactors = FALSE)
diagnosis_results <- read.csv("data/diagnosis_results.csv", stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 2. Diagnostic Performance — Sensitivity, Specificity, PPV, NPV (Table 3)
# -----------------------------------------------------------------------------
# Recode raw investigation values to: "Myopathy", "Non-myopathy", "Not done"
# Coding: 0 = not performed; 1 = non-myopathy result; 2 or 3 = myopathy result

investigation_performance_cols <- c("ck", "myositis_abs", "emg",
                                    "muscle_mri", "biopsy", "genetics")

results_coded <- data_raw %>%
  select(study_id, all_of(investigation_performance_cols), final_diagnosis) %>%
  mutate(across(
    all_of(investigation_performance_cols),
    ~ case_when(
      . == 0       ~ "Not done",
      . == 1       ~ "Non-myopathy",
      . %in% c(2, 3) ~ "Myopathy",
      TRUE         ~ NA_character_
    )
  ))

# Compute confusion matrix metrics for each investigation
compute_performance <- function(data, col) {
  df_filtered <- data %>%
    filter(.data[[col]] != "Not done") %>%
    mutate(
      pred = factor(.data[[col]],    levels = c("Non-myopathy", "Myopathy")),
      ref  = factor(final_diagnosis, levels = c("Non-myopathy", "Myopathy"))
    )

  confusionMatrix(
    data      = df_filtered$pred,
    reference = df_filtered$ref,
    positive  = "Myopathy"
  )
}

performance_results <- lapply(
  setNames(investigation_performance_cols, investigation_performance_cols),
  function(col) compute_performance(results_coded, col)
)

# Print sensitivity, specificity, PPV, NPV for each test
lapply(names(performance_results), function(test) {
  cm <- performance_results[[test]]
  cat(sprintf(
    "\n%s\n  n = %d | Sensitivity = %.1f%% | Specificity = %.1f%% | PPV = %.1f%% | NPV = %.1f%%\n",
    test,
    sum(cm$table),
    cm$byClass["Sensitivity"]  * 100,
    cm$byClass["Specificity"]  * 100,
    cm$byClass["Pos Pred Value"] * 100,
    cm$byClass["Neg Pred Value"] * 100
  ))
})

# -----------------------------------------------------------------------------
# 3. Test Combination Performance (Table 4)
# -----------------------------------------------------------------------------
# For each combination of tests, calculate the proportion of patients where
# at least one test result was correct relative to the final diagnosis (OR rule).
#
# The diagnosis_results dataset contains one row per patient; each test column
# holds 1 (correct result) or 0 (incorrect result), with NA if not performed.
# Only patients for whom all tests in a given combination were performed
# (no NA / 0-coded "not done") are included in the denominator.

# 3.1 Helper — compute OR-rule combination metrics
compute_combination_metrics <- function(combo, results_data) {
  subset_data <- results_data[, combo, drop = FALSE]

  # Retain only rows where all tests in the combination were performed
  rows_performed <- !apply(subset_data, 1, function(row) any(row == 0 | is.na(row)))
  subset_performed <- subset_data[rows_performed, , drop = FALSE]

  n_performed <- nrow(subset_performed)
  n_correct   <- sum(apply(subset_performed, 1, function(row) any(row == 1)))
  pct_correct <- if (n_performed > 0) (n_correct / n_performed) * 100 else NA_real_

  data.frame(
    combination  = paste(combo, collapse = " + "),
    n_performed  = n_performed,
    n_correct    = n_correct,
    pct_correct  = round(pct_correct, 1)
  )
}

# 3.2 Generate all combinations of 2 or more tests
combination_tests <- colnames(diagnosis_results)[2:8]

all_combinations <- lapply(2:length(combination_tests), function(k) {
  combn(combination_tests, k, simplify = FALSE)
}) %>% unlist(recursive = FALSE)

# 3.3 Compute metrics for all subjects
combination_results_all <- lapply(all_combinations, function(combo) {
  compute_combination_metrics(combo, diagnosis_results)
}) %>%
  bind_rows() %>%
  arrange(desc(n_performed))

combination_results_all

# 3.4 Compute metrics for myopathy subgroup
diagnosis_results_myopathy <- diagnosis_results %>%
  filter(final_diagnosis == "Myopathy")

combination_results_myopathy <- lapply(all_combinations, function(combo) {
  compute_combination_metrics(combo, diagnosis_results_myopathy)
}) %>%
  bind_rows() %>%
  arrange(desc(n_performed))

combination_results_myopathy

# 3.5 Compute metrics for non-myopathy subgroup
diagnosis_results_nonmyopathy <- diagnosis_results %>%
  filter(final_diagnosis == "Non-myopathy")

combination_results_nonmyopathy <- lapply(all_combinations, function(combo) {
  compute_combination_metrics(combo, diagnosis_results_nonmyopathy)
}) %>%
  bind_rows() %>%
  arrange(desc(n_performed))

combination_results_nonmyopathy

