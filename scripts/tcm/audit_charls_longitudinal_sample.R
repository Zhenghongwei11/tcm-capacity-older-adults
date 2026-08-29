#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

analysis_tsv <- argval("--analysis", "results/tcm/person_wave_tcm_core_density_analysis.tsv")
weight_out <- argval("--weight-output", "results/tcm/charls_weight_audit.tsv")
flow_out <- argval("--flow-output", "results/tcm/charls_participant_flow.tsv")
pattern_out <- argval("--pattern-output", "results/tcm/charls_longitudinal_observation_patterns.tsv")
retention_out <- argval("--retention-output", "results/tcm/charls_retention_by_wave.tsv")
comparison_out <- argval("--comparison-output", "results/tcm/charls_attrition_comparison.tsv")
decision_out <- argval("--decision-output", "results/tcm/charls_attrition_weight_decision.tsv")

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(weight_out), recursive = TRUE, showWarnings = FALSE)

d <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  mutate(person_id = as.character(person_id))

eligible <- d %>% filter(main_model_age60)

weight_audit <- eligible %>%
  group_by(year) %>%
  summarize(
    eligible_person_waves = n(),
    positive_unadjusted_weight = sum(is.finite(respondent_weight) & respondent_weight > 0),
    positive_adjusted_weight = sum(is.finite(respondent_weight_adjusted) & respondent_weight_adjusted > 0),
    adjusted_weight_coverage_percent = 100 * positive_adjusted_weight / eligible_person_waves,
    adjusted_weight_min = min(respondent_weight_adjusted[respondent_weight_adjusted > 0], na.rm = TRUE),
    adjusted_weight_p01 = quantile(respondent_weight_adjusted[respondent_weight_adjusted > 0], 0.01, na.rm = TRUE),
    adjusted_weight_median = median(respondent_weight_adjusted[respondent_weight_adjusted > 0], na.rm = TRUE),
    adjusted_weight_p99 = quantile(respondent_weight_adjusted[respondent_weight_adjusted > 0], 0.99, na.rm = TRUE),
    adjusted_weight_max = max(respondent_weight_adjusted[respondent_weight_adjusted > 0], na.rm = TRUE),
    adjusted_weight_effective_n = {
      w <- respondent_weight_adjusted[is.finite(respondent_weight_adjusted) & respondent_weight_adjusted > 0]
      sum(w)^2 / sum(w^2)
    },
    .groups = "drop"
  ) %>%
  mutate(
    weight_definition = "wave-specific individual cross-sectional weight with household and individual non-response adjustment",
    pooled_use_rule = "normalize to mean 1 within wave; use as sensitivity rather than a design-based prevalence estimator",
    compatibility_note = "Comparable adjusted cross-sectional weight concept is available in all four waves; no common four-wave longitudinal weight exists."
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

full_covariates <- c(
  "age_model", "female", "education_group", "married_or_partnered",
  "rural_community", "agricultural_hukou", "public_insurance",
  "chronic_count", "any_adl_limitation", "poor_self_rated_health",
  "household_income_quartile"
)

flow <- bind_rows(
  tibble(stage = "linked province-year TCM panel", n_person_waves = sum(d$main_core_density_linked_panel), n_persons = n_distinct(d$person_id[d$main_core_density_linked_panel])),
  tibble(stage = "age 60 years or older using harmonized age", n_person_waves = sum(d$main_core_density_linked_panel & !is.na(d$age_model) & d$age_model >= 60, na.rm = TRUE), n_persons = n_distinct(d$person_id[d$main_core_density_linked_panel & !is.na(d$age_model) & d$age_model >= 60])),
  tibble(stage = "non-missing primary outcome", n_person_waves = nrow(eligible), n_persons = n_distinct(eligible$person_id)),
  tibble(
    stage = "complete primary-model covariates",
    n_person_waves = eligible %>% filter(if_all(all_of(full_covariates), ~ !is.na(.x))) %>% nrow(),
    n_persons = eligible %>% filter(if_all(all_of(full_covariates), ~ !is.na(.x))) %>% summarize(n = n_distinct(person_id)) %>% pull(n)
  ),
  tibble(
    stage = "positive adjusted respondent weight and complete primary-model covariates",
    n_person_waves = eligible %>% filter(if_all(all_of(full_covariates), ~ !is.na(.x)), is.finite(respondent_weight_adjusted), respondent_weight_adjusted > 0) %>% nrow(),
    n_persons = eligible %>% filter(if_all(all_of(full_covariates), ~ !is.na(.x)), is.finite(respondent_weight_adjusted), respondent_weight_adjusted > 0) %>% summarize(n = n_distinct(person_id)) %>% pull(n)
  )
)

person_patterns <- eligible %>%
  distinct(person_id, year, primary_condition_tcm_any) %>%
  arrange(person_id, year) %>%
  group_by(person_id) %>%
  summarize(
    observation_pattern = paste(year, collapse = "-"),
    n_observed_waves = n_distinct(year),
    first_observed_year = min(year),
    last_observed_year = max(year),
    outcome_changed = n_distinct(primary_condition_tcm_any) > 1,
    .groups = "drop"
  )

patterns <- person_patterns %>%
  count(observation_pattern, n_observed_waves, first_observed_year, last_observed_year, name = "n_persons") %>%
  mutate(percent_persons = 100 * n_persons / sum(n_persons)) %>%
  arrange(desc(n_persons)) %>%
  mutate(percent_persons = round(percent_persons, 3))

wave_pairs <- tibble(index_year = c(2011L, 2013L, 2015L), next_year = c(2013L, 2015L, 2018L))
retention_rows <- lapply(seq_len(nrow(wave_pairs)), function(i) {
  y0 <- wave_pairs$index_year[[i]]
  y1 <- wave_pairs$next_year[[i]]
  index <- eligible %>% filter(year == y0)
  next_ids <- eligible %>% filter(year == y1) %>% distinct(person_id)
  index %>%
    mutate(observed_next_wave = person_id %in% next_ids$person_id) %>%
    summarize(
      index_year = y0,
      next_year = y1,
      n_index = n(),
      n_observed_next_wave = sum(observed_next_wave),
      retention_percent = 100 * mean(observed_next_wave)
    )
})
retention <- bind_rows(retention_rows) %>% mutate(retention_percent = round(retention_percent, 3))

comparison_vars <- c(
  age_model = "continuous", female = "binary", rural_community = "binary",
  chronic_count = "continuous", any_adl_limitation = "binary",
  poor_self_rated_health = "binary", primary_condition_tcm_any = "binary"
)

comparison_rows <- list()
row_id <- 1L
for (i in seq_len(nrow(wave_pairs))) {
  y0 <- wave_pairs$index_year[[i]]
  y1 <- wave_pairs$next_year[[i]]
  next_ids <- eligible %>% filter(year == y1) %>% distinct(person_id)
  index <- eligible %>% filter(year == y0) %>% mutate(observed_next_wave = person_id %in% next_ids$person_id)
  for (v in names(comparison_vars)) {
    x1 <- index[[v]][index$observed_next_wave]
    x0 <- index[[v]][!index$observed_next_wave]
    pooled_sd <- sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
    comparison_rows[[row_id]] <- tibble(
      index_year = y0,
      next_year = y1,
      variable = v,
      variable_type = unname(comparison_vars[[v]]),
      observed_next_n = sum(!is.na(x1)),
      not_observed_next_n = sum(!is.na(x0)),
      observed_next_mean = mean(x1, na.rm = TRUE),
      not_observed_next_mean = mean(x0, na.rm = TRUE),
      standardized_mean_difference = ifelse(is.finite(pooled_sd) && pooled_sd > 0, (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / pooled_sd, NA_real_)
    )
    row_id <- row_id + 1L
  }
}
comparisons <- bind_rows(comparison_rows) %>% mutate(across(where(is.numeric), ~ round(.x, 4)))

ipow_decision <- tibble(
  method = "inverse probability of subsequent observation weighting",
  status = "not_implemented",
  reason = "The age-eligible panel combines ageing into eligibility, death, non-response, and replenishment; the available analysis table does not identify a common censoring mechanism across all four waves.",
  retained_analysis = "wave-specific adjusted respondent-weight sensitivity",
  implication = "Observation patterns and retention comparisons are reported directly; no attrition weight is presented as if it recovered a single four-wave target population."
)

write_tsv(weight_audit, weight_out)
write_tsv(flow, flow_out)
write_tsv(patterns, pattern_out)
write_tsv(retention, retention_out)
write_tsv(comparisons, comparison_out)
write_tsv(ipow_decision, decision_out)

cat(sprintf("Wrote %s\n", weight_out))
cat(sprintf("Wrote %s\n", flow_out))
cat(sprintf("Wrote %s\n", pattern_out))
cat(sprintf("Wrote %s\n", retention_out))
cat(sprintf("Wrote %s\n", comparison_out))
cat(sprintf("Wrote %s\n", decision_out))
