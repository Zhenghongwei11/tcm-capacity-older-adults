#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

analysis_tsv <- argval("--analysis", "local CHARLS analytic table")
out_tsv <- argval("--output", "results/tcm/tcm_supply_longitudinal_identification_models.tsv")

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

d <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  mutate(
    outcome_primary = as.numeric(primary_condition_tcm_any),
    z_tcm_beds_per_10000 = as.numeric(scale(value_per_10000_population_tcm_hospital_beds)),
    time_index = match(year, sort(unique(year))) - 1,
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile)
  ) %>%
  group_by(person_id) %>%
  mutate(
    person_n_waves = n_distinct(year),
    person_outcome_changed = n_distinct(outcome_primary) > 1,
    mean_z_tcm_beds = mean(z_tcm_beds_per_10000, na.rm = TRUE),
    mean_married = mean(married_or_partnered, na.rm = TRUE),
    mean_rural = mean(rural_community, na.rm = TRUE),
    mean_agricultural_hukou = mean(agricultural_hukou, na.rm = TRUE),
    mean_public_insurance = mean(public_insurance, na.rm = TRUE),
    mean_chronic_count = mean(chronic_count, na.rm = TRUE),
    mean_adl = mean(any_adl_limitation, na.rm = TRUE),
    mean_poor_health = mean(poor_self_rated_health, na.rm = TRUE)
  ) %>%
  ungroup()

full_covariates <- paste(
  "age_model + I(age_model^2) + female + education_group + married_or_partnered +",
  "rural_community + agricultural_hukou + public_insurance + chronic_count +",
  "any_adl_limitation + poor_self_rated_health + household_income_quartile"
)

time_varying_covariates <- paste(
  "married_or_partnered + rural_community + agricultural_hukou + public_insurance +",
  "chronic_count + any_adl_limitation + poor_self_rated_health + household_income_quartile"
)

mundlak_means <- paste(
  "mean_z_tcm_beds + mean_married + mean_rural + mean_agricultural_hukou +",
  "mean_public_insurance + mean_chronic_count + mean_adl + mean_poor_health"
)

fit_specs <- list(
  list(
    id = "I0_primary_reproduction",
    label = "province and wave fixed effects",
    formula = as.formula(paste0("outcome_primary ~ z_tcm_beds_per_10000 + ", full_covariates, " | province_supply_key + year")),
    data = d,
    identifying_variation = "within-province change across four survey waves"
  ),
  list(
    id = "I1_individual_fixed_effects",
    label = "individual and wave fixed effects",
    formula = as.formula(paste0("outcome_primary ~ z_tcm_beds_per_10000 + ", time_varying_covariates, " | person_id + year")),
    data = d %>% filter(person_n_waves >= 2),
    identifying_variation = "within-person change among respondents observed in at least two eligible waves"
  ),
  list(
    id = "I2_mundlak_correlated_effects",
    label = "province and wave fixed effects with person means",
    formula = as.formula(paste0("outcome_primary ~ z_tcm_beds_per_10000 + ", full_covariates, " + ", mundlak_means, " | province_supply_key + year")),
    data = d %>% filter(person_n_waves >= 2),
    identifying_variation = "within-person deviations with observed person-level means entered explicitly"
  ),
  list(
    id = "I3_province_linear_trends",
    label = "province and wave fixed effects plus province-specific linear trends",
    formula = as.formula(paste0("outcome_primary ~ z_tcm_beds_per_10000 + ", full_covariates, " + i(province_supply_key, time_index) | province_supply_key + year")),
    data = d,
    identifying_variation = "deviations from province-specific linear change across four survey waves"
  )
)

fit_one <- function(spec) {
  needed <- all.vars(spec$formula)
  used_data <- spec$data %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- feols(spec$formula, data = used_data, warn = FALSE, notes = FALSE)
  raw_cluster_vcov <- suppressWarnings(vcov(fit, vcov = ~province_supply_key, vcov_fix = FALSE))
  minimum_vcov_eigenvalue <- min(eigen(raw_cluster_vcov, symmetric = TRUE, only.values = TRUE)$values)
  vcov_psd_adjusted <- minimum_vcov_eigenvalue < -1e-10
  reported_cluster_vcov <- if (vcov_psd_adjusted) {
    suppressWarnings(vcov(fit, vcov = ~province_supply_key, vcov_fix = TRUE))
  } else {
    raw_cluster_vcov
  }
  term <- "z_tcm_beds_per_10000"
  raw_effect <- unname(coef(fit)[[term]])
  se <- unname(sqrt(diag(reported_cluster_vcov))[[term]])
  df <- n_distinct(used_data$province_supply_key) - 1
  ci <- raw_effect + c(-1, 1) * qt(0.975, df = df) * se
  tibble(
    model = spec$id,
    model_label = spec$label,
    identifying_variation = spec$identifying_variation,
    effect_type = "percentage-point difference per 1 SD higher province-year TCM hospital bed density",
    effect = 100 * raw_effect,
    ci_lower = 100 * ci[[1]],
    ci_upper = 100 * ci[[2]],
    pvalue = 2 * pt(abs(raw_effect / se), df = df, lower.tail = FALSE),
    n_persons = n_distinct(used_data$person_id),
    n_person_waves = nrow(used_data),
    n_repeated_persons = n_distinct(used_data$person_id[used_data$person_n_waves >= 2]),
    n_outcome_changers = n_distinct(used_data$person_id[used_data$person_outcome_changed]),
    n_clusters = n_distinct(used_data$province_supply_key),
    raw_cluster_vcov_minimum_eigenvalue = minimum_vcov_eigenvalue,
    cluster_vcov_psd_adjusted = vcov_psd_adjusted,
    inference = "province-clustered standard errors with t reference using G-1 degrees of freedom",
    source_table = analysis_tsv
  )
}

results <- bind_rows(lapply(fit_specs, fit_one)) %>%
  mutate(across(c(effect, ci_lower, ci_upper, pvalue), ~ round(.x, 5)))

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results)
