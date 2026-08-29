#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

analysis_tsv <- argval(
  "--analysis",
  "results/tcm/person_wave_tcm_core_density_analysis.tsv"
)
out_tsv <- argval(
  "--output",
  "results/tcm/tcm_supply_robustness.tsv"
)

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

base <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile)
  )

supply_panel <- base %>%
  distinct(
    province_supply_key, year,
    value_per_10000_population_tcm_practicing_assistant_physicians,
    value_per_10000_population_tcm_hospital_beds
  ) %>%
  arrange(province_supply_key, year) %>%
  group_by(province_supply_key) %>%
  mutate(
    lag_tcm_physicians_per_10000 = lag(value_per_10000_population_tcm_practicing_assistant_physicians),
    lag_tcm_beds_per_10000 = lag(value_per_10000_population_tcm_hospital_beds)
  ) %>%
  ungroup() %>%
  mutate(
    z_lag_tcm_physicians_per_10000 = as.numeric(scale(lag_tcm_physicians_per_10000)),
    z_lag_tcm_beds_per_10000 = as.numeric(scale(lag_tcm_beds_per_10000))
  ) %>%
  select(province_supply_key, year, starts_with("z_lag"))

analysis <- base %>%
  left_join(supply_panel, by = c("province_supply_key", "year")) %>%
  mutate(
    z_tcm_physicians_per_10000 = z_py_tcm_physicians_per_10000,
    z_tcm_beds_per_10000 = z_py_tcm_beds_per_10000,
    primary = as.numeric(primary_condition_tcm_any),
    strict_visit = as.numeric(strict_tcm_hospital_visit),
    broader = as.numeric(broader_tcm_use_2011_2015)
  )

cluster_vcov <- function(model, cluster) {
  ok <- !is.na(cluster)
  cluster <- as.factor(cluster[ok])
  x <- model.matrix(model)[ok, , drop = FALSE]
  u <- residuals(model)[ok]
  xtx_inv <- solve(crossprod(x))
  meat <- matrix(0, ncol(x), ncol(x))
  for (g in levels(cluster)) {
    idx <- cluster == g
    xg <- x[idx, , drop = FALSE]
    ug <- u[idx]
    sg <- crossprod(xg, ug)
    meat <- meat + tcrossprod(sg)
  }
  n <- nrow(x)
  k <- ncol(x)
  g <- nlevels(cluster)
  correction <- (g / (g - 1)) * ((n - 1) / (n - k))
  correction * xtx_inv %*% meat %*% xtx_inv
}

fit_one <- function(data, outcome_var, supply_var, model_id, check_type, covariates, years, age_filter, multiplier = 100) {
  rhs <- c(supply_var, covariates, "province_fe", "wave_fe")
  f <- as.formula(paste(outcome_var, "~", paste(rhs, collapse = " + ")))
  needed <- all.vars(f)
  d <- data %>%
    filter(year %in% years) %>%
    filter(.data[[age_filter]]) %>%
    filter(if_all(all_of(needed), ~ !is.na(.x)))

  fit <- lm(f, data = d)
  vc <- cluster_vcov(fit, d$province_supply_key)
  coefs <- coef(fit)
  se <- sqrt(diag(vc))
  effect <- unname(coefs[[supply_var]])
  se_term <- unname(se[[supply_var]])
  df <- length(unique(d$province_supply_key)) - 1
  ci <- effect + c(-1, 1) * qt(0.975, df = df) * se_term
  p_value <- 2 * pt(abs(effect / se_term), df = df, lower.tail = FALSE)

  tibble(
    check_type = check_type,
    outcome = case_when(
      outcome_var == "primary" ~ "Any disease/condition-specific Chinese medicine treatment use",
      outcome_var == "strict_visit" ~ "TCM hospital visit in the past month",
      outcome_var == "broader" ~ "Broader TCM-related use in 2011-2015",
      TRUE ~ outcome_var
    ),
    supply_indicator = case_when(
      supply_var == "z_tcm_physicians_per_10000" ~ "TCM practicing/assistant physicians per 10,000 population",
      supply_var == "z_tcm_beds_per_10000" ~ "TCM hospital beds per 10,000 population",
      supply_var == "z_lag_tcm_physicians_per_10000" ~ "Lagged TCM practicing/assistant physicians per 10,000 population",
      supply_var == "z_lag_tcm_beds_per_10000" ~ "Lagged TCM hospital beds per 10,000 population",
      TRUE ~ supply_var
    ),
    model = model_id,
    effect_type = "percentage-point difference per 1 SD higher province-year supply",
    effect = multiplier * effect,
    ci_lower = multiplier * ci[[1]],
    ci_upper = multiplier * ci[[2]],
    pvalue = p_value,
    n_persons = n_distinct(d$person_id),
    n_person_waves = nrow(d),
    cluster_level = "province",
    n_clusters = n_distinct(d$province_supply_key),
    analysis_years = paste(sort(unique(d$year)), collapse = ";"),
    covariate_set = paste(covariates, collapse = "; "),
    fixed_effects = "province and survey year",
    source_table = analysis_tsv
  )
}

full_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile"
)
no_income_covariates <- setdiff(full_covariates, "household_income_quartile")

analysis <- analysis %>%
  mutate(
    age60_filter = main_model_age60,
    age45_filter = main_core_density_linked_panel & age_model >= 45
  )

results <- bind_rows(
  fit_one(analysis, "primary", "z_tcm_beds_per_10000", "R1_primary_age60_full", "main_reproduction", full_covariates, c(2011L, 2013L, 2015L, 2018L), "age60_filter"),
  fit_one(analysis, "primary", "z_tcm_physicians_per_10000", "R1_primary_age60_full", "alternative_supply_indicator", full_covariates, c(2011L, 2013L, 2015L, 2018L), "age60_filter"),
  fit_one(analysis, "primary", "z_tcm_beds_per_10000", "R2_primary_age45_full", "age_threshold", full_covariates, c(2011L, 2013L, 2015L, 2018L), "age45_filter"),
  fit_one(analysis, "primary", "z_tcm_beds_per_10000", "R3_primary_age60_no_income", "covariate_set", no_income_covariates, c(2011L, 2013L, 2015L, 2018L), "age60_filter"),
  fit_one(analysis, "primary", "z_lag_tcm_beds_per_10000", "R4_primary_lagged_supply", "lag_structure", full_covariates, c(2013L, 2015L, 2018L), "age60_filter"),
  fit_one(analysis, "primary", "z_lag_tcm_physicians_per_10000", "R4_primary_lagged_supply", "lag_structure", full_covariates, c(2013L, 2015L, 2018L), "age60_filter"),
  fit_one(analysis, "strict_visit", "z_tcm_beds_per_10000", "R5_strict_visit_age60_full", "alternative_outcome", full_covariates, c(2011L, 2013L, 2015L, 2018L), "age60_filter"),
  fit_one(analysis, "broader", "z_tcm_beds_per_10000", "R6_broader_use_age60_full", "alternative_outcome", full_covariates, c(2011L, 2013L, 2015L), "age60_filter")
) %>%
  mutate(across(c(effect, ci_lower, ci_upper, pvalue), ~ round(.x, 5)))

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results)
