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

analysis_tsv <- argval("--analysis", "local CHARLS analytic table")
out_tsv <- argval("--output", "results/tcm/tcm_supply_weighted_attrition_models.tsv")

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

full_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile"
)

d <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    outcome_primary = as.numeric(primary_condition_tcm_any),
    z_tcm_beds_per_10000 = as.numeric(scale(value_per_10000_population_tcm_hospital_beds))
  ) %>%
  group_by(year) %>%
  mutate(
    adjusted_weight_normalized = respondent_weight_adjusted / mean(respondent_weight_adjusted[respondent_weight_adjusted > 0], na.rm = TRUE),
    weight_p01 = quantile(adjusted_weight_normalized[adjusted_weight_normalized > 0], 0.01, na.rm = TRUE),
    weight_p99 = quantile(adjusted_weight_normalized[adjusted_weight_normalized > 0], 0.99, na.rm = TRUE),
    adjusted_weight_trimmed = pmin(pmax(adjusted_weight_normalized, weight_p01), weight_p99)
  ) %>%
  ungroup()

cluster_vcov_wls <- function(model, data, cluster_var, weight_var = NULL) {
  x <- model.matrix(model)
  u <- residuals(model)
  cluster <- as.factor(data[[cluster_var]])
  w <- if (is.null(weight_var)) rep(1, nrow(data)) else data[[weight_var]]
  bread <- solve(crossprod(x, w * x))
  meat <- matrix(0, ncol(x), ncol(x))
  for (g in levels(cluster)) {
    idx <- cluster == g
    score <- crossprod(x[idx, , drop = FALSE], w[idx] * u[idx])
    meat <- meat + tcrossprod(score)
  }
  n <- nrow(x)
  k <- ncol(x)
  g <- nlevels(cluster)
  (g / (g - 1)) * ((n - 1) / (n - k)) * bread %*% meat %*% bread
}

fit_one <- function(data, model_id, weight_var = NULL, weight_definition) {
  f <- as.formula(paste("outcome_primary ~", paste(c("z_tcm_beds_per_10000", full_covariates, "province_fe", "wave_fe"), collapse = " + ")))
  needed <- unique(c(all.vars(f), weight_var))
  dd <- data %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  if (!is.null(weight_var)) dd <- dd %>% filter(.data[[weight_var]] > 0, is.finite(.data[[weight_var]]))
  fit <- if (is.null(weight_var)) lm(f, data = dd) else lm(f, data = dd, weights = dd[[weight_var]])
  vc <- cluster_vcov_wls(fit, dd, "province_supply_key", weight_var)
  raw_effect <- unname(coef(fit)[["z_tcm_beds_per_10000"]])
  se <- unname(sqrt(diag(vc))[["z_tcm_beds_per_10000"]])
  df <- n_distinct(dd$province_supply_key) - 1
  ci <- raw_effect + c(-1, 1) * qt(0.975, df = df) * se
  w <- if (is.null(weight_var)) rep(1, nrow(dd)) else dd[[weight_var]]
  tibble(
    model = model_id,
    weighting = weight_definition,
    effect_type = "percentage-point difference per 1 SD higher province-year TCM hospital bed density",
    effect = 100 * raw_effect,
    ci_lower = 100 * ci[[1]],
    ci_upper = 100 * ci[[2]],
    pvalue = 2 * pt(abs(raw_effect / se), df = df, lower.tail = FALSE),
    n_persons = n_distinct(dd$person_id),
    n_person_waves = nrow(dd),
    effective_person_waves = sum(w)^2 / sum(w^2),
    n_clusters = n_distinct(dd$province_supply_key),
    fixed_effects = "province and survey year",
    attrition_weight = "not used; observation mechanisms are not jointly identifiable in the age-eligible pooled panel",
    source_table = analysis_tsv
  )
}

results <- bind_rows(
  fit_one(d, "W0_unweighted_reproduction", NULL, "none"),
  fit_one(d, "W1_adjusted_cross_sectional_weight", "adjusted_weight_normalized", "wave-specific individual weight with household and individual non-response adjustment, normalized within wave"),
  fit_one(d, "W2_trimmed_adjusted_cross_sectional_weight", "adjusted_weight_trimmed", "W1 weight truncated at wave-specific 1st and 99th percentiles")
) %>% mutate(across(where(is.numeric), ~ round(.x, 5)))

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results)
