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
component_out <- argval(
  "--component-output",
  "results/tcm/tcm_supply_primary_component_models.tsv"
)
cluster_out <- argval(
  "--cluster-output",
  "results/tcm/tcm_supply_cluster_sensitivity.tsv"
)
bootstrap_reps <- as.integer(argval("--bootstrap-reps", "499"))

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(component_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cluster_out), recursive = TRUE, showWarnings = FALSE)

full_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile"
)

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    z_tcm_beds_per_10000 = z_py_tcm_beds_per_10000,
    primary = as.numeric(primary_condition_tcm_any)
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

fit_lpm <- function(data, outcome_var, covariates = full_covariates) {
  rhs <- c("z_tcm_beds_per_10000", covariates, "province_fe", "wave_fe")
  f <- as.formula(paste(outcome_var, "~", paste(rhs, collapse = " + ")))
  needed <- all.vars(f)
  d <- data %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- lm(f, data = d)
  vc <- cluster_vcov(fit, d$province_supply_key)
  raw_effect <- unname(coef(fit)[["z_tcm_beds_per_10000"]])
  raw_se <- unname(sqrt(diag(vc))[["z_tcm_beds_per_10000"]])
  df <- n_distinct(d$province_supply_key) - 1
  ci <- raw_effect + c(-1, 1) * qt(0.975, df = df) * raw_se
  tibble(
    effect = 100 * raw_effect,
    ci_lower = 100 * ci[[1]],
    ci_upper = 100 * ci[[2]],
    pvalue = 2 * pt(abs(raw_effect / raw_se), df = df, lower.tail = FALSE),
    n_persons = n_distinct(d$person_id),
    n_person_waves = nrow(d),
    n_events = sum(d[[outcome_var]] == 1, na.rm = TRUE),
    event_percent = 100 * mean(d[[outcome_var]] == 1, na.rm = TRUE),
    n_clusters = n_distinct(d$province_supply_key)
  )
}

component_vars <- grep("^condition_tcm_", names(analysis), value = TRUE)

component_results <- bind_rows(lapply(component_vars, function(var) {
  d <- analysis %>% mutate(component_outcome = as.numeric(.data[[var]]))
  event_count <- sum(d$component_outcome == 1, na.rm = TRUE)
  if (event_count < 50) {
    return(tibble(
      outcome = var,
      status = "not_modeled_sparse_events",
      effect_type = "percentage-point difference per 1 SD higher TCM hospital bed density",
      effect = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      pvalue = NA_real_,
      n_persons = n_distinct(d$person_id),
      n_person_waves = nrow(d),
      n_events = event_count,
      event_percent = 100 * mean(d$component_outcome == 1, na.rm = TRUE),
      n_clusters = n_distinct(d$province_supply_key),
      fixed_effects = "province and survey year",
      covariate_set = paste(full_covariates, collapse = "; "),
      source_table = analysis_tsv
    ))
  }
  fit_lpm(d, "component_outcome") %>%
    mutate(
      outcome = var,
      status = "modeled",
      effect_type = "percentage-point difference per 1 SD higher TCM hospital bed density",
      fixed_effects = "province and survey year",
      covariate_set = paste(full_covariates, collapse = "; "),
      source_table = analysis_tsv,
      .before = effect
    )
})) %>%
  select(
    outcome, status, effect_type, effect, ci_lower, ci_upper, pvalue,
    n_persons, n_person_waves, n_events, event_percent, n_clusters,
    fixed_effects, covariate_set, source_table
  ) %>%
  mutate(across(c(effect, ci_lower, ci_upper, pvalue, event_percent), ~ round(.x, 5)))

write_tsv(component_results, component_out)

main_estimate <- fit_lpm(analysis %>% mutate(outcome_primary = primary), "outcome_primary") %>%
  mutate(check_type = "main_cluster_robust", detail = "all provinces", .before = effect)

province_keys <- sort(unique(analysis$province_supply_key))
loo <- bind_rows(lapply(province_keys, function(prov) {
  fit_lpm(
    analysis %>%
      filter(province_supply_key != prov) %>%
      mutate(outcome_primary = primary),
    "outcome_primary"
  ) %>%
    mutate(check_type = "leave_one_province_out", detail = prov, .before = effect)
}))

set.seed(20260718)
boot_data <- analysis %>% mutate(outcome_primary = primary)
base_formula <- as.formula(paste(
  "outcome_primary ~",
  paste(c("z_tcm_beds_per_10000", full_covariates, "province_fe", "wave_fe"), collapse = " + ")
))
needed <- all.vars(base_formula)
boot_data <- boot_data %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
clusters <- sort(unique(boot_data$province_supply_key))

boot_effects <- replicate(bootstrap_reps, {
  sampled <- sample(clusters, length(clusters), replace = TRUE)
  cluster_weights <- tibble(province_supply_key = sampled) %>%
    count(province_supply_key, name = "cluster_weight")
  d <- boot_data %>%
    left_join(cluster_weights, by = "province_supply_key") %>%
    mutate(cluster_weight = if_else(is.na(cluster_weight), 0L, cluster_weight))
  fit <- lm(base_formula, data = d, weights = cluster_weight)
  100 * unname(coef(fit)[["z_tcm_beds_per_10000"]])
})

cluster_results <- bind_rows(
  main_estimate,
  tibble(
    check_type = "leave_one_province_out_summary",
    detail = "effect range across omitted provinces",
    effect = median(loo$effect, na.rm = TRUE),
    ci_lower = min(loo$effect, na.rm = TRUE),
    ci_upper = max(loo$effect, na.rm = TRUE),
    pvalue = NA_real_,
    n_persons = main_estimate$n_persons,
    n_person_waves = main_estimate$n_person_waves,
    n_events = main_estimate$n_events,
    event_percent = main_estimate$event_percent,
    n_clusters = main_estimate$n_clusters
  ),
  tibble(
    check_type = "province_cluster_bootstrap",
    detail = paste0(bootstrap_reps, " weighted cluster-bootstrap replicates"),
    effect = median(boot_effects, na.rm = TRUE),
    ci_lower = unname(quantile(boot_effects, 0.025, na.rm = TRUE)),
    ci_upper = unname(quantile(boot_effects, 0.975, na.rm = TRUE)),
    pvalue = NA_real_,
    n_persons = n_distinct(boot_data$person_id),
    n_person_waves = nrow(boot_data),
    n_events = sum(boot_data$outcome_primary == 1, na.rm = TRUE),
    event_percent = 100 * mean(boot_data$outcome_primary == 1, na.rm = TRUE),
    n_clusters = n_distinct(boot_data$province_supply_key)
  ),
  loo
) %>%
  mutate(
    outcome = "Any disease/condition-specific Chinese medicine treatment use",
    supply_indicator = "TCM hospital beds per 10,000 population",
    effect_type = "percentage-point difference per 1 SD higher province-year supply",
    fixed_effects = "province and survey year",
    covariate_set = paste(full_covariates, collapse = "; "),
    source_table = analysis_tsv
  ) %>%
  select(
    check_type, detail, outcome, supply_indicator, effect_type, effect, ci_lower,
    ci_upper, pvalue, n_persons, n_person_waves, n_events, event_percent,
    n_clusters, fixed_effects, covariate_set, source_table
  ) %>%
  mutate(across(c(effect, ci_lower, ci_upper, pvalue, event_percent), ~ round(.x, 5)))

write_tsv(cluster_results, cluster_out)

cat(sprintf("Wrote %s\n", component_out))
print(component_results)
cat(sprintf("Wrote %s\n", cluster_out))
print(cluster_results %>% filter(check_type != "leave_one_province_out"))
