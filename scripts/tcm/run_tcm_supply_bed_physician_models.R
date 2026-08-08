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
context_tsv <- argval(
  "--context",
  "results/tcm/province_year_contextual_covariates.tsv"
)
out_tsv <- argval(
  "--output",
  "results/tcm/tcm_supply_bed_physician_models.tsv"
)
diagnostics_tsv <- argval(
  "--diagnostics-output",
  "results/tcm/tcm_supply_bed_physician_diagnostics.tsv"
)

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
if (!file.exists(context_tsv)) stop(sprintf("Missing contextual table: %s", context_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

cluster_vcov <- function(model, cluster) {
  cluster <- droplevels(as.factor(cluster))
  x <- model.matrix(model)
  u <- residuals(model)
  xtx_inv <- qr.solve(crossprod(x))
  meat <- matrix(0, ncol(x), ncol(x))
  for (g in levels(cluster)) {
    idx <- cluster == g
    score_g <- crossprod(x[idx, , drop = FALSE], u[idx])
    meat <- meat + tcrossprod(score_g)
  }
  n <- nrow(x)
  k <- ncol(x)
  groups <- nlevels(cluster)
  correction <- (groups / (groups - 1)) * ((n - 1) / (n - k))
  out <- correction * xtx_inv %*% meat %*% xtx_inv
  dimnames(out) <- list(colnames(x), colnames(x))
  out
}

context <- read_tsv(context_tsv, show_col_types = FALSE) %>%
  transmute(
    province_supply_key = province,
    year = as.integer(year),
    comprehensive_beds = value_per_10000_population_comprehensive_hospital_beds,
    log_gdp_per_capita = log_gdp_per_capita_current_yuan,
    urban_population_percent
  )

supply_panel <- read_tsv(context_tsv, show_col_types = FALSE) %>%
  transmute(
    province_supply_key = province,
    year = as.integer(year),
    tcm_beds = value_per_10000_population_tcm_hospital_beds,
    tcm_physicians = value_per_10000_population_tcm_practicing_assistant_physicians
  ) %>%
  arrange(province_supply_key, year) %>%
  group_by(province_supply_key) %>%
  mutate(
    lag_tcm_beds = lag(tcm_beds),
    lag_tcm_physicians = lag(tcm_physicians)
  ) %>%
  ungroup()

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_core_density_linked_panel, age_60plus == 1) %>%
  left_join(
    supply_panel %>% select(province_supply_key, year, lag_tcm_beds, lag_tcm_physicians),
    by = c("province_supply_key", "year")
  ) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    outcome_primary = as.numeric(primary_condition_tcm_any),
    z_tcm_beds = z_py_tcm_beds_per_10000,
    z_tcm_physicians = z_py_tcm_physicians_per_10000,
    z_lag_tcm_beds = as.numeric(scale(lag_tcm_beds)),
    z_lag_tcm_physicians = as.numeric(scale(lag_tcm_physicians)),
    z_comprehensive_beds = z_py_comprehensive_hospital_beds,
    z_log_gdp = z_py_log_gdp,
    z_urbanization = z_py_urbanization
  )

individual_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile"
)

model_specs <- list(
  J1 = c("z_tcm_beds", "z_tcm_physicians"),
  J2 = c(
    "z_tcm_beds", "z_tcm_physicians", "z_comprehensive_beds",
    "z_log_gdp", "z_urbanization"
  ),
  J3 = c("z_lag_tcm_beds", "z_lag_tcm_physicians"),
  J4 = c(
    "z_lag_tcm_beds", "z_lag_tcm_physicians", "z_comprehensive_beds",
    "z_log_gdp", "z_urbanization"
  )
)

model_labels <- c(
  J1 = "Concurrent TCM beds and physicians",
  J2 = "Concurrent TCM resources with contextual adjustment",
  J3 = "Lagged TCM beds and physicians",
  J4 = "Lagged TCM resources with contextual adjustment"
)

term_labels <- c(
  z_tcm_beds = "Concurrent TCM hospital beds per 10,000 population",
  z_tcm_physicians = "Concurrent TCM physicians per 10,000 population",
  z_lag_tcm_beds = "Lagged TCM hospital beds per 10,000 population",
  z_lag_tcm_physicians = "Lagged TCM physicians per 10,000 population",
  z_comprehensive_beds = "Comprehensive-hospital beds per 10,000 population",
  z_log_gdp = "Log provincial GDP per capita",
  z_urbanization = "Urban population share"
)

fit_spec <- function(model_id, terms) {
  rhs <- c(terms, individual_covariates, "province_fe", "wave_fe")
  formula <- as.formula(paste("outcome_primary ~", paste(rhs, collapse = " + ")))
  needed <- all.vars(formula)
  d <- analysis %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- lm(formula, data = d)
  vcov <- cluster_vcov(fit, d$province_supply_key)
  df <- n_distinct(d$province_supply_key) - 1

  bind_rows(lapply(terms, function(term) {
    effect <- unname(coef(fit)[[term]])
    se <- sqrt(vcov[term, term])
    ci <- effect + c(-1, 1) * qt(0.975, df = df) * se
    pvalue <- 2 * pt(abs(effect / se), df = df, lower.tail = FALSE)
    tibble(
      model = model_id,
      model_label = unname(model_labels[[model_id]]),
      term = term,
      term_label = unname(term_labels[[term]]),
      effect_type = "percentage-point difference per 1 SD higher province-year measure",
      effect = 100 * effect,
      ci_lower = 100 * ci[[1]],
      ci_upper = 100 * ci[[2]],
      cluster_pvalue = pvalue,
      n_persons = n_distinct(d$person_id),
      n_person_waves = nrow(d),
      n_clusters = n_distinct(d$province_supply_key),
      analysis_years = paste(sort(unique(d$year)), collapse = ";"),
      fixed_effects = "province and survey year"
    )
  }))
}

results <- bind_rows(Map(fit_spec, names(model_specs), model_specs)) %>%
  mutate(across(c(effect, ci_lower, ci_upper, cluster_pvalue), ~ round(.x, 5)))

panel <- analysis %>%
  distinct(
    province_supply_key, year, z_tcm_beds, z_tcm_physicians,
    z_lag_tcm_beds, z_lag_tcm_physicians, z_comprehensive_beds,
    z_log_gdp, z_urbanization
  )

diagnose_spec <- function(model_id, terms) {
  d <- panel %>% filter(if_all(all_of(terms), ~ !is.na(.x)))
  residualized <- lapply(terms, function(term) {
    residuals(lm(
      as.formula(paste(term, "~ factor(province_supply_key) + factor(year)")),
      data = d
    ))
  })
  x <- scale(do.call(cbind, residualized))
  colnames(x) <- terms
  correlation <- cor(x)
  eigenvalues <- eigen(crossprod(x), symmetric = TRUE, only.values = TRUE)$values
  positive <- eigenvalues[eigenvalues > max(eigenvalues) * .Machine$double.eps]
  condition_index <- sqrt(max(positive) / min(positive))

  bind_rows(lapply(seq_along(terms), function(j) {
    other <- setdiff(seq_along(terms), j)
    vif <- if (length(other) == 0) 1 else
      1 / (1 - summary(lm(x[, j] ~ x[, other, drop = FALSE]))$r.squared)
    tibble(
      model = model_id,
      term = terms[[j]],
      term_label = unname(term_labels[[terms[[j]]]]),
      two_way_fe_vif = vif,
      maximum_absolute_pairwise_correlation = if (length(other) == 0) 0 else max(abs(correlation[j, other])),
      design_condition_index = condition_index,
      province_years = nrow(d),
      provinces = n_distinct(d$province_supply_key),
      waves = n_distinct(d$year),
      diagnostic_scale = "province-year exposure residuals after province and survey-year fixed effects"
    )
  }))
}

diagnostics <- bind_rows(Map(diagnose_spec, names(model_specs), model_specs)) %>%
  mutate(across(c(two_way_fe_vif, maximum_absolute_pairwise_correlation, design_condition_index), ~ round(.x, 5)))

write_tsv(results, out_tsv)
write_tsv(diagnostics, diagnostics_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
cat(sprintf("Wrote %s\n", diagnostics_tsv))
print(results %>% filter(grepl("tcm_beds|tcm_physicians", term)) %>% select(model, term, effect, ci_lower, ci_upper, cluster_pvalue, n_person_waves))
print(diagnostics %>% filter(grepl("tcm_beds|tcm_physicians", term)) %>% select(model, term, two_way_fe_vif, maximum_absolute_pairwise_correlation, design_condition_index))
