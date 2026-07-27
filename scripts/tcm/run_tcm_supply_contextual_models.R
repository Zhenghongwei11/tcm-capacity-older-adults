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
  "local CHARLS analytic table"
)
context_tsv <- argval(
  "--context",
  "results/tcm/province_year_contextual_covariates.tsv"
)
out_tsv <- argval(
  "--output",
  "results/tcm/tcm_supply_contextual_models.tsv"
)
diagnostics_tsv <- argval(
  "--diagnostics-output",
  "results/tcm/tcm_supply_contextual_diagnostics.tsv"
)
bootstrap_reps <- as.integer(argval("--bootstrap-reps", "9999"))
bootstrap_seed <- as.integer(argval("--bootstrap-seed", "20260720"))

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
  g <- nlevels(cluster)
  correction <- (g / (g - 1)) * ((n - 1) / (n - k))
  out <- correction * xtx_inv %*% meat %*% xtx_inv
  dimnames(out) <- list(colnames(x), colnames(x))
  out
}

wild_cluster_pvalue <- function(model, term, cluster, reps, seed) {
  cluster <- droplevels(as.factor(cluster))
  cluster_index <- as.integer(cluster)
  g <- nlevels(cluster)
  x <- model.matrix(model)
  y <- model.response(model.frame(model))
  term_col <- match(term, colnames(x))
  if (is.na(term_col)) stop(sprintf("Term not found in model matrix: %s", term))

  x0 <- x[, -term_col, drop = FALSE]
  beta0 <- qr.solve(crossprod(x0), crossprod(x0, y))
  fitted0 <- as.numeric(x0 %*% beta0)
  residual0 <- y - fitted0
  xtx_inv <- qr.solve(crossprod(x))

  observed_vcov <- cluster_vcov(model, cluster)
  observed_t <- unname(coef(model)[[term]] / sqrt(observed_vcov[term, term]))
  webb_weights <- c(-sqrt(3 / 2), -1, -sqrt(1 / 2), sqrt(1 / 2), 1, sqrt(3 / 2))

  set.seed(seed)
  boot_t <- rep(NA_real_, reps)
  for (b in seq_len(reps)) {
    weights <- sample(webb_weights, g, replace = TRUE)
    y_star <- fitted0 + residual0 * weights[cluster_index]
    beta_star <- as.numeric(xtx_inv %*% crossprod(x, y_star))
    residual_star <- y_star - as.numeric(x %*% beta_star)

    meat <- matrix(0, ncol(x), ncol(x))
    for (cluster_id in seq_len(g)) {
      idx <- cluster_index == cluster_id
      score_g <- crossprod(x[idx, , drop = FALSE], residual_star[idx])
      meat <- meat + tcrossprod(score_g)
    }
    n <- nrow(x)
    k <- ncol(x)
    correction <- (g / (g - 1)) * ((n - 1) / (n - k))
    vcov_star <- correction * xtx_inv %*% meat %*% xtx_inv
    se_star <- sqrt(vcov_star[term_col, term_col])
    if (is.finite(se_star) && se_star > 0) boot_t[[b]] <- beta_star[[term_col]] / se_star
  }

  valid <- is.finite(boot_t)
  tibble(
    wild_cluster_pvalue = (1 + sum(abs(boot_t[valid]) >= abs(observed_t))) /
      (1 + sum(valid)),
    wild_cluster_reps_requested = reps,
    wild_cluster_reps_valid = sum(valid),
    wild_cluster_weight = "Webb six-point",
    wild_cluster_null_imposed = TRUE
  )
}

context <- read_tsv(context_tsv, show_col_types = FALSE) %>%
  transmute(
    province_supply_key = province,
    year = as.integer(year),
    comprehensive_beds = value_per_10000_population_comprehensive_hospital_beds,
    total_hospital_beds = value_per_10000_population_hospital_beds,
    log_gdp_per_capita = log_gdp_per_capita_current_yuan,
    urban_population_percent,
    population_age_65_plus_percent
  )

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    outcome_primary = as.numeric(primary_condition_tcm_any),
    z_tcm_beds = as.numeric(scale(value_per_10000_population_tcm_hospital_beds)),
    z_comprehensive_beds = as.numeric(scale(value_per_10000_population_comprehensive_hospital_beds)),
    z_total_hospital_beds = as.numeric(scale(value_per_10000_population_hospital_beds)),
    z_log_gdp = as.numeric(scale(log_gdp_per_capita_current_yuan)),
    z_urbanization = as.numeric(scale(urban_population_percent)),
    z_age_65_plus = as.numeric(scale(population_age_65_plus_percent))
  )

if (anyNA(analysis[c(
  "value_per_10000_population_comprehensive_hospital_beds",
  "value_per_10000_population_hospital_beds", "log_gdp_per_capita_current_yuan",
  "urban_population_percent", "population_age_65_plus_percent"
)])) {
  stop("Contextual join produced missing province-year covariates")
}

individual_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile"
)

model_specs <- list(
  C0 = c("z_tcm_beds"),
  C1 = c("z_tcm_beds", "z_comprehensive_beds"),
  C2 = c("z_tcm_beds", "z_log_gdp", "z_urbanization"),
  C3 = c("z_tcm_beds", "z_comprehensive_beds", "z_log_gdp", "z_urbanization"),
  C4 = c("z_tcm_beds", "z_comprehensive_beds", "z_log_gdp", "z_urbanization", "z_age_65_plus"),
  C5 = c("z_tcm_beds", "z_total_hospital_beds", "z_log_gdp", "z_urbanization")
)

model_labels <- c(
  C0 = "Individual covariates",
  C1 = "Add comprehensive-hospital bed density",
  C2 = "Add provincial GDP per capita and urbanization",
  C3 = "Joint health-system and socioeconomic adjustment",
  C4 = "Add province-level population aging",
  C5 = "Alternative adjustment using total hospital beds"
)

term_labels <- c(
  z_tcm_beds = "TCM hospital beds per 10,000 population",
  z_comprehensive_beds = "Comprehensive-hospital beds per 10,000 population",
  z_total_hospital_beds = "All hospital beds per 10,000 population",
  z_log_gdp = "Log provincial GDP per capita",
  z_urbanization = "Urban population share",
  z_age_65_plus = "Population aged 65 years or older"
)

fit_spec <- function(model_id, exposures) {
  rhs <- c(exposures, individual_covariates, "province_fe", "wave_fe")
  formula <- as.formula(paste("outcome_primary ~", paste(rhs, collapse = " + ")))
  needed <- all.vars(formula)
  d <- analysis %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- lm(formula, data = d)
  vcov <- cluster_vcov(fit, d$province_supply_key)
  df <- n_distinct(d$province_supply_key) - 1

  estimates <- bind_rows(lapply(exposures, function(term) {
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
      wild_cluster_pvalue = NA_real_,
      wild_cluster_reps_requested = NA_integer_,
      wild_cluster_reps_valid = NA_integer_,
      wild_cluster_weight = NA_character_,
      wild_cluster_null_imposed = NA,
      n_persons = n_distinct(d$person_id),
      n_person_waves = nrow(d),
      n_clusters = n_distinct(d$province_supply_key),
      fixed_effects = "province and survey year",
      individual_covariates = paste(individual_covariates, collapse = "; "),
      contextual_covariates = paste(setdiff(exposures, "z_tcm_beds"), collapse = "; ")
    )
  }))

  if (model_id %in% c("C1", "C3")) {
    wild <- wild_cluster_pvalue(
      fit, "z_tcm_beds", d$province_supply_key,
      reps = bootstrap_reps,
      seed = bootstrap_seed + match(model_id, names(model_specs))
    )
    estimates[estimates$term == "z_tcm_beds", names(wild)] <- wild[1, ]
  }
  estimates
}

results <- bind_rows(Map(fit_spec, names(model_specs), model_specs)) %>%
  mutate(across(c(effect, ci_lower, ci_upper, cluster_pvalue, wild_cluster_pvalue), ~ round(.x, 5)))

panel_for_diagnostics <- analysis %>%
  distinct(
    province_supply_key, year, z_tcm_beds, z_comprehensive_beds,
    z_total_hospital_beds, z_log_gdp, z_urbanization, z_age_65_plus
  )

diagnose_spec <- function(model_id, exposures) {
  residualized <- lapply(exposures, function(term) {
    residuals(lm(
      as.formula(paste(term, "~ factor(province_supply_key) + factor(year)")),
      data = panel_for_diagnostics
    ))
  })
  x <- do.call(cbind, residualized)
  colnames(x) <- exposures
  x <- scale(x)
  correlation <- cor(x)
  eigenvalues <- eigen(crossprod(x), symmetric = TRUE, only.values = TRUE)$values
  positive <- eigenvalues[eigenvalues > max(eigenvalues) * .Machine$double.eps]
  condition_index <- sqrt(max(positive) / min(positive))

  bind_rows(lapply(seq_along(exposures), function(j) {
    other <- setdiff(seq_along(exposures), j)
    vif <- if (length(other) == 0) {
      1
    } else {
      1 / (1 - summary(lm(x[, j] ~ x[, other, drop = FALSE]))$r.squared)
    }
    tibble(
      model = model_id,
      term = exposures[[j]],
      term_label = unname(term_labels[[exposures[[j]]]]),
      two_way_fe_vif = vif,
      maximum_absolute_pairwise_correlation = if (length(other) == 0) 0 else max(abs(correlation[j, other])),
      design_condition_index = condition_index,
      province_years = nrow(panel_for_diagnostics),
      provinces = n_distinct(panel_for_diagnostics$province_supply_key),
      waves = n_distinct(panel_for_diagnostics$year),
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
print(results %>% filter(term == "z_tcm_beds") %>% select(model, effect, ci_lower, ci_upper, cluster_pvalue, wild_cluster_pvalue))
print(diagnostics %>% filter(term == "z_tcm_beds") %>% select(model, two_way_fe_vif, maximum_absolute_pairwise_correlation, design_condition_index))
