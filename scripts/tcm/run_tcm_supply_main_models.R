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
  "results/tcm/tcm_supply_main_models.tsv"
)
bootstrap_reps <- as.integer(argval("--bootstrap-reps", "9999"))
bootstrap_seed <- as.integer(argval("--bootstrap-seed", "20260720"))

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    household_wealth_quartile = factor(household_wealth_quartile),
    z_tcm_physicians_per_10000 = z_py_tcm_physicians_per_10000,
    z_tcm_beds_per_10000 = z_py_tcm_beds_per_10000,
    outcome_primary = as.numeric(primary_condition_tcm_any)
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
    wild_cluster_pvalue = (1 + sum(abs(boot_t[valid]) >= abs(observed_t))) / (1 + sum(valid)),
    wild_cluster_reps_requested = reps,
    wild_cluster_reps_valid = sum(valid),
    wild_cluster_weight = "Webb six-point",
    wild_cluster_null_imposed = TRUE,
    wild_cluster_seed = seed
  )
}

fit_one <- function(data, supply_var, model_id, covariates, run_wild = FALSE) {
  rhs <- c(supply_var, covariates, "province_fe", "wave_fe")
  f <- as.formula(paste("outcome_primary ~", paste(rhs, collapse = " + ")))
  needed <- all.vars(f)
  d <- data %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- lm(f, data = d)
  vc <- cluster_vcov(fit, d$province_supply_key)
  coefs <- coef(fit)
  se <- sqrt(diag(vc))
  term <- supply_var
  effect <- unname(coefs[[term]])
  se_term <- unname(se[[term]])
  t_value <- effect / se_term
  df <- length(unique(d$province_supply_key)) - 1
  p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
  ci <- effect + c(-1, 1) * qt(0.975, df = df) * se_term

  result <- tibble(
    outcome = "Any disease/condition-specific Chinese medicine treatment use",
    supply_indicator = case_when(
      supply_var == "z_tcm_physicians_per_10000" ~ "TCM practicing/assistant physicians per 10,000 population",
      supply_var == "z_tcm_beds_per_10000" ~ "TCM hospital beds per 10,000 population",
      TRUE ~ supply_var
    ),
    model = model_id,
    effect_type = "percentage-point difference per 1 SD higher province-year supply",
    effect = 100 * effect,
    ci_lower = 100 * ci[[1]],
    ci_upper = 100 * ci[[2]],
    pvalue = p_value,
    wild_cluster_pvalue = NA_real_,
    wild_cluster_reps_requested = NA_integer_,
    wild_cluster_reps_valid = NA_integer_,
    wild_cluster_weight = NA_character_,
    wild_cluster_null_imposed = NA,
    wild_cluster_seed = NA_integer_,
    n_persons = n_distinct(d$person_id),
    n_person_waves = nrow(d),
    cluster_level = "province",
    n_clusters = n_distinct(d$province_supply_key),
    covariate_set = paste(covariates, collapse = "; "),
    fixed_effects = "province and survey year",
    source_table = analysis_tsv
  )

  if (run_wild) {
    wild <- wild_cluster_pvalue(
      fit, term, d$province_supply_key,
      reps = bootstrap_reps,
      seed = bootstrap_seed
    )
    result[names(wild)] <- wild[1, ]
  }
  result
}

minimal_covariates <- c("age_model", "I(age_model^2)", "female", "rural_community")
full_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile"
)

results <- bind_rows(
  fit_one(analysis, "z_tcm_physicians_per_10000", "M1_minimal_adjusted_lpm", minimal_covariates),
  fit_one(analysis, "z_tcm_beds_per_10000", "M1_minimal_adjusted_lpm", minimal_covariates),
  fit_one(analysis, "z_tcm_physicians_per_10000", "M2_covariate_adjusted_lpm", full_covariates),
  fit_one(
    analysis, "z_tcm_beds_per_10000", "M2_covariate_adjusted_lpm",
    full_covariates, run_wild = TRUE
  )
) %>%
  mutate(across(c(effect, ci_lower, ci_upper, pvalue, wild_cluster_pvalue), ~ round(.x, 5)))

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results)
