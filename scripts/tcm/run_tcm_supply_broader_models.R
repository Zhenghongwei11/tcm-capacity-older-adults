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
  "results/tcm/tcm_supply_broader_secondary_models.tsv"
)

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60, year %in% c(2011L, 2013L, 2015L)) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    z_tcm_physicians_per_10000 = z_py_tcm_physicians_per_10000,
    z_tcm_beds_per_10000 = z_py_tcm_beds_per_10000,
    outcome_broader = as.numeric(broader_tcm_use_2011_2015)
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

fit_one <- function(data, supply_var, model_id, covariates) {
  rhs <- c(supply_var, covariates, "province_fe", "wave_fe")
  f <- as.formula(paste("outcome_broader ~", paste(rhs, collapse = " + ")))
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

  tibble(
    outcome = "Broader TCM-related use in 2011-2015",
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

minimal_covariates <- c("age_model", "I(age_model^2)", "female", "rural_community")
full_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile"
)

results <- bind_rows(
  fit_one(analysis, "z_tcm_physicians_per_10000", "B1_minimal_adjusted_lpm", minimal_covariates),
  fit_one(analysis, "z_tcm_beds_per_10000", "B1_minimal_adjusted_lpm", minimal_covariates),
  fit_one(analysis, "z_tcm_physicians_per_10000", "B2_covariate_adjusted_lpm", full_covariates),
  fit_one(analysis, "z_tcm_beds_per_10000", "B2_covariate_adjusted_lpm", full_covariates)
) %>%
  mutate(across(c(effect, ci_lower, ci_upper, pvalue), ~ round(.x, 5)))

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results)
