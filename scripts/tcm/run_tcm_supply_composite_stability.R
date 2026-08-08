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

analysis_tsv <- argval("--analysis", "results/tcm/person_wave_tcm_core_density_analysis.tsv")
out_tsv <- argval("--output", "results/tcm/tcm_supply_composite_stability.tsv")

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
    z_tcm_beds_per_10000 = z_py_tcm_beds_per_10000
  )

component_vars <- grep("^condition_tcm_", names(d), value = TRUE)

cluster_vcov <- function(model, cluster) {
  cluster <- as.factor(cluster)
  x <- model.matrix(model)
  u <- residuals(model)
  bread <- solve(crossprod(x))
  meat <- matrix(0, ncol(x), ncol(x))
  for (g in levels(cluster)) {
    idx <- cluster == g
    score <- crossprod(x[idx, , drop = FALSE], u[idx])
    meat <- meat + tcrossprod(score)
  }
  n <- nrow(x)
  k <- ncol(x)
  g <- nlevels(cluster)
  (g / (g - 1)) * ((n - 1) / (n - k)) * bread %*% meat %*% bread
}

fit_omission <- function(data, omitted, scenario) {
  retained <- setdiff(component_vars, omitted)
  dd <- data %>%
    mutate(stability_outcome = as.numeric(rowSums(across(all_of(retained)), na.rm = TRUE) > 0))
  f <- as.formula(paste("stability_outcome ~", paste(c("z_tcm_beds_per_10000", full_covariates, "province_fe", "wave_fe"), collapse = " + ")))
  needed <- all.vars(f)
  dd <- dd %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- lm(f, data = dd)
  vc <- cluster_vcov(fit, dd$province_supply_key)
  raw_effect <- unname(coef(fit)[["z_tcm_beds_per_10000"]])
  se <- unname(sqrt(diag(vc))[["z_tcm_beds_per_10000"]])
  df <- n_distinct(dd$province_supply_key) - 1
  ci <- raw_effect + c(-1, 1) * qt(0.975, df = df) * se
  tibble(
    scenario = scenario,
    omitted_components = ifelse(length(omitted) == 0, "none", paste(omitted, collapse = "; ")),
    retained_component_count = length(retained),
    outcome_prevalence_percent = 100 * mean(dd$stability_outcome),
    effect_type = "percentage-point difference per 1 SD higher province-year TCM hospital bed density",
    effect = 100 * raw_effect,
    ci_lower = 100 * ci[[1]],
    ci_upper = 100 * ci[[2]],
    pvalue = 2 * pt(abs(raw_effect / se), df = df, lower.tail = FALSE),
    n_persons = n_distinct(dd$person_id),
    n_person_waves = nrow(dd),
    n_clusters = n_distinct(dd$province_supply_key),
    fixed_effects = "province and survey year",
    source_table = analysis_tsv
  )
}

family_omissions <- list(
  omit_cardiometabolic = intersect(component_vars, c(
    "condition_tcm_dyslipidemia", "condition_tcm_heart", "condition_tcm_hypertension",
    "condition_tcm_diabetes", "condition_tcm_stroke_chinese_medicine",
    "condition_tcm_stroke_acupuncture_moxibustion"
  )),
  omit_musculoskeletal = intersect(component_vars, "condition_tcm_arthritis"),
  omit_digestive_hepatorenal = intersect(component_vars, c(
    "condition_tcm_liver", "condition_tcm_kidney", "condition_tcm_digestive"
  ))
)

results <- bind_rows(
  fit_omission(d, character(), "primary_composite_reproduction"),
  bind_rows(lapply(component_vars, function(v) fit_omission(d, v, paste0("leave_one_out_", sub("^condition_tcm_", "", v))))),
  bind_rows(lapply(names(family_omissions), function(nm) fit_omission(d, family_omissions[[nm]], nm)))
) %>% mutate(across(where(is.numeric), ~ round(.x, 5)))

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results)
