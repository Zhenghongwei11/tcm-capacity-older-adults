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
  "results/tcm/tcm_supply_falsification_models.tsv"
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

future_supply <- read_tsv(context_tsv, show_col_types = FALSE) %>%
  transmute(
    province_supply_key = province,
    year = as.integer(year),
    tcm_beds = value_per_10000_population_tcm_hospital_beds
  ) %>%
  arrange(province_supply_key, year) %>%
  group_by(province_supply_key) %>%
  mutate(
    next_year = lead(year),
    future_tcm_beds = lead(tcm_beds)
  ) %>%
  ungroup() %>%
  select(province_supply_key, year, next_year, future_tcm_beds)

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  left_join(future_supply, by = c("province_supply_key", "year")) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    primary = as.numeric(primary_condition_tcm_any),
    outpatient = as.numeric(general_outpatient_visit_last_month),
    inpatient = as.numeric(general_hospital_stay_last_year),
    z_tcm_beds = z_py_tcm_beds_per_10000,
    z_future_tcm_beds = as.numeric(scale(future_tcm_beds)),
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

fit_one <- function(outcome, terms, model_id, check_type, interpretation) {
  rhs <- c(terms, individual_covariates, "province_fe", "wave_fe")
  formula <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
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
      check_type = check_type,
      model = model_id,
      outcome = case_when(
        outcome == "primary" ~ "Any disease/condition-specific Chinese medicine treatment use",
        outcome == "outpatient" ~ "Any doctor visit or outpatient care in the past month",
        outcome == "inpatient" ~ "Any hospital stay in the past year",
        TRUE ~ outcome
      ),
      term = term,
      term_label = case_when(
        term == "z_tcm_beds" ~ "Concurrent TCM hospital beds per 10,000 population",
        term == "z_future_tcm_beds" ~ "Next-wave TCM hospital beds per 10,000 population",
        term == "z_comprehensive_beds" ~ "Comprehensive-hospital beds per 10,000 population",
        term == "z_log_gdp" ~ "Log provincial GDP per capita",
        term == "z_urbanization" ~ "Urban population share",
        TRUE ~ term
      ),
      effect_type = "percentage-point difference per 1 SD higher province-year measure",
      effect = 100 * effect,
      ci_lower = 100 * ci[[1]],
      ci_upper = 100 * ci[[2]],
      cluster_pvalue = pvalue,
      n_persons = n_distinct(d$person_id),
      n_person_waves = nrow(d),
      n_clusters = n_distinct(d$province_supply_key),
      analysis_years = paste(sort(unique(d$year)), collapse = ";"),
      fixed_effects = "province and survey year",
      interpretation_boundary = interpretation
    )
  }))
}

future_terms <- c("z_tcm_beds", "z_future_tcm_beds")
specificity_terms <- c(
  "z_tcm_beds", "z_comprehensive_beds", "z_log_gdp", "z_urbanization"
)

results <- bind_rows(
  fit_one(
    "primary", future_terms, "F1_future_exposure",
    "future-exposure falsification",
    "A next-wave association conditional on concurrent supply would weaken temporal interpretation; absence is supportive but not proof against time-varying confounding."
  ),
  fit_one(
    "outpatient", specificity_terms, "F2_general_outpatient",
    "general-utilization specificity outcome",
    "This all-cause outcome can include TCM care and is not a strict non-TCM negative-control outcome."
  ),
  fit_one(
    "inpatient", specificity_terms, "F3_general_inpatient",
    "general-utilization specificity outcome",
    "This all-cause outcome can include TCM care and is not a strict non-TCM negative-control outcome."
  )
) %>%
  mutate(across(c(effect, ci_lower, ci_upper, cluster_pvalue), ~ round(.x, 5)))

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results %>% filter(term %in% c("z_tcm_beds", "z_future_tcm_beds")) %>% select(check_type, model, outcome, term, effect, ci_lower, ci_upper, cluster_pvalue, n_person_waves))
