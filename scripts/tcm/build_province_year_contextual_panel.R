#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

tcm_path <- argval("--tcm", "results/tcm/province_year_tcm_core_density.tsv")
health_path <- argval("--health", "results/tcm/province_year_health_system_covariates.tsv")
socio_path <- argval("--socioeconomic", "results/tcm/province_year_socioeconomic_covariates.tsv")
out_path <- argval("--output", "results/tcm/province_year_contextual_covariates.tsv")
corr_path <- argval(
  "--correlations",
  "results/tcm/province_year_contextual_covariate_correlations.tsv"
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

tcm <- read_tsv(tcm_path, show_col_types = FALSE) %>%
  select(
    province, year, denominator_population,
    value_tcm_hospital_beds,
    value_tcm_physicians,
    value_per_10000_population_tcm_hospital_beds,
    value_per_10000_population_tcm_physicians,
    tcm_source_id = source, tcm_source_url = source_url,
    tcm_table_id = table_id
  )

health <- read_tsv(health_path, show_col_types = FALSE) %>%
  select(
    province, year, health_denominator_population = denominator_population,
    value_total_health_institution_beds,
    value_hospital_beds,
    value_comprehensive_hospital_beds,
    value_per_10000_population_total_health_institution_beds,
    value_per_10000_population_hospital_beds,
    value_per_10000_population_comprehensive_hospital_beds,
    health_source_id = source_id, health_source_url = source_url,
    health_table_id = table_id, health_quality_flag = quality_flag
  )

socio <- read_tsv(socio_path, show_col_types = FALSE)

panel <- tcm %>%
  inner_join(health, by = c("province", "year")) %>%
  inner_join(socio, by = c("province", "year")) %>%
  mutate(
    denominator_match = denominator_population == health_denominator_population
  ) %>%
  arrange(year, province)

if (nrow(panel) != 124 || anyDuplicated(panel[c("province", "year")])) {
  stop("Expected a unique 31-province contextual panel for four years")
}
if (anyNA(panel)) stop("Unexpected missing value in contextual panel")
if (!all(panel$denominator_match)) stop("Health-system and TCM density denominators differ")

numeric_covariates <- c(
  "value_per_10000_population_tcm_hospital_beds",
  "value_per_10000_population_tcm_physicians",
  "value_per_10000_population_comprehensive_hospital_beds",
  "value_per_10000_population_hospital_beds",
  "log_gdp_per_capita_current_yuan",
  "urban_population_percent",
  "population_age_65_plus_percent"
)

panel <- panel %>%
  mutate(across(
    all_of(numeric_covariates),
    ~ as.numeric(scale(.x)),
    .names = "z_{.col}"
  ))

long <- panel %>%
  select(province, year, all_of(numeric_covariates)) %>%
  pivot_longer(all_of(numeric_covariates), names_to = "variable", values_to = "value") %>%
  group_by(province, variable) %>%
  mutate(within_value = value - mean(value)) %>%
  ungroup() %>%
  group_by(year, variable) %>%
  mutate(year_mean = mean(value)) %>%
  ungroup() %>%
  group_by(province, variable) %>%
  mutate(province_mean = mean(value)) %>%
  ungroup() %>%
  group_by(variable) %>%
  mutate(
    grand_mean = mean(value),
    two_way_fe_value = value - province_mean - year_mean + grand_mean
  ) %>%
  ungroup()

correlations <- tidyr::crossing(
  variable_1 = numeric_covariates,
  variable_2 = numeric_covariates
) %>%
  rowwise() %>%
  mutate(
    raw_correlation = cor(
      long$value[long$variable == variable_1],
      long$value[long$variable == variable_2]
    ),
    within_province_correlation = cor(
      long$within_value[long$variable == variable_1],
      long$within_value[long$variable == variable_2]
    ),
    two_way_fe_residual_correlation = cor(
      long$two_way_fe_value[long$variable == variable_1],
      long$two_way_fe_value[long$variable == variable_2]
    ),
    n_province_years = nrow(panel)
  ) %>%
  ungroup()

write_tsv(panel, out_path)
write_tsv(correlations, corr_path)

cat("Wrote ", out_path, "\n", sep = "")
cat("Wrote ", corr_path, "\n", sep = "")
print(panel %>% count(year, name = "province_years"))
print(correlations %>%
  filter(
    variable_1 == "value_per_10000_population_tcm_hospital_beds",
    variable_2 != variable_1
  ) %>%
  select(
    variable_2, raw_correlation, within_province_correlation,
    two_way_fe_residual_correlation
  ))
