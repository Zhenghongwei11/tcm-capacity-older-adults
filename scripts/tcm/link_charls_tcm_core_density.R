#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

outcomes_tsv <- argval("--outcomes", "results/tcm/person_wave_tcm_outcomes.tsv")
supply_tsv <- argval("--supply", "results/tcm/province_year_tcm_core_density.tsv")
out_tsv <- argval("--output", "results/tcm/person_wave_tcm_core_density_linked.tsv")

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

outcomes <- read_tsv(outcomes_tsv, show_col_types = FALSE)
supply <- read_tsv(supply_tsv, show_col_types = FALSE)

required_outcome_cols <- c("ID", "communityID", "province", "year")
required_supply_cols <- c(
  "province", "year",
  "value_tcm_hospital_beds",
  "value_tcm_practicing_assistant_physicians",
  "value_per_10000_population_tcm_hospital_beds",
  "value_per_10000_population_tcm_practicing_assistant_physicians",
  "denominator_population", "quality_flag", "resource_definition"
)

missing_outcome_cols <- setdiff(required_outcome_cols, names(outcomes))
missing_supply_cols <- setdiff(required_supply_cols, names(supply))
if (length(missing_outcome_cols) > 0) {
  stop(sprintf("Outcome table missing columns: %s", paste(missing_outcome_cols, collapse = ", ")))
}
if (length(missing_supply_cols) > 0) {
  stop(sprintf("Supply table missing columns: %s", paste(missing_supply_cols, collapse = ", ")))
}

supply_for_join <- supply %>%
  mutate(province_supply_key = province) %>%
  rename(
    supply_province = province,
    core_density_source = source,
    core_density_source_name = source_name,
    core_density_source_url = source_url,
    core_density_source_edition = source_edition,
    core_density_extraction_status = extraction_status,
    core_density_quality_flag = quality_flag,
    core_density_notes = notes
  )

outcomes_for_join <- outcomes %>%
  mutate(
    province_supply_key = recode(
      province,
      "北京" = "北京市",
      "天津" = "天津市",
      "广西省" = "广西壮族自治区",
      .default = province
    )
  )

linked <- outcomes_for_join %>%
  left_join(supply_for_join, by = c("province_supply_key", "year")) %>%
  mutate(
    core_density_link_status = case_when(
      !is.na(value_per_10000_population_tcm_practicing_assistant_physicians) ~ "linked_core_density_ready",
      TRUE ~ "missing_core_density_supply"
    ),
    main_core_density_linked_panel = year %in% c(2011L, 2013L, 2015L, 2018L) &
      core_density_link_status == "linked_core_density_ready"
  )

write_tsv(linked, out_tsv)

qc <- linked %>%
  count(year, core_density_link_status, name = "n_person_waves")

print(qc)
cat(sprintf("Wrote %s\n", out_tsv))
