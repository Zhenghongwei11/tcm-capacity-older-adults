#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

image_tsv <- argval(
  "--image-transcription",
  "data/metadata/nbs_province_covariates_image_transcription.tsv"
)
out_tsv <- argval(
  "--output",
  "results/tcm/province_year_socioeconomic_covariates.tsv"
)

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

province_map <- tribble(
  ~province_en, ~province,
  "Beijing", "北京市",
  "Tianjin", "天津市",
  "Hebei", "河北省",
  "Shanxi", "山西省",
  "Inner Mongolia", "内蒙古自治区",
  "Liaoning", "辽宁省",
  "Jilin", "吉林省",
  "Heilongjiang", "黑龙江省",
  "Shanghai", "上海市",
  "Jiangsu", "江苏省",
  "Zhejiang", "浙江省",
  "Anhui", "安徽省",
  "Fujian", "福建省",
  "Jiangxi", "江西省",
  "Shandong", "山东省",
  "Henan", "河南省",
  "Hubei", "湖北省",
  "Hunan", "湖南省",
  "Guangdong", "广东省",
  "Guangxi", "广西壮族自治区",
  "Hainan", "海南省",
  "Chongqing", "重庆市",
  "Sichuan", "四川省",
  "Guizhou", "贵州省",
  "Yunnan", "云南省",
  "Tibet", "西藏自治区",
  "Shaanxi", "陕西省",
  "Gansu", "甘肃省",
  "Qinghai", "青海省",
  "Ningxia", "宁夏回族自治区",
  "Xinjiang", "新疆维吾尔自治区"
)

read_legacy_xls <- function(path) {
  direct <- tryCatch(
    read_excel(path, col_names = FALSE, .name_repair = "minimal"),
    error = function(e) NULL
  )
  if (!is.null(direct)) return(direct)

  if (Sys.which("soffice") == "") {
    stop("LibreOffice soffice is required to read legacy NBS file: ", path)
  }
  conversion_dir <- tempfile("nbs_xlsx_")
  dir.create(conversion_dir)
  status <- system2(
    "soffice",
    c("--headless", "--convert-to", "xlsx", "--outdir", conversion_dir, path),
    stdout = FALSE,
    stderr = FALSE
  )
  converted <- file.path(
    conversion_dir,
    paste0(tools::file_path_sans_ext(basename(path)), ".xlsx")
  )
  if (status != 0 || !file.exists(converted)) stop("Failed to convert ", path)
  read_excel(converted, col_names = FALSE, .name_repair = "minimal")
}

province_rows <- function(raw) {
  names(raw) <- paste0("v", seq_len(ncol(raw)))
  raw %>%
    mutate(province_en = as.character(v1)) %>%
    inner_join(province_map, by = "province_en")
}

find_year_column <- function(raw, year) {
  probe <- raw[seq_len(min(10, nrow(raw))), , drop = FALSE]
  hit <- which(as.matrix(probe) == as.character(year), arr.ind = TRUE)
  if (nrow(hit) == 0) {
    numeric_probe <- suppressWarnings(matrix(as.numeric(as.matrix(probe)), nrow = nrow(probe)))
    hit <- which(numeric_probe == year, arr.ind = TRUE)
  }
  columns <- unique(hit[, "col"])
  if (length(columns) != 2 && year %in% c(2011L, 2013L)) {
    # GDP tables repeat year labels for level and index blocks; urbanization tables do not.
    columns <- unique(hit[, "col"])
  }
  columns
}

urban_gdp_2011_2013 <- function(year) {
  urban_path <- "data/raw/province_covariates/nbs/2014/Z0206E.xls"
  gdp_path <- "data/raw/province_covariates/nbs/2014/Z0315E.xls"
  urban_raw <- read_legacy_xls(urban_path)
  gdp_raw <- read_legacy_xls(gdp_path)

  urban_cols <- find_year_column(urban_raw, year)
  gdp_cols <- find_year_column(gdp_raw, year)
  if (length(urban_cols) != 1 || length(gdp_cols) < 1) {
    stop("Could not identify NBS columns for ", year)
  }

  urban <- province_rows(urban_raw) %>%
    transmute(
      province,
      urban_population_percent = as.numeric(.data[[paste0("v", urban_cols[[1]])]])
    )
  gdp <- province_rows(gdp_raw) %>%
    transmute(
      province,
      gdp_per_capita_current_yuan = as.numeric(.data[[paste0("v", gdp_cols[[1]])]])
    )
  inner_join(gdp, urban, by = "province") %>% mutate(year = as.integer(year))
}

age_2011_2013 <- function(year) {
  path <- if (year == 2011L) {
    "data/raw/province_covariates/nbs/2012/D0311e.xls"
  } else {
    "data/raw/province_covariates/nbs/2014/Z0211E.xls"
  }
  raw <- province_rows(read_legacy_xls(path))
  raw %>%
    transmute(
      province,
      age_composition_sample_population = as.numeric(v2),
      age_composition_sample_age_65_plus = as.numeric(v5),
      year = as.integer(year)
    )
}

tabular_years <- map_dfr(c(2011L, 2013L), function(year) {
  inner_join(
    urban_gdp_2011_2013(year),
    age_2011_2013(year),
    by = c("province", "year")
  ) %>%
    mutate(extraction_method = "official_nbs_excel_parsed")
})

image_years <- read_tsv(image_tsv, show_col_types = FALSE) %>%
  mutate(year = as.integer(year))

panel <- bind_rows(tabular_years, image_years) %>%
  mutate(
    population_age_65_plus_percent =
      age_composition_sample_age_65_plus / age_composition_sample_population * 100,
    log_gdp_per_capita_current_yuan = log(gdp_per_capita_current_yuan),
    gdp_source_id = case_when(
      year %in% c(2011L, 2013L) ~ "nbs_yearbook_2014_z0315e",
      year == 2015L ~ "nbs_yearbook_2016_0310en",
      year == 2018L ~ "nbs_yearbook_2019_e0309"
    ),
    gdp_table_id = case_when(
      year %in% c(2011L, 2013L) ~ "3-15",
      year == 2015L ~ "3-10",
      year == 2018L ~ "3-9"
    ),
    gdp_source_url = case_when(
      year %in% c(2011L, 2013L) ~ "https://www.stats.gov.cn/sj/ndsj/2014/zk/html/Z0315E.xls",
      year == 2015L ~ "https://www.stats.gov.cn/sj/ndsj/2016/html/0310EN.jpg",
      year == 2018L ~ "https://www.stats.gov.cn/sj/ndsj/2019/html/E0309.jpg"
    ),
    urbanization_source_id = case_when(
      year %in% c(2011L, 2013L) ~ "nbs_yearbook_2014_z0206e",
      year == 2015L ~ "nbs_yearbook_2016_0207en",
      year == 2018L ~ "nbs_yearbook_2019_e0207"
    ),
    urbanization_table_id = if_else(year %in% c(2011L, 2013L), "2-6", "2-7"),
    urbanization_source_url = case_when(
      year %in% c(2011L, 2013L) ~ "https://www.stats.gov.cn/sj/ndsj/2014/zk/html/Z0206E.xls",
      year == 2015L ~ "https://www.stats.gov.cn/sj/ndsj/2016/html/0207EN.jpg",
      year == 2018L ~ "https://www.stats.gov.cn/sj/ndsj/2019/html/E0207.jpg"
    ),
    aging_source_id = case_when(
      year == 2011L ~ "nbs_yearbook_2012_d0311e",
      year == 2013L ~ "nbs_yearbook_2014_z0211e",
      year == 2015L ~ "nbs_yearbook_2016_0212en",
      year == 2018L ~ "nbs_yearbook_2019_e0212"
    ),
    aging_table_id = case_when(
      year == 2011L ~ "3-11",
      year == 2013L ~ "2-11",
      TRUE ~ "2-12"
    ),
    aging_source_url = case_when(
      year == 2011L ~ "https://www.stats.gov.cn/sj/ndsj/2012/html/D0311e.xls",
      year == 2013L ~ "https://www.stats.gov.cn/sj/ndsj/2014/zk/html/Z0211E.xls",
      year == 2015L ~ "https://www.stats.gov.cn/sj/ndsj/2016/html/0212EN.jpg",
      year == 2018L ~ "https://www.stats.gov.cn/sj/ndsj/2019/html/E0212.jpg"
    ),
    quality_flag = if_else(
      extraction_method == "official_nbs_excel_parsed",
      "official_excel_parsed",
      "official_image_transcribed_twice"
    ),
    aging_measure_notes = paste(
      "Province age structure is estimated from the NBS annual population sample survey for the corresponding year;",
      "sampling fractions differ across waves, so this variable is reserved for sensitivity adjustment."
    ),
    gdp_measure_notes = paste(
      "Per capita gross regional product at current prices; models include wave fixed effects and use the log value."
    ),
    urbanization_measure_notes = "Proportion of year-end usual-resident population living in urban areas (%)."
  ) %>%
  select(
    province, year, gdp_per_capita_current_yuan, log_gdp_per_capita_current_yuan,
    urban_population_percent, age_composition_sample_population,
    age_composition_sample_age_65_plus, population_age_65_plus_percent,
    gdp_source_id, gdp_table_id, gdp_source_url, gdp_measure_notes,
    urbanization_source_id, urbanization_table_id, urbanization_source_url,
    urbanization_measure_notes, aging_source_id, aging_table_id,
    aging_source_url, aging_measure_notes, extraction_method, quality_flag
  ) %>%
  arrange(year, match(province, province_map$province))

if (nrow(panel) != 124 || anyDuplicated(panel[c("province", "year")])) {
  stop("Expected a unique 31-province panel for four years")
}
if (anyNA(panel)) stop("Unexpected missing value in province-year socioeconomic panel")
if (any(panel$urban_population_percent <= 0 | panel$urban_population_percent >= 100)) {
  stop("Urbanization outside plausible range")
}
if (any(panel$population_age_65_plus_percent <= 0 | panel$population_age_65_plus_percent >= 30)) {
  stop("Age 65+ share outside plausible range")
}
if (any(panel$gdp_per_capita_current_yuan < 5000)) stop("GDP per capita outside plausible range")

write_tsv(panel, out_tsv)
cat("Wrote ", out_tsv, "\n", sep = "")
print(panel %>% group_by(year) %>% summarise(
  provinces = n(),
  median_gdp_per_capita = median(gdp_per_capita_current_yuan),
  mean_urban_population_percent = mean(urban_population_percent),
  mean_age_65_plus_percent = mean(population_age_65_plus_percent),
  .groups = "drop"
))
