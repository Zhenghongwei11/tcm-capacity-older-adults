#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

out_tsv <- argval(
  "--output",
  "results/tcm/province_year_tcm_core_density.tsv"
)

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

province_map <- tibble::tribble(
  ~province_raw, ~province,
  "北京市", "北京市",
  "天津市", "天津市",
  "河北省", "河北省",
  "山西省", "山西省",
  "内蒙古自治区", "内蒙古自治区",
  "辽宁省", "辽宁省",
  "吉林省", "吉林省",
  "黑龙江省", "黑龙江省",
  "上海市", "上海市",
  "江苏省", "江苏省",
  "浙江省", "浙江省",
  "安徽省", "安徽省",
  "福建省", "福建省",
  "江西省", "江西省",
  "山东省", "山东省",
  "河南省", "河南省",
  "湖北省", "湖北省",
  "湖南省", "湖南省",
  "广东省", "广东省",
  "广西壮族自治区", "广西壮族自治区",
  "海南省", "海南省",
  "重庆市", "重庆市",
  "四川省", "四川省",
  "贵州省", "贵州省",
  "云南省", "云南省",
  "西藏自治区", "西藏自治区",
  "陕西省", "陕西省",
  "甘肃省", "甘肃省",
  "青海省", "青海省",
  "宁夏回族自治区", "宁夏回族自治区",
  "新疆维吾尔自治区", "新疆维吾尔自治区"
)

sources <- tibble::tribble(
  ~year, ~table_id, ~path, ~url, ~expected_beds, ~expected_physicians,
  2011L, "A86", "data/raw/tcm_supply/official_excels/natcm_2011_extract/A86.xlsx",
  "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2011/A86.xlsx",
  529349, 309272,
  2013L, "A95", "data/raw/tcm_supply/official_excels/natcm_core_density/A95_2013.xls",
  "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2013/A95.xls",
  686793, 381682,
  2015L, "A95", "data/raw/tcm_supply/official_excels/natcm_core_density/A95_2015.xls",
  "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2015/A95.xls",
  819412, 452190,
  2018L, "A95", "data/raw/tcm_supply/official_excels/natcm_core_density/A95_2018.xlsx",
  "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2018/A95.xlsx",
  1021548, 575454
)

clean_province <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("　", "") %>%
    str_trim()
}

extract_one <- function(year, table_id, path, url, expected_beds, expected_physicians) {
  if (!file.exists(path)) {
    stop(sprintf("Missing NATCM source file for %s: %s", year, path))
  }

  raw <- read_excel(path, col_names = FALSE)

  out <- raw %>%
    transmute(
      province_raw = clean_province(...1),
      population_10k = suppressWarnings(as.numeric(...2)),
      value_tcm_hospital_beds = suppressWarnings(as.numeric(...3)),
      value_per_10000_population_tcm_hospital_beds = suppressWarnings(as.numeric(...4)),
      rank_tcm_hospital_beds_per_10000 = suppressWarnings(as.integer(...5)),
      value_tcm_physicians = suppressWarnings(as.numeric(...6)),
      value_per_10000_population_tcm_physicians = suppressWarnings(as.numeric(...7)),
      rank_tcm_physicians_per_10000 = suppressWarnings(as.integer(...8))
    ) %>%
    inner_join(province_map, by = "province_raw") %>%
    transmute(
      province,
      year = as.integer(year),
      source = paste0("natcm_", year, "_extract_", tolower(table_id)),
      source_name = paste0(
        "National Administration of Traditional Chinese Medicine ",
        year,
        " National TCM Statistical Extract"
      ),
      source_url = url,
      source_edition = as.character(year),
      table_id = table_id,
      table_title = paste0(
        year,
        "年各地区万人口中医类医院床位数及万人口全国中医执业(助理)医师数"
      ),
      extraction_status = "official_excel_numeric_extracted",
      quality_flag = "official_natcm_excel_parsed",
      resource_definition = "TCM hospital beds and TCM physicians per population.",
      denominator_population = population_10k * 10000,
      denominator_unit = "persons",
      denominator_source_id = paste0("natcm_", year, "_extract_", tolower(table_id)),
      denominator_notes = "NATCM table reports province population in 10,000 persons; converted to persons.",
      value_tcm_hospital_beds,
      value_tcm_physicians,
      value_per_10000_population_tcm_hospital_beds,
      value_per_10000_population_tcm_physicians,
      rank_tcm_hospital_beds_per_10000,
      rank_tcm_physicians_per_10000,
      unit_tcm_hospital_beds = "beds",
      unit_tcm_physicians = "persons",
      notes = "Official NATCM Excel table; main core-density supply exposure."
    ) %>%
    arrange(province)

  if (nrow(out) != 31) {
    stop(sprintf("Expected 31 provinces for %s, got %s", year, nrow(out)))
  }
  if (sum(out$value_tcm_hospital_beds) != expected_beds) {
    stop(sprintf("Unexpected bed total for %s", year))
  }
  if (sum(out$value_tcm_physicians) != expected_physicians) {
    stop(sprintf("Unexpected physician total for %s", year))
  }
  out
}

panel <- bind_rows(purrr::pmap(sources, extract_one)) %>%
  arrange(year, province)

stopifnot(nrow(panel) == 124)
stopifnot(!any(duplicated(panel[c("province", "year")])))

write_tsv(panel, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(panel %>% count(year, source, table_id))
print(panel %>% group_by(year) %>% summarise(
  provinces = n(),
  beds_total = sum(value_tcm_hospital_beds),
  physicians_total = sum(value_tcm_physicians),
  .groups = "drop"
))
