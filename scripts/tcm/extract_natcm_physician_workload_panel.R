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
  "results/tcm/province_year_tcm_physician_workload.tsv"
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
  ~year, ~table_id, ~path, ~url, ~table_title,
  2011L, "B35", "data/raw/tcm_supply/official_excels/natcm_service_volume/B35_2011.xls",
  "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2011/B35.xls",
  "2011年政府办中医类医院按地区分平均每一职工、医师产出情况",
  2013L, "B38", "data/raw/tcm_supply/official_excels/natcm_service_volume/B38_2013.xls",
  "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2013/B38.xls",
  "2013年政府办中医类医院按地区分医院医师工作效率",
  2015L, "B38", "data/raw/tcm_supply/official_excels/natcm_service_volume/B38_2015.xls",
  "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2015/B38.xls",
  "2015年政府办中医类医院按地区分医院医师工作效率",
  2018L, "B38", "data/raw/tcm_supply/official_excels/natcm_service_volume/B38_2018.xls",
  "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2018/B38.xls",
  "2018年政府办中医类医院按地区分医院医师工作效率"
)

clean_province <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("　", "") %>%
    str_trim()
}

extract_one <- function(year, table_id, path, url, table_title) {
  if (!file.exists(path)) {
    stop(sprintf("Missing NATCM service-volume source for %s: %s", year, path))
  }

  raw <- read_excel(path, col_names = FALSE)
  names(raw) <- paste0("v", seq_len(ncol(raw)))

  out <- raw %>%
    transmute(
      province_raw = clean_province(v1),
      physician_annual_outpatient_emergency_visits = suppressWarnings(as.numeric(v2)),
      physician_annual_inpatient_bed_days = suppressWarnings(as.numeric(v3)),
      physician_daily_outpatient_emergency_visits = suppressWarnings(as.numeric(v4)),
      physician_daily_inpatient_bed_days = suppressWarnings(as.numeric(v5))
    ) %>%
    inner_join(province_map, by = "province_raw") %>%
    transmute(
      province,
      year = as.integer(year),
      source = paste0("natcm_", year, "_", tolower(table_id)),
      source_name = paste0(
        "National Administration of Traditional Chinese Medicine ",
        year,
        " National TCM Statistical Extract"
      ),
      source_url = url,
      source_edition = as.character(year),
      table_id = table_id,
      table_title = table_title,
      extraction_status = "official_excel_numeric_extracted",
      quality_flag = "official_natcm_excel_parsed",
      institution_scope = "government-run TCM-category hospitals",
      mechanism_role = "province-level physician workload and service-intensity mediator candidate",
      physician_annual_outpatient_emergency_visits,
      physician_annual_inpatient_bed_days,
      physician_daily_outpatient_emergency_visits,
      physician_daily_inpatient_bed_days,
      notes = paste0(
        "This table does not measure patient-level TCM use; it describes whether ",
        "registered government-run TCM hospital physicians translate headcount capacity ",
        "into outpatient/emergency and inpatient service workload."
      )
    ) %>%
    arrange(province)

  if (nrow(out) != 31) {
    stop(sprintf("Expected 31 provinces for %s, got %s", year, nrow(out)))
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
  median_annual_outpatient_emergency_visits =
    median(physician_annual_outpatient_emergency_visits, na.rm = TRUE),
  median_annual_inpatient_bed_days =
    median(physician_annual_inpatient_bed_days, na.rm = TRUE),
  .groups = "drop"
))
