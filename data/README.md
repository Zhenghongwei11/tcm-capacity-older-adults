# Data Access

## CHARLS

CHARLS microdata are not included. To rerun models from the beginning, obtain the required CHARLS regular-wave files and Harmonized CHARLS Version D through their official access systems, then place them under:

```text
data/raw/charls/downloads/
data/raw/charls/downloads/harmonized/regular_waves/
```

The analysis scripts expect the regular 2011, 2013, 2015, and 2018 waves and the Harmonized CHARLS file used for covariate construction. Restricted geography files are not needed for the province-year analysis.

## Official Province-Year Sources

Public official source URLs and derived non-identifying province-year tables are documented in `docs/DATA_MANIFEST.tsv`.
