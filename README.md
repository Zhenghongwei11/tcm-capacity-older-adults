# Traditional Chinese Medicine Service Capacity and Realized Use Among Older Adults in China

This repository contains code, derived aggregate tables, and figure-generation files for a longitudinal study linking CHARLS 2011, 2013, 2015, and 2018 survey waves with official province-year traditional Chinese medicine (TCM) service-resource data.

## Study Scope

The analysis examines whether provincial TCM service capacity was reflected in reported TCM treatment use among older adults in China before COVID-19. The public repository is designed to reproduce the reported tables and figures after authorized users place the restricted CHARLS files in the expected local directory.

## Data Access

CHARLS microdata are not redistributed here. Users must obtain CHARLS regular-wave and Harmonized CHARLS files through the official data portals and place them locally as described in `data/README.md`.

Official province-year TCM service-resource tables are public National Administration of Traditional Chinese Medicine statistical extracts. Derived provincial supply tables are included under `results/tcm/`.

## Quick Start

Install the R packages listed in `environment.yml`, then run:

```bash
bash scripts/reproduce_one_click.sh
```

The script checks for the expected CHARLS inputs, rebuilds the official supply and province-context tables, reruns the longitudinal models, and regenerates reported tables and publication figures.

## Repository Contents

- `scripts/`: analysis, model, table, and figure scripts.
- `results/tcm/`: aggregate derived tables and model summaries used by the figures and tables.
- `results/tcm/figure_source_data/`: source data for each figure panel.
- `plots/publication/`: current figure exports.
- `tables/`: current reported table set.
- `docs/`: public data manifest, figure provenance, table provenance, statistical decision rules, and compute notes.
- `data/`: access instructions and small public metadata tables used by extraction scripts.

## Reproduction Boundary

The repository excludes raw CHARLS files, row-level person-wave analytic tables, protected geography, and local-only linkage files. Reproducing individual-level model estimates requires access to CHARLS under the original data-use terms.

## Citation

If using this repository, cite the archived version 1.0.0 release at https://doi.org/10.5281/zenodo.21622933 and cite CHARLS and the official statistical sources according to their source requirements.
