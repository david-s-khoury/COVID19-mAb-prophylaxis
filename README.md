# COVID-19 mAb prophylaxis

This repository provides the data and code to reproduce the analysis for COVID-19 mAb prophylaxis. 

Key features include:

- `driver.R` which can run the analysis from start to finish
- `setup.R` which loads all libraries and other variables used throughout
- `.gitignore` which list files and directories which should be ignored by `git`
- `raw-data/` which holds raw-data files. This directory can be added to the `.gitignore` file for sensitive files that should not be tracked
- `processing/` which hold scripts to clean and process data before analysis
- `analysis/` which hold scripts used to analyse data and generate outputs
- `output/` which holds outputs. This directory should be added to the `.gitignore` so outputs are not stored in the `git` history
- `helper-functions` which holds any custom functions used throughout the project

