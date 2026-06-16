# glmmPackCompare

This repository contains R scripts for the paper "A Comparison of `R` Packages for Estimating Generalized Linear Mixed Models". 
The preprint link is https://arxiv.org/abs/2606.15933.

## Environmental requirements
`R` (>=4.0) with the following packages: `lme4`, `GLMMadaptive`, `glmmTMB`, `MASS`, `hglm`, `brms`, `rstanarm`, `loo`, `jsonlite`.

## Usage
Run the scripts in folder `scripts/` **in numerical order** from an `R` session or from the command line. Adjust parameters as needed before running.

### Step 0 – Set hyperparameters and generate configuration files
```r
source("scripts/00_SetHyperparameters.R")
```
This creates JSON files `config_1.json` … `config_24.json` inside the `config/` folder. 
Each file defines a combination of:

- `n`: sample size (50, 100, 200)

- `M`: maximum number of repeated measurements (3, 9)

- `re_dist`: response distribution (`'Bernoulli'` or `'Poisson'`)

- `with_rd_slope`: random‑intercept only (`FALSE`) or random‑intercept‑plus‑slope (`TRUE`)

### Step 1 - Generate simulated datasets
```r
source("scripts/01_GenerateData.R")
```
Reads each configuration file and create simulated datasets.
The data are saved as `.RData` files in the `data/` folder, named e.g., `n50_M3_Bernoulli_rdi.RData`.

### Step 2 - Estimate models
```r
source('scripts/02_Estimate.R')
```
Before running, you must edit the script to specify:

- `c` - which configuration file to run (default is 15). To run all configurations, you would need to wrap the script in a loop or submit array jobs.
- `n_sim` - which specify the number of simulated datasets to test. To test all simulated datasets, you would specify `n_sim = length(simulate_data_ls)`.

The experimental results are saved into `results/` (e.g., `results/n50_M3_Bernoulli_rdi.RData`).

### Step 3 - Gather results
```r
source('scripts/03_GatherResults.R')
```
The final output is saved in `results/final.RData`, which includes:
- `final_result_df` - detailed results for each configuration and package
- `mean_abs_bias_df` - average absolute bias per package across all configurations
- `mean_rmse_df` - average RMSE per package across all configurations

## Contact
For questions, please open an issue in the repository or contact the author via x.li@math.leidenuniv.nl.