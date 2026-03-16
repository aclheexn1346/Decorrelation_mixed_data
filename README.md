## Mixed data de-correlation methodology

This repository implements a **de-correlation methodology for mixed data (continuous and discrete variables)**. The goal of the procedure is to transform correlated observations into approximately **independent latent variables**, which can then be used in downstream tasks such as causal discovery.

The main implementation is contained in: `main.R`


---

## Overview

The method estimates latent variables and a covariance structure for mixed data using an **iterative EM-style procedure**. It handles both continuous and discrete variables by introducing latent Gaussian variables for the discrete observations.

The main function:

```r
mixed_data_decor(data, block_sizes)
```

Primary functions used within **mixed_data_decor** from `helperFuncMixed.R`
- **EM_tau_beta_mixed:**  
  Pre-estimation of thresholds T and coefficients $\beta$.

- **Sig_Estimate_DAG_mixed:**  
  Pre-estimation of row-covariance Sigma.

- **epsilon_draw_mixed:**  
  Draws of epsilon subject to estimated thresholds.

- **obtain_hidden_Z_mixed:**  
  Recovery of latent Z for discrete variables.

A user need only input the dataset and block sizes among the units of the data. For example, if you have a dataset with $n = 20$ observations and every 5 units are dependent, the block_sizes argument should be:

```r
block_sizes = c(5,5,5,5)
```

The algorithm will detect the discrete columns if the number of levels in the discrete variables is < 5. If there are more, please input the discrete columns as a vector of 1's and 0's of the same length as the number of total variables (1 corresponds to a discrete variable, 0 corresponds to a continuous variable). For example, if there are 5 total variables where column 1 and 4 are discrete variables, the discColumns argument should be:

```r
discColumns = c(1,0,0,1,0)
```

---

## Required files

Ensure that the packages are installed within `libraries.R`

```r
source("libraries.R")
source("helperFuncMixed.R")
```

`helperFuncMixed.R` contains the individual functions used in `main.R`

---

## Example usage

Download `test_data.RData`, which contains a mixed variable dataset with $n = 100$ observations and $p = 100$ variables. The example block structure suggests every 5 units are dependent. (Obs. 1 - 5, Obs. 6 - 10, etc. are dependent)

```r
load("test_data.RData")
test_decor = mixed_data_decor(test_data, block_sizes = rep(5,20))
```

---

## Output

The function returns a list with two elements:

- **Element 1:**  
  A list of de-correlated datasets generated during the iterative algorithm.  
  These datasets can be used for the consensus approach.

- **Element 2:**  
  The average of the de-correlated datasets, which can be used for the average approach.

---

## Single-cell RNA-seq data for GRN estimation

The code for the single-cell experiments is located in:

- `single_cell_funcs.R`
- `single_cell_sim_mixed.R`

To reproduce these experiments, users must download the gene expression dataset from:

Chu et al. (2016), *Single-cell RNA-seq reveals novel regulators of human embryonic stem cell differentiation to definitive endoderm*, **Genome Biology**.

DOI: https://doi.org/10.1186/s13059-016-1033-x

