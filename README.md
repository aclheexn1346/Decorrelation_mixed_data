## Mixed Data De-correlation Methodology

This repository implements a **de-correlation methodology for mixed data (continuous and discrete variables)**. The goal of the procedure is to transform correlated observations into approximately **independent latent variables**, which can then be used in downstream tasks such as causal discovery.

The main implementation is contained in: main.R


---

## Overview

The method estimates latent variables and a covariance structure for mixed data using an **iterative EM-style procedure**. It handles both continuous and discrete variables by introducing latent Gaussian variables for the discrete observations.

The main function:

```r
mixed_data_decor(data, block_sizes)
```

A user need only input the dataset and block sizes among the units of the data. For example, if you have a dataset with $n = 20$ observations and every 5 units are dependent,
```r
block_sizes = c(5,5,5,5)
```

---

## Required Files

Ensure that the packages are installed within libraries.R

```r
source("libraries.R")
source("helperFuncMixed.R")
```

helperFuncMixed.R contains the individual functions used in main.R

---

## Example Usage

Download test_data.RData, which contains a dataset with $n = 100$ observations and $p = 100$ variables. The example block structure suggests every 5 units are dependent. (Obs. 1 - 5, Obs. 6 - 10, etc. are dependent)

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

