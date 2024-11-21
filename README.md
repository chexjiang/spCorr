# spCorr

**spCorr** models spatially variable gene co-expression patterns in spatial transcriptomics.

---

## Installation<a name="installation-"></a>

To install the development version from GitHub, please run:

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("chexjiang/spCorr")
```

## Quick Start<a name="quick-start"></a>

The following code is a quick example of running spCorr. The function `spCorr()` takes in `count_mat`, `gene_list`, `gene_pair_list`, and `cov_mat` (which should include spatial coordinates and any confounders).


``` r
model_list <- spCorr(count_mat = test_data$count_mat,
                     gene_list = test_data$gene_list,
                     gene_pair_list = test_data$gene_pair_list,
                     cov_mat = test_data$cov_mat,
                     formula1 = "layer_annotations",
                     family1 = 'nb',
                     formula2 = "s(x1, x2, bs='tp', k=50)",
                     family2 = quasiproductr(),
                     DT = TRUE,
                     return_models = FALSE,
                     ncores = 2,
                     control = list(),
                     seed = 123,
                     local_testing = FALSE,
                     preconstruct_smoother = TRUE)
```

The parameters of `spCorr()` are:

- `count_mat`: A matrix of counts where rows represent genes and columns represent observations.

- `gene_list`: A list of gene names for which the conditional margins are to be fit.

- `gene_pair_list`: A data frame or matrix containing pairs of gene names (or indices) to be analyzed.

- `cov_mat`: A matrix of covariates, including spatial coordinates and any confounders.

- `formula1`: A formula specifying the model to be used for fitting the marginal distributions.  

- `family1`: The distribution family for marginals (e.g., `'gaussian'`, `'poisson'`, `'nb'`).  

- `formula2`: A formula specifying the model for fitting the product distributions.  
    - *Examples*: `"s(x1, x2, bs='tp', k=50)"`, `"s(x1, x2, bs='gp', k=50)"`.

- `family2`:  The distribution family for the product model.  Default is `quasiproductr()`.

- `DT`: Boolean indicating whether to apply a discrete transformation during marginal fitting. Default is `TRUE`.

- `return_models`: Boolean indicating whether to return the model objects along with results. Default is `FALSE`.

- `ncores`: The number of cores for parallel processing. Default is 2. 

- `control`: A list of control parameters for the fitting functions.

- `seed`: A seed for reproducibility. Default is `123`.

- `local_testing`: Boolean indicating whether to perform local significance testing for each gene pair. Default is `FALSE`.

- `preconstruct_smoother`: Boolean indicating whether to use a cached smoother for faster computation. Default is `TRUE`.

