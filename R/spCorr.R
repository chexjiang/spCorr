library(mgcv)
library(MASS)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(parallel)
library(ggplot2)
library(patchwork)



source('mgcv_modified.R')
source('construct_smoother.R')
source('quasiproductr.R')
source('fit_marginals.R')
source('fit_products.R')
source('test_models.R')
source('summary.R')
source('plot_helper.R')

spCorr <- function(count_mat,
                   gene_list,
                   gene_pair_list,
                   cov_mat,
                   formula1 = "1",
                   family1 = 'nb',
                   formula2 = "s(x1, x2, bs='tp', k=50)",
                   family2 = quasiproductr(),
                   DT = TRUE,
                   return_models = FALSE,
                   ncores = 4,
                   control = list(),
                   seed = 123,
                   local_testing = FALSE,
                   preconstruct_smoother = TRUE){

  set.seed(seed)

  # Fit marginal to gene_list
  message("Start Marginal Fitting for ", length(gene_list), " genes")
  marginals <- covert_cond_marginals(gene_list = gene_list,
                                     count_mat = count_mat,
                                     cov_mat = cov_mat,
                                     formula1 = formula1,
                                     family = family1,
                                     to = 'gaussian',
                                     DT = DT,
                                     ncores = ncores)


  # Fit product to gene_pair_list
  message("Start Product Fitting for ", nrow(gene_pair_list), " gene pairs")
  smoother_env <- new.env()
  model_list <- fit_products(gene_pair_list = gene_pair_list,
                             marginals = marginals,
                             cov_mat = cov_mat,
                             formula2 = formula2,
                             family2 = quasiproductr(),
                             control = control,
                             ncores = ncores,
                             preconstruct_smoother = preconstruct_smoother)

  # Extract the gene expression list
  gene_pair_expr_list <- lapply(1:nrow(gene_pair_list), function(i) {
    # Extract the expressions for the gene pair
    y1 <- count_mat[gene_pair_list[i, 1], ]
    y2 <- count_mat[gene_pair_list[i, 2], ]

    # Combine with cov_mat
    gene_pair_expr <- cbind(y1=y1, y2=y2)
    gene_pair_expr
  })
  names(gene_pair_expr_list) <- row.names(gene_pair_list)


  # Testing for gene_pair_list
  message("Start Testing for ", nrow(gene_pair_list), " gene pairs")
  test_res <- test_models(model_list, ncores, local_testing)

  if(return_models){
    return(list(test_res=test_res,
                gene_expr=gene_pair_expr_list,
                cov_mat=cov_mat,
                model_list=model_list))
  }else{
    return(list(test_res=test_res,
                gene_expr=gene_pair_expr_list,
                cov_mat=cov_mat))
  }
}

