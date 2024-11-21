#' Fit product distributions to standardized gene expressions for a list of gene pairs.
#'
#' This function fits product distributions to gene pairs using their standardized expression values.
#' It applies a Generalized Additive Model (GAM) to each pair of genes from `gene_pair_list` to model
#' the joint effect while incorporating covariates.
#'
#' @param gene_pair_list A data frame or matrix where each row represents a gene pair to be modeled.
#' @param marginals A matrix containing standardized expression values for each gene.
#' @param cov_mat A matrix or data frame of covariates used in fitting the GAM models.
#' @param formula2 A string or formula specifying the GAM structure for modeling the product distribution.
#' @param family2 The distribution family used in fitting the model. Default is `quasiproductr()`.
#' @param control A list of control parameters passed to the GAM fitting process.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param preconstruct_smoother Logical; if `TRUE`, modifies the smoother to enable caching for faster computation. Default is `TRUE`.
#' @return A list of fitted models, where each model corresponds to a gene pair from `gene_pair_list`.
#' @examples
#' data(test_data)
#' # Fit standardized marginals for gene expressions
#' marginals <- fit_marginals(
#'   gene_list = test_data$gene_list,
#'   count_mat = test_data$count_mat,
#'   cov_mat = test_data$cov_mat,
#'   formula1 = "layer_annotations",
#'   family1 = "nb",
#'   to = "gaussian",
#'   DT = TRUE,
#'   ncores = 2
#' )
#'
#' # Fit product distributions for gene pairs
#' model_list <- fit_products(
#'   gene_pair_list = test_data$gene_pair_list,
#'   marginals = marginals,
#'   cov_mat = test_data$cov_mat,
#'   formula2 = "s(x1, x2, bs='tp', k=50)",
#'   family2 = quasiproductr(),
#'   control = list(),
#'   ncores = 2,
#'   preconstruct_smoother = TRUE
#' )
#'
#' @importFrom parallel mclapply
#' @export

fit_products <- function(gene_pair_list,
                         marginals,
                         cov_mat,
                         formula2,
                         family2 = quasiproductr(),
                         control = list(),
                         ncores = ncores,
                         preconstruct_smoother = TRUE) {

  gene_pair_list_name <- row.names(gene_pair_list)
  gene_pair_list_split <- split(gene_pair_list, seq(nrow(gene_pair_list)))
  names(gene_pair_list_split) <- gene_pair_list_name


  model_list <- parallel::mclapply(gene_pair_list_split, function(gene_pair) {
    fit_product(gene_pair = gene_pair,
                marginals = marginals,
                cov_mat = cov_mat,
                formula2 = formula2,
                family2 = family2,
                control = control,
                preconstruct_smoother = preconstruct_smoother)
  }, mc.cores = ncores)

  return(model_list)
}


fit_product <- function(gene_pair,
                        marginals,
                        cov_mat,
                        formula2,
                        family2 = quasiproductr(),
                        control = list(),
                        preconstruct_smoother = TRUE){


  # Extract the standardized expressions
  gene_pair <- unlist(gene_pair)
  gj <- marginals[gene_pair[1],]
  gk <- marginals[gene_pair[2],]

  # Model the product distribution
  z <- gj*gk
  dat <- cbind(z, cov_mat)

  # Create the formula
  #mgcv_formula <- stats::reformulate(formula2, response = "z")

  if(preconstruct_smoother){
    # Original formula with the response variable
    bs_type <- sub(".*bs='(.*?)'.*", "\\1", formula2)

    # Modify the formula based on bs_type
    if (bs_type == "tp") {
      formula2_cached <- gsub("bs='tp'", "bs='tpcached'", formula2)
    } else if (bs_type == "gp") {
      formula2_cached <- gsub("bs='gp'", "bs='gpcached'", formula2)
    } else {
      stop("bs must be either 'tp' or 'gp'.")
    }
    mgcv_formula <- as.formula(paste("z ~", formula2_cached))
  }else{
    mgcv_formula <- as.formula(paste("z ~", formula2))
  }


  # Fit the model
  model <- suppressWarnings({
    mgcv::gam(formula = mgcv_formula, family = family2,
              data = dat, control = control)
  })

  model
}
