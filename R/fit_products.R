## Fit product distributions to standardized expressions for a gene pair list
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


  model_list <- mclapply(gene_pair_list_split, function(gene_pair) {
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


## Fit product distribution to standardized expressions for a gene pair (j,k)
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
