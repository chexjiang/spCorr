#' Check and Subset Spatially Varying Cross Products
#'
#' Computes cross-products of standardized gene expression values and optionally filters
#' gene pairs based on spatial autocorrelation (Moran's I).
#'
#' @param gene_pair_list A data frame with two columns representing pairs of genes (by name or index).
#' @param marginals A matrix of standardized gene expression values (e.g., output from `fit_marginals()`).
#' @param cov_mat A data frame containing spatial coordinates (e.g., columns `x1` and `x2`).
#' @param ncores Number of CPU cores to use for parallel computation.
#' @param check_morani Logical indicating whether to filter based on Moran's I statistic.
#' @param p_thresh P-value threshold for filtering gene pairs by Moran's I. Default is 0.05.
#'
#' @return A list with:
#' \describe{
#'   \item{gene_pair_list_subset}{Filtered gene pair list (if `check_morani = TRUE`)}
#'   \item{product_list}{List of gene-wise expression cross-products}
#' }
#' #'
#' @examples
#' \dontrun{
#' data(test_data)
#' # Fit standardized marginals for gene expressions
#' marginal_res <- fit_marginals(
#'   gene_list = test_data$gene_list,
#'   count_mat = test_data$count_mat,
#'   cov_mat = test_data$cov_mat,
#'   formula1 = "layer_annotations",
#'   family1 = "nb",
#'   DT = TRUE,
#'   epsilon = 1e-6,
#'   ncores = 2
#' )
#' # Check and subset spatially varying cross product
#' check_result <- check_products(
#'   gene_pair_list = test_data$gene_pair_list,
#'   marginals = marginal_res$marginal,
#'   cov_mat = test_data$cov_mat,
#'   check_morani = FALSE, ,
#'   ncores = 2
#' )
#' }
#'
#' @importFrom parallel mclapply
#' @importFrom stats dist
#' @importFrom ape Moran.I
#' @export
check_products <- function(gene_pair_list,
                           marginals,
                           cov_mat,
                           ncores,
                           check_morani,
                           p_thresh = 0.05) {
  # Compute cross product
  product_list <- return_product_list(
    gene_pair_list = gene_pair_list,
    marginals = marginals,
    cov_mat = cov_mat,
    ncores = ncores
  )

  if (check_morani) {
    # Calculate morani for each product
    weight_mat <- as.matrix(dist(cov_mat[, c("x1", "x2")]))
    weight_mat <- 1 / weight_mat
    diag(weight_mat) <- 0

    morani_list <- fit_products_morani(
      product_list = product_list,
      weight_mat = weight_mat,
      ncores = ncores
    )

    morani_pval_list <- lapply(morani_list, function(morani) morani$p.value)
    subset_name <- names(which(morani_pval_list < p_thresh))
    gene_pair_list_subset <- gene_pair_list[subset_name, ]
  } else {
    gene_pair_list_subset <- gene_pair_list
  }

  return(list(
    gene_pair_list_subset = gene_pair_list_subset,
    product_list = product_list
  ))
}


return_product_list <- function(gene_pair_list,
                                marginals,
                                cov_mat,
                                ncores) {
  # Split gene pair list into individual rows
  gene_pair_list_name <- row.names(gene_pair_list)
  gene_pair_list_split <- split(gene_pair_list, seq(nrow(gene_pair_list)))
  names(gene_pair_list_split) <- gene_pair_list_name
  # Fit models for each gene pair in parallel
  product_list <- parallel::mclapply(gene_pair_list_split, function(gene_pair) {
    # Extract standardized expressions for the gene pair
    gene_pair <- unlist(gene_pair)
    g1 <- marginals[gene_pair[1], ]
    g2 <- marginals[gene_pair[2], ]

    # Compute product of expressions
    z <- g1 * g2
  }, mc.cores = ncores)

  return(product_list)
}


fit_products_morani <- function(product_list,
                                weight_mat,
                                ncores) {
  # Fit models for each gene pair in parallel
  morani_list <- parallel::mclapply(product_list, function(z) {
    # Run Moran's I
    morani_list <- ape::Moran.I(x = z, weight = weight_mat, scaled = FALSE, na.rm = TRUE, alternative = "g")
    morani_list
  }, mc.cores = ncores)

  return(morani_list)
}
