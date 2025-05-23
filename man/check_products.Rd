% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_products.R
\name{check_products}
\alias{check_products}
\title{Check and Subset Spatially Varying Cross Products}
\usage{
check_products(
  gene_pair_list,
  marginals,
  cov_mat,
  ncores,
  check_morani,
  p_thresh = 0.05
)
}
\arguments{
\item{gene_pair_list}{A data frame with two columns representing pairs of genes (by name or index).}

\item{marginals}{A matrix of standardized gene expression values (e.g., output from \code{fit_marginals()}).}

\item{cov_mat}{A data frame containing spatial coordinates (e.g., columns \code{x1} and \code{x2}).}

\item{ncores}{Number of CPU cores to use for parallel computation.}

\item{check_morani}{Logical indicating whether to filter based on Moran's I statistic.}

\item{p_thresh}{P-value threshold for filtering gene pairs by Moran's I. Default is 0.05.}
}
\value{
A list with:
\describe{
\item{gene_pair_list_subset}{Filtered gene pair list (if \code{check_morani = TRUE})}
\item{product_list}{List of gene-wise expression cross-products}
}
#'
}
\description{
Computes cross-products of standardized gene expression values and optionally filters
gene pairs based on spatial autocorrelation (Moran's I).
}
\examples{
\dontrun{
data(test_data)
# Fit standardized marginals for gene expressions
marginal_res <- fit_marginals(
  gene_list = test_data$gene_list,
  count_mat = test_data$count_mat,
  cov_mat = test_data$cov_mat,
  formula1 = "layer_annotations",
  family1 = "nb",
  DT = TRUE,
  epsilon = 1e-6,
  ncores = 2
)
# Check and subset spatially varying cross product
check_result <- check_products(
  gene_pair_list = test_data$gene_pair_list,
  marginals = marginal_res$marginal,
  cov_mat = test_data$cov_mat,
  check_morani = FALSE, ,
  ncores = 2
)
}

}
