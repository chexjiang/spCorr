test_that("spCorr function runs without error", {
  data(test_data)
  result <- spCorr(
    count_mat = test_data$count_mat,
    gene_list = test_data$gene_list,
    gene_pair_list = test_data$gene_pair_list,
    cov_mat = test_data$cov_mat,
    formula1 = "layer_annotations",
    family1 = "nb",
    formula2 = "s(x1, x2, bs='tp', k=50)",
    family2 = quasiproductr(),
    DT = TRUE,
    global_test = "wald",
    ncores = 2,
    control = list(),
    seed = 123,
    preconstruct_smoother = TRUE,
    return_models = FALSE,
    return_coefs = FALSE,
    check_morani = FALSE
  )
  expect_type(result, "list")
})
