library(devtools)
devtools::build()
devtools::install()


library(spCorr)
data(test_data)

time1 <- Sys.time()
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
                     ncores = 1,
                     control = list(),
                     seed = 123,
                     local_testing = FALSE,
                     preconstruct_smoother = TRUE)
time2 <- Sys.time()
time2 - time1


time1 <- Sys.time()
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
                     ncores = 1,
                     control = list(),
                     seed = 123,
                     local_testing = FALSE,
                     preconstruct_smoother = FALSE)
time2 <- Sys.time()
time2 - time1






