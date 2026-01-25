remove.packages("spCorr")
unlink("man", recursive = TRUE)
unlink("NAMESPACE")
devtools::document()
styler::style_pkg()
devtools::document()
devtools::build()
devtools::install()
devtools::check()
pkgdown::build_site()

library(spCorr)
data(test_data)

time1 <- Sys.time()
res1 <- spCorr(count_mat = test_data$count_mat,
               gene_list = test_data$gene_list,
               gene_pair_list = test_data$gene_pair_list,
               cov_mat = test_data$cov_mat,
               formula1 = "layer_annotations",
               family1 = 'nb',
               formula2 = "s(x1, x2, bs='tp', k=50)",
               family2 = quasiproductr(),
               DT = TRUE,
               global_test = "lrt",
               ncores = 2,
               control = list(),
               seed = 123,
               preconstruct_smoother = TRUE,
               return_models = FALSE,
               return_coefs = FALSE,
               check_morani = FALSE)
time2 <- Sys.time()
time2 - time1


time1 <- Sys.time()
res2 <- spCorr(count_mat = test_data$count_mat,
               gene_list = test_data$gene_list,
               gene_pair_list = test_data$gene_pair_list,
               cov_mat = test_data$cov_mat,
               formula1 = "layer_annotations",
               family1 = 'nb',
               formula2 = "s(x1, x2, bs='tp', k=50)",
               family2 = quasiproductr(),
               DT = TRUE,
               global_test = "LRT",
               ncores = 1,
               control = list(),
               seed = 123,
               preconstruct_smoother = FALSE,
               return_models = FALSE,
               return_coefs = FALSE,
               check_morani = FALSE)
time2 <- Sys.time()
time2 - time1

