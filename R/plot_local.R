#' Plot local correlation and gene expression for significant spots
#'
#' This function generates different types of visualizations for local correlation and gene expression
#' data, including correlation plots, expression plots, and local p-value significance plots.
#'
#' @param res_list A list containing results, including `test_res`, `gene_expr`, and `cov_mat`.
#' @param name A character string specifying the gene pair's name for which to plot the data.
#' @param area_name An optional character string indicating the annotation column name for subsetting (e.g., region name).
#' @param area_id An optional vector specifying the annotation values for filtering (`area_name` column). Should be categorical.
#' @param plot_type A character string specifying the type of plot to create. Options are `"rho"` (local correlations), `"all"` (both correlation and expression), and `"pval"` (local p-values). Default is `"all"`.
#' @param cov Optional data frame for covariates used in the plot. If not specified, `res_list$cov_mat` is used.
#' @return A ggplot object representing the selected type of visualization.
#' @examples
#' # Assuming `res_list` contains test results and gene expressions:
#' plot_rho <- plot_local(res_list, name = "geneA_geneB", plot_type = "rho")
#' plot_all <- plot_local(res_list, name = "geneA_geneB", area_name = "region", area_id = "region1", plot_type = "all")
#' @export
plot_local <- function(res_list,
                       name,      # gene pair's name
                       area_name=NULL, # annotation column name
                       area_id=NULL,   # annotation values (categorical)
                       plot_type='all',
                       cov=NULL){

  if(is.null(cov)){
    cov <- res_list$cov_mat
  }

  n <- nrow(res_list$gene_expr[[name]])
  res <- res_list$test_res[[name]]


  plot_id <- rep(TRUE, n)
  if(!is.null(area_name)){
    # subset area
    area_list <- res_list$gene_expr[[name]][,area_name]
    plot_id_false <- which(!area_list %in% area_id)
    plot_id[plot_id_false] <- FALSE
  }


  if(plot_type=='rho'){
    plot_dat <- data.frame(rho=res$local_test$fitted_rho, x1 = cov$x1, x2 = cov$x2, plot_id = plot_id)

    plot <- plot_dat %>%
      ggplot(aes(x = x1, y = -x2, color = rho)) +
      geom_point(size = 1.1) +
      scale_colour_gradientn(colors = viridis_pal(option = "magma")(12))  +
      labs(x = "x1", y = "x2", title = NULL, color = expression(hat(rho))) +
      theme_minimal() +
      coord_fixed(ratio = 1) +
      theme(legend.position = NULL,
            legend.direction = "vertical",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.2, "cm"),
            axis.text = element_blank() )


  }else if(plot_type=='all'){
    genes <- strsplit(name, "_")[[1]]
    gene_pair_expr <- res_list$gene_expr[[name]]
    plot_dat <- data.frame(gene_pair_expr, x1 = cov$x1, x2 = cov$x2,
                           rho=res$local_test$fitted_rho, plot_id = plot_id)

    # Plot rho
    plot1 <- plot_dat %>%
      ggplot(aes(x = x1, y = -x2, color = ifelse(plot_id == TRUE, rho, NA))) +
      geom_point(size = 0.5) +
      # Use two color scales: one for plot_id == TRUE and another for FALSE
      scale_colour_gradientn(colors = viridis_pal(option = 'magma')(10), na.value = "gray") +
      coord_fixed(ratio = 1) +
      labs(x = NULL, y = NULL, title = "Correlation", subtitle = name, color = expression(hat(rho))) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.65, "cm"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(1, "lines")
      )

    # Plot expression
    plot_long <- plot_dat %>%
      pivot_longer(cols = c(y1, y2), names_to = "variable", values_to = "value") %>%
      mutate(value_plot = ifelse(plot_id == TRUE, value, NA))  # Create a new column for conditional coloring

    facet_labels <- c("y1" = paste(genes[1]), "y2" = paste(genes[2]))

    plot2 <- plot_long %>%
      ggplot(aes(x = x1, y = -x2, color = log1p(value_plot))) +
      geom_point(size = 0.5) +
      scale_colour_gradientn(colors = viridis_pal(option = "plasma")(10), na.value = "gray") +
      coord_fixed(ratio = 1) +
      labs(x = NULL, y = NULL, title = "Expression", color = expression(log1p(y))) +
      facet_wrap(~variable, labeller = as_labeller(facet_labels)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.55, "cm"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(1, "lines")
      )

    plot <- plot2 + plot1 + plot_layout(ncol = 2, widths = c(1, 0.5))


  }else if(plot_type=='pval'){
    rho_p <- ifelse(res$local_test$rho_p>0.05, NA, res$local_test$rho_p)
    plot_dat <- data.frame(rho_p=-log1p(rho_p), x1 = cov$x1, x2 = cov$x2, plot_id = plot_id)

    plot <- plot_dat %>%
      ggplot(aes(x = x1, y = -x2, color = ifelse(plot_id == TRUE, rho_p, NA))) +
      geom_point(size = 0.8) +
      # Use two color scales: one for plot_id == TRUE and another for FALSE
      scale_colour_gradientn(colors = viridis_pal(option = 'magma')(10), na.value = "gray") +
      coord_fixed(ratio = 1) +
      labs(x = NULL, y = NULL, title = "Local p-value (significant)", subtitle = name, color = "-log1p(pval)") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.position = "bottom",
        axis.text = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.75, "cm"),
        panel.spacing = unit(1, "lines")
      )

  }

  plot
}
