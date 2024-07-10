#' SHAP values visualization
#'
#' This function takes the output from the `mlearn` function and creates either waterfall plot or force plot
#' showing the top n SHAP values for the specified number of features.
#'
#' @param data.shap Model output containing SHAP values from `mlearn`
#' @param data.descp A table with feature description
#' @param top_n The number of top features to display.
#' @param plot_type 'waterfall' or 'force'
#' @param num The nth of incidence picked
#' @return
#' @export
#'
#' @examples
#' # Assuming `shap_value` is your data frame with SHAP values
#' visualize_shap(shap_value, dbmart.concepts, plot_type = "waterfall", top_n = 5, num = 6)
mshapviz <- function(shap_value, dbmart.concepts, plot_type = "force", top_n = 5, num = 1) {
  shp1 <- shap_value %>%
    filter(B == num) %>%
    mutate(S = as.numeric(contribution),
           abs_S = abs(S))
  
  shp1.neg = shp1 %>%
    filter(sign == -1) %>%
    arrange(desc(abs_S)) %>%
    top_n(top_n) %>%
    dplyr::mutate(
      from = cumsum(lag(contribution, default = 0)),
      to = cumsum(contribution),
      label = variable) %>% data.frame()
  
  shp1.pos = shp1 %>%
    filter(sign == 1) %>%
    arrange(desc(abs_S)) %>%
    top_n(top_n) %>% 
    dplyr::mutate(
      from = cumsum(lag(contribution, default = 0)),
      to = cumsum(contribution),
      label = variable) %>% data.frame()
  
  shap.dat <- rbind(shp1.pos, shp1.neg)
  
  if (plot_type == "force") {
    shap.dat$id = "1"
  } else if (plot_type == "waterfall") {
    shap.dat <- shap.dat %>%
      arrange(desc(abs_S)) %>%
      dplyr::mutate(id = row_number())
  }
  
  shap.dat$contribution = round(shap.dat$contribution, 2)
  dbmart.concepts$phenx = as.character(dbmart.concepts$phenx)
  shap.dat <- merge(shap.dat, dbmart.concepts, by.x = "variable_name", by.y = "phenx")
  
  fill_colors <- c('TRUE' = '#f7d13d', 'FALSE' = '#a52c60')  
  
  plot = ggplot2::ggplot(
    shap.dat,
    ggplot2::aes(
      xmin = from, xmax = to, y = id, label = contribution,forward = F,
      fill = factor(S < 0, levels = c(FALSE, TRUE))
    )
  ) +
    gggenes::geom_gene_arrow(
      show.legend = FALSE,
      arrowhead_width = grid::unit(2, "mm")
    ) + gggenes::geom_gene_label() +
    ggrepel::geom_text_repel(
      ggplot2::aes(x = (from + to) / 2, y = as.numeric(id) + 0.08, label = DESCRIPTION),
      nudge_y = 0.3,
      segment.size = 0.1,
      segment.alpha = 0.5,
      direction = "both"
    ) +
    ggplot2::coord_cartesian() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.13)) +
    ggplot2::scale_fill_manual(values = fill_colors, drop = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1 / 4,
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(y = ggplot2::element_blank(), x = "SHAP value") 
  
  print(plot)
}










