#' SHAP values visualization
#'
#' This function takes the output from the `mlearn` function and creates either waterfall plot or force plot
#' showing the top n SHAP values for the specified number of features.
#'
#' @param shap_value SHAP values from `mlearn`
#' @param dbmart.concepts A table with feature description
#' @param top_n The number of top features to display.
#' @return
#' @export
#'
#' @examples
#' # Assuming `shap_value` is your data frame with SHAP values
#' mshapviz_all(shap_value, dbmart.concepts, top_n = 10)

position_bee <- function(width = NULL, adjust = NULL) {
  ggplot2::ggproto(NULL, PositionBee, width = width, adjust = adjust)
}

PositionBee <- ggplot2::ggproto(
  "PositionBee",
  ggplot2::Position,
  required_aes = c("x", "y"),
  
  setup_params = function(self, data) {
    list(
      width = if (!is.null(self$width)) self$width else
        ggplot2::resolution(data$y, zero = FALSE) * 0.4,
      adjust = if (!is.null(self$adjust)) self$adjust else 0.5
    )
  },
  
  compute_panel = function(self, data, params, scales) {
    data <- ggplot2::flip_data(data, params$flipped_aes)
    y_jit <- ave2(data$x, g = data$y, FUN = shifter, adjust = params$adjust)
    data <- ggplot2::transform_position(
      data, trans_y = function(y) y + y_jit * params$width
    )
    ggplot2::flip_data(data, params$flipped_aes)
  }
)

shifter <- function(y, ...) {
  if (length(y) == 1L) {
    return(0)
  }
  dens <- stats::density(y, ...)
  dens_y <- dens[["y"]] / max(dens[["y"]])
  shift <- halton_sequence(length(y))[rank(y, ties.method = "first")] - 0.5
  2 * shift * stats::approx(dens[["x"]], dens_y, xout = y)[["y"]]
}

ave2 <- function(x, g = NULL, FUN = mean, ...) {
  if (is.null(g)) {
    x[] <- FUN(x, ...)
  } else {
    split(x, g) <- lapply(split(x, g), FUN, ...)
  }
  x
}

halton_sequence <- function(n, b = 2) {
  vapply(seq_len(n), halton, FUN.VALUE = 0.0)
}

halton <- function(i, b = 2) {
  f <- 1
  r <- 0
  while (i > 0) {
    f <- f / b
    r <- r + f * (i %% b)
    i <- trunc(i / b)
  }
  r
}
.get_color_scale <- function(viridis_args, bar = TRUE, ncol = 2L) {
  if (bar) {
    viridis_args_plus <-
      list(
        breaks = if (ncol >= 2L) 0:1 else 0.5,
        labels = if (ncol >= 2L) c("Low", "High") else "Avg",
        guide = ggplot2::guide_colorbar(
          barwidth = 0.4,
          barheight = 8,
          title.theme = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0),
          title.position = "left"
        )
      )
  } else {
    viridis_args_plus <- list(guide = "none")
  }
  return(do.call(ggplot2::scale_color_viridis_c, c(viridis_args, viridis_args_plus)))
}

mshapviz_all <- function(shap_value, dbmart.concepts, top_n = 10) {
  top_feature_neg = shap_value %>%
    filter(sign == -1) %>%
    group_by(variable) %>%
    dplyr::summarise(mean_S = mean(contribution)) %>%
    arrange(desc(mean_S)) %>%
    top_n(top_n) %>%
    dplyr::mutate(
      from = cumsum(lag(mean_S, default = 0)),
      to = cumsum(mean_S),
      label = variable) %>% data.frame()
  
  top_feature_pos = shap_value %>%
    filter(sign == 1) %>%
    group_by(variable) %>%
    dplyr::summarise(mean_S = mean(contribution)) %>%
    arrange(desc(mean_S)) %>%
    top_n(top_n)  %>%
    dplyr::mutate(
      from = cumsum(lag(mean_S, default = 0)),
      to = cumsum(mean_S),
      label = variable) %>% data.frame()
  
  shap.feature <- rbind(top_feature_neg, top_feature_pos)
  shp <- shap_value[shap_value$variable %in% shap.feature$variable,]
  
  shap.dat = merge(shp, shap.feature, by = "variable")
  
  dbmart.concepts$phenx = as.character(dbmart.concepts$phenx)
  shap.dat <- merge(shap.dat, dbmart.concepts, by.x = "variable_name", by.y = "phenx")
  
  
  plot = ggplot2::ggplot(shap.dat, 
                         ggplot2::aes(x = contribution, y = variable_name)) +
    #ggplot2::geom_bar(ggplot2::aes(x =contribution/200, y = variable_name),stat = "identity")+
    ggplot2::geom_vline(xintercept = 0, color = "darkgray") +
    ggplot2::geom_point(ggplot2::aes(color = contribution),
                        position = position_bee(width = 0.4, adjust = 0.5)) +
    labs(color = "SHAP value") +
    ggplot2::theme_bw() +
    ggplot2::labs(y = ggplot2::element_blank(), x = "SHAP value") 
  
  print(plot)
}










