#' Plotting method for \code{gm_survData} objects
#'
#' This is the generic \code{plot} S3 method for the \code{gm_survData} class.
#' It plots the number of survivors as a function of time.
#'
#' @param x an object of class \code{gm_survData}
#' @param xlab a title for the \eqn{x}-axis (optional)
#' @param ylab a label for the \eqn{y}-axis
#' @param main main title for the plot
#' @param facetting a logical value to use facetting on profile
#' @param level_facetting a vector of character...
#'
#' @note the function calls function \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#'
#' @keywords plot
#'
#' @examples
#'
#'
#' @import ggplot2
#' @importFrom graphics plot axis legend lines par points title
#' @importFrom methods is
#' @importFrom stats aggregate
#'
#' @export

plot.gm_survData <- function(data,
                             lab.x = "Time",
                             lab.y = "Number of surviving individuals",
                             lab.main = NULL,
                             facetting = TRUE,
                             facetting.level = NULL){


  data_plt = filter( data , !is.na(Nsurv))

  if(!is.null(facetting.level)){
    data_plt$profile = factor(data_plt$profile, levels = facetting.level)
  }


  plt = ggplot(data_plt, aes(x = time, y = Nsurv, group = profile, color = conc)) +
          theme_minimal() +
          theme(
            ## facet elements
            # strip.background = element_rect(fill="orange", colour = "orange"),
            # strip.text = element_text(colour = "grey30"),
            ## legend
            legend.title = element_text(size = 9),
            legend.text=element_text(size = 7),
            legend.key.width = unit(0.3, "in"),
            legend.key.height = unit(0.15, "in"),
            legend.position = "top"
          ) +
          labs(title = lab.main,
               x = lab.x,
               y = lab.y,
               colour = "Concentration" # legend title
          ) +
          scale_colour_gradient(
            name="Concentration",
            low="grey20", high="orange"
          ) +
          scale_fill_gradient(
            name="Concentration",
            low="grey20", high="orange"
          ) +
          expand_limits(x = 0, y = 0) +
          geom_point() +
          geom_line()

  if(facetting == TRUE){
    plt = plt + facet_wrap(~ profile, ncol = 2)
  } else if(facetting != TRUE){
    plt
  } else stop("'facetting' is a logical 'TRUE' or 'FALSE'.")

  return(plt)
}

