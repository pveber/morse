#' Plotting method for \code{survDataVarExp} objects
#'
#' This is the generic \code{plot} S3 method for the \code{survDataVarC} class.
#' It plots the number of survivors as a function of time.
#'
#' @param x an object of class \code{survDataVarExp}
#' @param xlab a label for the \eqn{X}-axis, by default \code{Time}
#' @param ylab a label for the \eqn{Y}-axis, by default \code{Number of survivors}
#' @param main main title for the plot
#' @param one.plot if \code{TRUE}, draws all the points in one plot instead of one per \code{replicate}
#' @param facetting_level a vector of \code{characters} to rank \code{replicates} in the multi plot (i.e. \code{one.plot == FALSE})
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return an object of class \code{ggplot}, see function \code{\link[ggplot2]{ggplot}}
#'
#' @keywords plot
#'
#' @import ggplot2
#' @importFrom graphics plot axis legend lines par points title
#' @importFrom methods is
#' @importFrom stats aggregate
#'
#' @export


plot.survDataVarExp <- function(x,
                                xlab = "Time",
                                ylab = "Number of survivors",
                                main = NULL,
                                one.plot = FALSE,
                                facetting_level = NULL,
                                ...){
  
  data_plt = filter( x , !is.na(Nsurv))
  
  if(!is.null(facetting_level)){
    data_plt$replicate = factor(data_plt$replicate, levels = facetting_level)
  }
  
  plt = ggplot(data_plt, aes(x = time, y = Nsurv, group = replicate, color = conc)) +
    theme_minimal() +
    theme(
      ## facet elements
      strip.background = element_rect(fill="white", colour = "white"),
      strip.text = element_text(colour = "black"),
      ## legend
      legend.title = element_text(size = 9),
      legend.text=element_text(size = 7),
      legend.key.width = unit(0.3, "in"),
      legend.key.height = unit(0.15, "in"),
      legend.position = "top"
    ) +
    labs(title = main,
         x = xlab,
         y = ylab,
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
  
  if(one.plot == FALSE){
    plt <- plt + facet_wrap(~ replicate, ncol = 2)
  } else if(one.plot == TRUE){
    plt
  } else stop("'one.plot' is a logical 'TRUE' or 'FALSE'.")
  
  return(plt)
}

