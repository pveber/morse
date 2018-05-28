#' Predict the Multiplication Factor leading to x\% of reduction in survival
#' at a specific time.
#' 
#' Generic method for \code{MFx}, a function denoted \eqn{MF(x,t)} for 
#' \eqn{x}\% Multiplication Factor at time \eqn{t}.
#' 
#' The function \code{MFx}, \eqn{x}\% Multiplication Factor at time \eqn{t}, (\eqn{MF(x,t)}),
#' is used to compute the multiplication factor
#' applied to the concentration exposure profile in order to
#' reduce by \eqn{x}\% (argument \code{X}) the survival rate at a
#'  specified test duration \eqn{t} (argument \code{time_MFx}) (default is the maximum
#'  time point of the experiment).
#'  
#'  Mathematical definition of \eqn{x}\% Multiplication Factor at time \eqn{t}
#'  (at the end of a time series \eqn{T = \{0, \dots, t\}}),
#'  denoted \eqn{MF(x,t)}, is given by:
#'  
#'  \eqn{S(MF(x,t) * C_w(\tau \in T), t) = S( C_w(\tau \in T), t)*(1- x/100)},
#'  
#'  where \eqn{C_w(\tau \in T)} is the initial exposure profile without
#'  multiplication factor. And so the expression \eqn{S(MF(x,t)* C_w(\tau \in T), t)}
#'  is the survival rate after an exposure profile
#'  \eqn{MF(x,t)* C_w(\tau \in T)} at time \eqn{t}.
#' 
#' @param object An object used to select a method
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#' 
MFx <- function(object, ...){
  UseMethod("MFx")
}

