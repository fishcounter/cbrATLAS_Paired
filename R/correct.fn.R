#' @title Restrains uncorrected probability estimates between 0 and 1.
#'
#' @description The Cormack-Jolly-Seber (CJS) likelihood estimates are currently bounded between 0 and
#' 1, which may bias survival estimates low in cases of low detection probability.
#' This is not usually the case with active tag technology with typically high
#' detection rates. The CJS estimate may estimated greater than 100%, and this may
#' be included in later iterations.
#'
#' @param x Value to be rounded up/down to be between 0 to 1.
#'
#' @return Returns original value if between 0 and 1, otherwise rounds up to 1e-10 or down to 1-1e-10.
#' @export
#'
correct.fn=function(x){
  # keep probabilities between 0 and 1
  x[x<0.0000001]=1e-10
  x[x>0.9999999]=1-1e-10
  return(x)
}

