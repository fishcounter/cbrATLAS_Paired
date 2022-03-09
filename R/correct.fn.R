#' @title Restrains uncorrected probability estimates between 0 and 1
#'
#' @description Estimation of Cormack-Jolly-Seber (CJS) survival parameters can
#'   be greater than 1 when detection probabilities are low. To reduce computational
#'   failure when maximizing the likelihood, correct.fn bounds the CJS parameter
#'   estimates between 0 and 1 in the likelihood.  This may result in survival
#'   estimates biased low, but low detection rates do not usually occur in studies
#'   with active tag technology.
#'
#' @param x Value to be rounded up/down to be between 0 to 1
#'
#' @return Returns original value if between 0 and 1, otherwise rounds up to 1e-10 or down to 1-1e-10
#' @export
#'
correct.fn=function(x){
  # keep probabilities between 0 and 1
  x[x<0.0000001]=1e-10
  x[x>0.9999999]=1-1e-10
  return(x)
}

