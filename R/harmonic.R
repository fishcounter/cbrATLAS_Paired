#' @title Estimates harmonic mean
#'
#' @description This function estimates harmonic mean.
#'
#' @param t.in Numeric vector (no NAs).
#'
#' @return This function returns a character vector of harmonic mean and standard error
#'
#' @export
#'
harmonic=function(t.in){
   w=1/t.in
   t.out=length(w)/sum(w)
   t.var=(stats::var(w)/length(w))*(t.out^4)
   t.se=sqrt(t.var)
   return(paste(c(sprintf("%.4f",t.out),"  (",sprintf("%.4f",t.se),")      ",length(t.in)),collapse=""))
}
