#' @title Estimate harmonic mean
#'
#' @param t.in numeric vector (no NA's)
#'
#' @return character vector of harmonic mean and se
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
