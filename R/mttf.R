#' @title Estimate mean time to failure based on 3-parameter Weibull
#'
#' @param x vector of beta, gamma, eta parameters of Weibull distribution
#'
#' @return mean and sd of time to failure from Weibull curve
#' @export
#'
mttf=function(x){
  b.in=x[1]
  g.in=x[2]
  n.in=x[3]

  mean.out=round(g.in+n.in*gamma((1/b.in)+1),4)
  sd.out=round(n.in*sqrt(gamma((2/b.in)+1)-gamma((1/b.in)+1)^2),4)
  return(paste(c(mean.out," (",sd.out,")"),collapse=""))
}
