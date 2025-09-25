#' @title Estimates a ratio of 2 survival estimates w/standard errors
#'
#' @description This function uses the point and stardard estimates to get the adjust survival and standard error. default is no covariance between survivals
#'
#' @param s1 survival upper reach
#' @param s2 survival lower reach
#' @param se1 survival upper reach standard error
#' @param se2 survival lowerr reach standard error
#' @param cov.s1s2=0 covariance between s1 & s2, default = 0
#' @param se.in=T  are standard errors or variances being used for inputs, default = SE
#' @param se.output=T  is output of survival ratio point estimate with standard error or variance, default = SE
#'
#' @return This function returns a vector with:
#' \describe{
#'    \item{point estimate with standard error or variance}{s1 = upper release survival, s2=lower release survival, se's are corresponding standard errors}
#' }
#'
#' @export
#'
vs2.fn=function(s1, s2, se1, se2,se.in=T,se.output=T){
## gives ratio of two estimates with s.e.s as a number

	if(se.in) {
	        v1 <- se1^2
			v2 <- se2^2
			}else{
			v1=se1
			v2=se2
			}
	s.out <- s1/s2
	v.out <- (s.out^2) * ((v1/(s1^2)) +(v2/(s2^2)) + (v1*v2)/((s1*s2)^2))
	se.out <- round(sqrt(v.out), 6)
	v.out=round(v.out,6)
	s.out <- round(s.out, 6)
	ifelse(se.output,return(cbind(s.out, se.out)),return(cbind(s.out, v.out)))
}

