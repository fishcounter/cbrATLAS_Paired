#' @title  Summarizes detection histories and provides Cormack-Jolly-Seber
#' 	estimates with or without adjustment of potential tag-failure
#'
#' @description CJS.fn summarizes the individual detection histories and
#' 	estimates seeds for Cormack-Jolly-Seber (CJS) likelhood function. It is based
#' 	on Skalski et al. (1998). This likelihood is also used to adjust survival for
#' 	estimated tag life as described in Townsend et al. (2006).  If no tag-life
#' 	data provided, only the unadjusted CJS survival will be returned.
#'
#' @param detect.in Detection history in wide ATLAS format without rel.group or bin columns
#' @param L.in Vector of mean probabilities tag is working, given detection at each
#'  site. If NULL, function will estimate the CJS parameters without tag failure
#' @param seeds.in Optional starting estimates for CJS model parameters Si, pi,li 
#' 	where i is 1:(total number of sites-1)
#' @param se.out (T|F) Should standard error of cjs parameters be estimated.
#' 	Uses variance matrix based on Hessian from MLE of cjs likelihood
#' @param d Proportion of detected tags censored at a location
#'
#' @return Returns a two-column matrix with CJS estimates and standard errors (if estimated)
#' 
#' @export
#'
cjs.fn=function(detect.in,L.in=NULL,seeds.in=NULL,se.out=F,d=NULL, parscale.in=1e-04, reltol.in=1e-30, maxit.in=2e+06){
  num.period=dim(detect.in)[2]
  if(is.null(L.in)){L.in=rep(1,num.period)} # if no tag failure, L = 100% for all sites

  # summary of counts.in/history and unique detection matrix in same order
  counts.in = thist0(detect.in)$count
  obs.hist=thist0(detect.in)$hist.matrix

  if(is.null(seeds.in)){
    # seed estimates do not estimate survival/detection with censoring present.
    seeds.in = rep(0.5, num.period * 2 - 1)
    names(seeds.in)=c(paste0(rep(c("S","p"),rep(num.period-1,2)),rep(1:(num.period-1),2)),"l")

    if (num.period > 2) {
      seeds.in[num.period] = sum(detect.in[apply(detect.in[,-1]==1, 1, sum) > 0, 1] == 1)/sum(apply(detect.in[, -1] == 1, 1, sum) > 0)
      for (i in (2:(num.period - 2))) {
        seeds.in[num.period + i - 1] = sum(detect.in[apply(detect.in[,-c(1:i)] == 1, 1, sum) > 0, i] == 1)/
          sum(apply(detect.in[,-c(1:i)] == 1, 1, sum) > 0)
      }
      seeds.in[num.period * 2 - 2] = sum(detect.in[detect.in[,num.period] == 1, num.period - 1] == 1)/sum(detect.in[,num.period] == 1)
    }else{seeds.in[num.period] = sum(detect.in[detect.in[, 2] == 1,1] == 1)/sum(detect.in[, 2] == 1)}

    seeds.in[(2 * num.period) - 1] = sum(detect.in[, num.period] == 1)/(sum(detect.in[, num.period - 1] == 1)/seeds.in[2 * (num.period - 1)])
    seeds.in[1] = (sum(detect.in[, 1] == 1)/dim(detect.in)[1])/seeds.in[num.period]
    if (num.period > 2) {
      for (i in 2:(num.period - 1)) {
        seeds.in[i] = (sum(detect.in[, i] == 1)/seeds.in[i + num.period - 1])/(sum(detect.in[, i - 1] == 1)/seeds.in[i + num.period - 2])
      }
    }
    # "p(censored)")
    d=rep(0,num.period-1)
    for(i in (1:(num.period-1))){
      d[i]= sum(detect.in[,i]==2)/sum(detect.in[,i]>0)
      names(d)[i]=paste0("d",i)
    }
  }

  # adjust S,lambda by expected tag-life.  L.in
  seeds.in[c(1:(num.period - 1), 2 * num.period - 1)] = seeds.in[c(1:(num.period - 1), 2 * num.period - 1)]/L.in
  # print("seeds done ")

  ## maximum likelihood estimate of CJS parameters

  cjs.out = stats::optim(seeds.in, cjs.lik, counts.in = counts.in, num.period = num.period, use.hist=obs.hist,L=L.in,d.in=d,
                  method = "BFGS", hessian = F, control = list(parscale = rep(parscale.in, num.period * 2 - 1), reltol = reltol.in, maxit=maxit.in,
				  ndeps=rep(parscale.in, num.period * 2 - 1)))
  cjs.out$par=correct.fn(cjs.out$par)


  #print(cjs.out$par)
  if(se.out){
    # rerun to get hessian
    cjs.out = stats::optim(cjs.out$par, cjs.lik, counts.in = counts.in, num.period = num.period, use.hist=obs.hist,L=L.in,d.in=d,
                    method = "BFGS", hessian = T, control = list(parscale = rep(parscale.in, num.period * 2 - 1), reltol = reltol.in,
					maxit=maxit.in,ndeps=rep(parscale.in, num.period * 2 - 1)))
    cjs.out$par=correct.fn(cjs.out$par)

    # calculate standard errors on CJS estimates, based on Hessian
    temp = cjs.out$hessian
    singular = colSums(temp) == 0
    var.vec = rep(0, num.period * 2 - 1)
    var.vec[!singular] = diag(solve(temp[!singular, !singular]))

    # format estimates into a matrix, 2 columns. 1 row per parameter
    cjs.param = matrix(c(cjs.out$par, sqrt(var.vec)), ncol = 2)
    row.names(cjs.param)=names(seeds.in)
    dimnames(cjs.param)[[2]] = c("Estimate", "s.e.")
  }else{
    # format estimates into a matrix, 1 columns. 1 row per parameter
    cjs.param = matrix(c(cjs.out$par), ncol = 1)
    row.names(cjs.param)=names(seeds.in)
  }
	if(is.null(d)){return(cjs.param)}else{
		return(list(cjs.param=cjs.param, d=d))
		}
}		
