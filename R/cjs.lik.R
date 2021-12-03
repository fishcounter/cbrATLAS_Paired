#' @title Likelihood function for unadjusted- or adjusted-for-taglife CJS
#'
#' @param params the likelihood MLE estimators (S,p,l) in a single vector
#' @param counts.in vector of summed detections by history (names are history)
#' @param num.period number of detection periods in history
#' @param use.hist matrix of unique detection histories
#' @param L.in The probability that a tag is working at detection site time. Set to 1 if active tags not indicated. p(tag working) is cumulative to each site, not multiplicate from site-to-site
#' @param d.in the probability a tag is censored, given that it's detected
#'
#' @return negative log-likelihood
#' @export
#'
cjs.lik=function(params,counts.in,num.period,use.hist,L.in=NULL,d.in=NULL){
  if(is.null(L.in)){L.in=rep(1,num.period)}

  if(!(num.period>1)){cat("Number of periods must be greater than 1")}
  stopifnot(num.period>1)

  params=correct.fn(params)

  # split mle vector into separate types
  s=params[1:(num.period-1)]
  p=params[(num.period):((num.period-1)*2)]
  l=params[1+(num.period-1)*2]
  if(is.null(d.in)){d=rep(0,num.period-1)}else{d=d.in}

  # combinatoric
  out=lgamma(sum(counts.in)+1)-sum(lgamma(counts.in+1))

  # fish detected in last period
  for(i in c(1:dim(use.hist)[1])[use.hist[,dim(use.hist)[2]]==1]){
    out.hist=log(l)+log(L.in[num.period])
    for(j in 1:(num.period-1)){
      if(use.hist[i,j]==1){out.hist=out.hist+log(s[j])+log(p[j])+log(1-d[j])}
      if(use.hist[i,j]==0){out.hist=out.hist+log(s[j])+log(1-p[j])}
    }

    out=out+counts.in[i]*out.hist

  }

  # fish not detected last period and no censures
  in.hist=c(1:dim(use.hist)[1])[(use.hist[,dim(use.hist)[2]]==0)&(apply(use.hist,1,function(x) {!(2 %in% x)}))]
  for(i in in.hist){
    last.detect=c(1:num.period)[use.hist[i,]==1]
    if(length(last.detect)==0){last.detect=0}else{last.detect=max(last.detect)}
    out.hist=1-l*L.in[num.period]
    if(!((num.period-1)<(last.detect+1))){
      for (k in (num.period-1):(last.detect+1)){out.hist=(1-L.in[k])+L.in[k]*(1-s[k]+s[k]*(1-p[k])*out.hist)}
    }
    out.hist=log(out.hist)
    if(last.detect>0){
      out.hist=out.hist+log(L.in[last.detect]) # p(tag active at last detected site)
      for (k in last.detect:1){
        if(use.hist[i,k]==1){out.hist=out.hist+log(s[k])+log(p[k])+log(1-d[k])}
        if(use.hist[i,k]==0){out.hist=out.hist+log(s[k])+log(1-p[k])}
      }
    }
    out=out+counts.in[i]*out.hist
  }

  # fish not detected last period and censured at some point
  in.hist=c(1:dim(use.hist)[1])[(use.hist[,dim(use.hist)[2]]==0)&(apply(use.hist,1,function(x) {(2 %in% x)}))]
  for(i in in.hist){
    out.hist=1
    last.detect=c(1:num.period)[use.hist[i,]==2]
    out.hist=log(out.hist)
    for (k in last.detect:1){
      if(use.hist[i,k]==2){out.hist=out.hist+log(s[k])+log(p[k])+log(d[k])+log(L.in[k])}
      if(use.hist[i,k]==1){out.hist=out.hist+log(s[k])+log(p[k])+log(1-d[k])}
      if(use.hist[i,k]==0){out.hist=out.hist+log(s[k])+log(1-p[k])}
    }
    out=out+counts.in[i]*out.hist
  }
  return(-out)
}
