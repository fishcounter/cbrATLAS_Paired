#' @title A likelihood function for multiple unadjusted- or adjusted-for-taglife Cormack-Jolly-Seber (CJS) estimates 
#'	that have common parameters (S, p, l)
#'
#' @description Based on the likelihood function for unadjusted- or adjusted-for-taglife
#'  Cormack-Jolly-Seber (CJS) estimates and Skalski et al. (1998) 
#'  Adjustment to survival for estimated tag life as described in Townsend et al.
#'  (2006).  User input indicates which parameters should be set to the same value. This
#'  function is set above individual CJS likelihood calls
#'
#' @param params List of likelihood MLE estimators (S,p,l) in a single vector (1 each for upper, lower release)
#' @param counts.in List of Vectors of summed detections by history (names are history)
#' @param num.period Number of detection periods in history
#' @param use.hist List of Matrices of unique detection histories
#' @param L.in List of probability that a tag is working at a detection site at a given time. Set to 1 if active tags not indicated. 
#'        p(tag working) is cumulative to each site, not a multiplicate from site to site
#' @param d.in List of probability a tag is censored given that it's detected
#' @param common.start numeric indicating point survival/detection p's are equated
#'
#' @return This function returns a negative log-likelihood
#' @export
#'
paired.cjs.lik=function(params.in,counts.in,num.period,use.hist,L.in=NULL,d.in=NULL,common.start=NULL,like.check=F,hessian=F,fix.1=NULL){
  # fix.1 = character vector naming the parameters to be fixed to 1 (e.g., fix.1=c("p11"))
# print("paired.cjs.lik")
	params=list(params.in[1:(length(params.in)/2)],params.in[(length(params.in)/2+1):length(params.in)])

 if(!is.null(fix.1))
  { 
    npars<-length(params.in) + length(fix.1)
    num.period <- (npars/2-1)/2+1
    fix.1<-sort(fix.1)
    
    # okay as long as number of sites < 11
    fix1.rel1<-grep("p.1",fix.1,value=TRUE)
    if(length(fix1.rel1)>0)
    {
      fr1.sites<-as.numeric(substr(fix1.rel1,2,2))
      for(i in fr1.sites) {params.in<-c(params.in[1:(2*(i-1)+1)],1,params.in[(2*(i-1)+2):length(params.in)])
						   names(params.in)[2*fr1.sites]<-fix1.rel1}
    } else fr1.sites<-c()
    
    #print(params.in)
    
    fix1.rel2<-grep("p.2",fix.1,value=TRUE)
    if(length(fix1.rel2)>0)
    {
      fr2.sites<-as.numeric(substr(fix1.rel2,2,2))
      st.site.rel2<- 2*num.period -1 + 1 # length(params.in)/2 + 1 + length(fr1.sites)
      #print(st.site.rel2)
      #print(mode(params.in))
      #print(fr2.sites)
      for(i in fr2.sites) params.in<-c(params.in[1:(st.site.rel2 + 2*(i-1))],1,params.in[(st.site.rel2 + 2*(i-1) + 1):length(params.in)])
      #print(params.in)
      names(params.in)[(2*num.period):length(params.in)][2*fr2.sites]<-fix1.rel2
    }
    #print(params.in)
    
  }
  
  
  params=list(params.in[1:(length(params.in)/2)],params.in[(length(params.in)/2+1):length(params.in)])
  # browser()
  #  if(is.null(L.in[[1]])){L.tmp=list(rep(1,num.period),rep(1,num.period))}else{print(L.in)} # both upper & lower must either be 0 or both tag corrected 
  #  print(round(params.in,10));print("start paired.cjs.lik")
  if(!(num.period>1)){cat("Number of periods must be greater than 1\n")}
  stopifnot(num.period>1)
  
  #  params=lapply(params.in,correct.fn)
  
  # set common parameters in lower release to equal upper release params
  if(is.null(common.start)){common.start=0} 
  if(common.start>0){
    #browser()
    common.param=c(seq(1,(length(params[[1]])-1),2),seq(2,(length(params[[1]])-1),2),length(params[[1]]))
    params[[2]][!(common.param<(common.start+1))]=params[[1]][!(common.param<(common.start+1))]
  }  
  
#print("paired.cjs.lik_rab: input to paired.cjs.lik")
# this undoes the equating of common parameters (still has the fixing to 1) - use for testing
#params=list(params.in[1:(length(params.in)/2)],
#            params.in[((length(params.in)/2)+1):length(params.in)])



  # if(like.check){browser()}
  
  # print(params)
  for(g in 1:2){ #1=upper release params 2=lower
    # split mle vector into separate types
    if(like.check){browser()}
    s=params[[g]][1:(num.period-1)]
    p=params[[g]][(num.period):((num.period-1)*2)]
    l=params[[g]][1+(num.period-1)*2]
    if(!is.null(L.in)){L.tmp=L.in[[g]]}
    if(is.null(d.in[[g]])){d.lik=rep(0,num.period-1)}else{d.lik=d.in[[g]]}
    
    # combinatoric
    if(g==1){out=lgamma(sum(counts.in[[g]])+1)-sum(lgamma(counts.in[[g]]+1))}else{
      out=out+lgamma(sum(counts.in[[g]])+1)-sum(lgamma(counts.in[[g]]+1))}
    
#cat("\n paired.cjs.lik_rab: out after combinatoric for rel", g, ": ", out, "\n")    
    
    # fish detected in last period
    for(i in c(1:dim(use.hist[[g]])[1])[use.hist[[g]][,dim(use.hist[[g]])[2]]==1]){
      out.hist=log(l)+log(L.tmp[num.period])
      for(j in 1:(num.period-1)){
        if(use.hist[[g]][i,j]==1){out.hist=out.hist+log(s[j])+log(p[j])+log(1-d.lik[j])}
        if(use.hist[[g]][i,j]==0){out.hist=out.hist+log(s[j])+log(1-p[j])}
      }
      
#cat("\n paired.cjs.lik_rab: first out after detected in last period for rel", g, ": ", out, "\n")  
      
      
      out=out+counts.in[[g]][i]*out.hist
      
#cat("\n paired.cjs.lik_rab: final out after detected in last period for rel", g, ": ", out, "\n")  
      
    }
    
    # fish not detected last period and no censures
    in.hist=c(1:dim(use.hist[[g]])[1])[(use.hist[[g]][,dim(use.hist[[g]])[2]]==0)&(apply(use.hist[[g]],1,function(x) {!(2 %in% x)}))]
    for(i in in.hist){
      last.detect=c(1:num.period)[use.hist[[g]][i,]==1]
      if(length(last.detect)==0){last.detect=0}else{last.detect=max(last.detect)}
      out.hist=1-l*L.tmp[num.period]
      if(!((num.period-1)<(last.detect+1))){
        for (k in (num.period-1):(last.detect+1)){out.hist=(1-L.tmp[k])+L.tmp[k]*(1-s[k]+s[k]*(1-p[k])*out.hist)}
      }
      
      out.hist=log(out.hist)
      if(last.detect>0){
        out.hist=out.hist+log(L.tmp[last.detect]) # p(tag active at last detected site)
        for (k in last.detect:1){
          if(use.hist[[g]][i,k]==1){out.hist=out.hist+log(s[k])+log(p[k])+log(1-d.lik[k])}
          if(use.hist[[g]][i,k]==0){out.hist=out.hist+log(s[k])+log(1-p[k])}
        }
      }
      out=out+counts.in[[g]][i]*out.hist
    }
    
    
    # fish not detected last period and censured at some point
    in.hist=c(1:dim(use.hist[[g]])[1])[(use.hist[[g]][,dim(use.hist[[g]])[2]]==0)&(apply(use.hist[[g]],1,function(x) {(2 %in% x)}))]
    for(i in in.hist){
      out.hist=1
      last.detect=c(1:num.period)[use.hist[[g]][i,]==2]
      out.hist=log(out.hist)
      
#cat("\n paired.cjs.lik_rab: out.hist from not detected in last period and censored at some point for rel", g, " before updating: ", out.hist, "\n") 
#cat("d.lik = ", d.lik,"\n")      

      for (k in last.detect:1){
        if(use.hist[[g]][i,k]==2){out.hist=out.hist+log(s[k])+log(p[k])+log(d.lik[k])+log(L.tmp[k])}
#cat("\n paired.cjs.lik_rab: out.hist from not detected in last period and censored at some point for rel", g, " and det.code = 2: ", out.hist, "\n") 
#cat("log(s[k]) = ",log(s[k]),"\n")
#cat("log(p[k]) = ",log(p[k]),"\n")
#cat("log(d.lik[k]) = ",log(d.lik[k]),"\n")
#cat("log(L.tmp[k]) = ",log(L.tmp[k]),"\n")

        
        if(use.hist[[g]][i,k]==1){out.hist=out.hist+log(s[k])+log(p[k])+log(1-d.lik[k])}
#cat("\n paired.cjs.lik_rab: out.hist from not detected in last period and censored at some point for rel", g, " and det.code = 1: ", out.hist, "\n")     
        
        if(use.hist[[g]][i,k]==0){out.hist=out.hist+log(s[k])+log(1-p[k])} 
#cat("\n paired.cjs.lik_rab: out.hist from not detected in last period and censored at some point for rel", g, " and det.code = 0: ", out.hist, "\n")     
        
      }
#cat("\n paired.cjs.lik_rab: out.hist from not detected in last period and censored at some point for rel", g, ": ", out.hist, "\n")     
#cat("\n paired.cjs.lik_rab: counts.in[[g]][i] from not detected in last period and censored at some point for rel", g, ": ", counts.in[[g]][i], "\n")     

      out=out+counts.in[[g]][i]*out.hist
    }
    
#cat("\n paired.cjs.lik_rab: final out after not detected in last period and censored at some point for rel", g, ": ", out, "\n")  
    
  }
  
#cat(paste0("paired.cjs.lik_rab: params = ",params, "; function value = ", -out,"\n"))

  return(-out)
}

