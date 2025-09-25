#' @title Estimates the length of time a tag is active (tag activation to first detection time)
#'
#' @description This function uses the length of time a tag has been active to
#'   estimate the probability it would fail when detected. The average P(failure)
#'   is estimated at each site. The estimated survival is divided by the mean
#'   P(failure) at that site to adjust for the predicted tag failure in the
#'   study (Townsend et al. 2006).
#'
#' @param activetime.matrix Matrix of time elapsed after tag activation to first detection at site. One column per site
#' @param site.names Vector of site designations
#' @param num.period Number of detection sites
#' @param taglife.fit Results from fitting tag-life study tags: Model name, estimated parameters, mean time to fail
#' @param num.bootstrap Number of Bootstrap interations for variance on estimated P(Li) for each site. Uses 0 for initial fittings
#' @param cjs.est 2-column Matrix with unadjusted Cormack-Jolly-Seber estimates and standard errors
#'
#' @importFrom failCompare fc_pred
#'
#' @return returns a list with:
#' \describe{
#'    \item{model.out}{name of model used to fit tag-life curve}
#'    \item{L}{vector of probabilty tag active to each detection site}
#'    \item{L.se}{vector of standard errors on each estimated L}
#'    \item{adj.Si.se}{vector of estimated standared errors for adjusted survival estimates}
#'    }
#'
#' @export
cjs.taglife.corr=function(activetime.matrix,site.names=NULL,num.period=num.period,taglife.fit=taglife.fit,num.bootstrap=0,cjs.est=NULL){
#print("cjs.taglife.corr")
 model.used=taglife.fit$mod_choice

  var.ratio.fn=function(Ti,Tj,var.Ti,var.Tj,cov.TiTj){
    # B.25 ATLAS 1.4 manual
    out=((Ti/Tj)^2)*((var.Ti/(Ti^2))+(var.Tj/(Tj^2))-2*(cov.TiTj/(Ti*Tj)))
    return(out)
  }

  L.out = rep(1,num.period)
  names(L.out)=site.names

  for (i in 1:num.period) {
    temp.time = as.numeric(activetime.matrix[!is.na(activetime.matrix[, i]),i])

    fail.times=try(failCompare::fc_pred(mod_obj = taglife.fit,times=temp.time))
    L.out[i] = mean(fail.times)
  }

  if(num.bootstrap>0){
    stopifnot('Estimation of var(Li) requires CJS estimates' = !is.null(cjs.est))

    est.p=cjs.est[-c(1:(num.period-1)),1] # only need pi and lambda
    est.s=cjs.est[c(1:(num.period-1),dim(cjs.est)[1]),]

    boot.Ls=boot.L(at.time.matrix=activetime.matrix,model.in=taglife.fit,num.bootstrap=num.bootstrap)

    L.cov1 = stats::cov(boot.Ls$L.matrix) # only taglife study resampled
    L.cov2 = stats::cov(boot.Ls$L2.matrix) # taglife and active time resampled

    L.var1 = diag(L.cov1)
    L.var2 = diag(L.cov2)

    L.var= (1-est.p)*(L.var2-L.var1)+L.var1 # B.26 ATLAS 1.4 manual
    L.se = sqrt(L.var)
    names(L.se)=names(L.out)

    # estimate variance on Si including var(Ti) and covar(Li,L(i-1)) B.25 ATLAS 1.4 manual
    varL.ratio=rep(NA,num.period)
    var.Si=rep(NA,num.period)

    varL.ratio[1]=var.ratio.fn(1,L.out[1],0,L.var[1],0) # L_0=1, var(L_0)=0, cov(L0,L1)=0
    var.Si[1]=(est.s[1,1]^2)*varL.ratio[1]+(1/L.out[1]^2)*(est.s[1,2]^2)-varL.ratio[1]*(est.s[1,2]^2)


   for(i in 2:num.period){
      cov.LiLj=(1-est.p[i-1])*(1-est.p[i])*(L.cov2[i-1,i]-L.cov1[i-1,i])+L.cov1[i-1,i]
      varL.ratio[i]=var.ratio.fn(L.out[i-1],L.out[i],L.var[i-1],L.var[i],cov.LiLj)
      var.Si[i]=(est.s[i,1]^2)*varL.ratio[i]+((L.out[i-1]/L.out[i])^2)*(est.s[i,2]^2)-varL.ratio[i]*(est.s[i,2]^2)
    }

    adj.Si.se=sqrt(var.Si)
#	print(adj.Si.se)


  }else{L.se=NULL; adj.Si.se=NULL}

  return(list(model.out=model.used,L=L.out,L.se=L.se,adj.Si.se=adj.Si.se))
}


