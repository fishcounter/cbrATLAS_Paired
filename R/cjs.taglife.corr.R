#' @title Estimate time tag active (activation to 1st detection time)
#'
#' @param activetime.matrix matrix of time elapsed after tag activation to 1st detection at site. 1 column per site
#' @param site.names vector of site designations
#' @param num.period number of detection sites
#' @param taglife.fit results from fitting tag-life study tags. model name, estimated parameters, mean time to fail
#' @param num.boots bootstrap for variance on estimated P(Li) for each site? Use 0 for initial fittings
#' @param cjs.est 2-column matrix with unadjusted CJS estimates and standard errors
#'
#' @return model.out: name of model fit to tag-life curve
#' @export
#'
cjs.taglife.corr=function(activetime.matrix,site.names=NULL,num.period=num.period,taglife.fit=taglife.fit,num.boots=0,cjs.est=NULL){
  # estimate time tag active (activation to 1st detection time)
  # activetime.matrix: matrix of time elapsed after tag activation to 1st detection at site. 1 column per site
  # site.names: vector of site designations
  # num.period: number of detection sites
  # taglife.fit: results from fitting tag-life study tags. model name, estimated parameters, mean time to fail
  # num.boots: bootstrap for variance on estimated P(Li) for each site? Use 0 for initial fittings.
  # cjs.params: 2-column matrix with unadjusted CJS estimates and standard errors
  # OUTPUT
  # model.out: name of model fit to tag-life curve
  # L: [if num.boots=0] vector of P(Li), one for each detection site
  # L.se: [if num.boots=0] NULL
  #    [if num.boots>0] matrix of 2 columns (L, L.se), and 1 row for each detection site.
  # boot.L.matrix: [if num.boots=0] NULL
  #  [if num.boots>0] matrix of P(Li) columns and num.boots rows using bootstrapped taglife values only.

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

    fail.times=fc_pred(times=temp.time,pars=taglife.fit$mod_objs[,1],model=model.used)
    L.out[i] = mean(fail.times)
  }

  if(num.boots>0){
    stopifnot('Estimation of var(Li) requires CJS estimates' = !is.null(cjs.est))
    est.p=cjs.est[-c(1:(num.period-1)),1] # only need pi and lambda
    est.s=cjs.est[c(1:(num.period-1),dim(cjs.est)[1]),]

    boot.Ls=boot.L(at.time.matrix=activetime.matrix,model.in=taglife.fit,num.boots=num.boots)
    L.cov1 = cov(boot.Ls$L.matrix) # only taglife study resampled
    L.cov2 = cov(boot.Ls$L2.matrix) # taglife and active time resampled

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

  }else{L.se=NULL; adj.Si.se=NULL}

  return(list(model.out=model.used,L=L.out,L.se=L.se,adj.Si.se=adj.Si.se))
}

boot.L=function(at.time.matrix,model.in,num.boots){
  # bootstrap "activation to 1st detection per site" to get standard errors on estimated p(Li)
  # INPUT
  # at.time.matrix: matrix of time from activation to detection at site i
  # model.in: output from taglife.fn
  # num.boots: number of desired resampling bootstraps to estimes se on p(Li)
  # OUTPUT
  # L.matrix: matrix of bootstrapped Li w/resampled taglife tags
  # L2.matrix: matrix of bootstrapped Li 2/resampled taglife tags and active times to detection
  num.period=dim(at.time.matrix)[2]

  L.matrix=matrix(0,nrow=num.boots,ncol=num.period)
  L2.matrix=matrix(0,nrow=num.boots,ncol=num.period)

  # estimate s^2[Ti|d] based on resampled taglife study and observed active time
  for (i in 1:num.boots) {
    t.in=sort(sample(model.in$times[,1],replace=T))
    x=fc_fit(time=t.in,model=model.in$mod_choice)

    for (j in 1:dim(at.time.matrix)[2]){
      at.temp=as.numeric(at.time.matrix[!is.na(at.time.matrix[, j]),j])
      L.matrix[i,j]=mean(fc_pred(time=at.temp,pars=x$mod_objs[,1],model=model.in$mod_choice))
    }
  }

  # estimate s^2[Ti] based on resampled taglife study and resampled observed active time
  for (i in 1:num.boots) {
    t.in=sort(sample(model.in$times[,1],replace=T))
    x=fc_fit(time=t.in,model=model.in$mod_choice)


    for (j in 1:dim(at.time.matrix)[2]){
      at.sample=sample(1:dim(at.time.matrix)[1],replace=T)
      at.tmp=at.time.matrix[at.sample,]
      at.temp=as.numeric(at.tmp[!is.na(at.tmp[, j]),j])
      L2.matrix[i,j]=mean(fc_pred(time=at.temp,pars=x$mod_objs[,1],model=model.in$mod_choice))
    }
  }

  return(list(L.matrix=L.matrix,L2.matrix=L2.matrix))
}


