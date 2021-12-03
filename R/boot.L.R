#' @title bootstrap "activation to 1st detection per site" to get standard errors on estimated p(Li)
#'
#' @description To better estimate the variance on the adjusted-for-tag-failure survival estimates,
#' aa bootstrap is conducted on both the data used estimate tag-life and on the observed times the study tags were active.
#' Statistical methods are described in the ATLAS 1.4 manual, Appendix B2.
#'
#' @param at.time.matrix matrix of time from activation to detection at site i
#' @param model.in output from taglife.fn
#' @param num.boots number of desired resampling bootstraps to estimate the standard error for each taglife estimate
#'
#' @import failCompare
#'
#' @return list of L.matrix matrix of bootstrapped Li w/resampled taglife tags;L2.matrix matrix of bootstrapped Li w/resampled taglife tags and active times to detection
#' @export
#'
boot.L=function(at.time.matrix,model.in,num.boots){
    num.period=dim(at.time.matrix)[2]

  L.matrix=matrix(0,nrow=num.boots,ncol=num.period)
  L2.matrix=matrix(0,nrow=num.boots,ncol=num.period)

  # estimate s^2[Ti|d] based on resampled taglife study and observed active time
  for (i in 1:num.boots) {
    t.in=sort(sample(model.in$times[,1],replace=T))
    x=failCompare::fc_fit(time=t.in,model=model.in$mod_choice)

    for (j in 1:dim(at.time.matrix)[2]){
      at.temp=as.numeric(at.time.matrix[!is.na(at.time.matrix[, j]),j])
      L.matrix[i,j]=mean(failCompare::fc_pred(mod_obj=model.in,times=at.temp))
    }
  }

  # estimate s^2[Ti] based on resampled taglife study and resampled observed active time
  for (i in 1:num.boots) {
    t.in=sort(sample(model.in$times[,1],replace=T))
    x=failCompare::fc_fit(time=t.in,model=model.in$mod_choice)


    for (j in 1:dim(at.time.matrix)[2]){
      at.sample=sample(1:dim(at.time.matrix)[1],replace=T)
      at.tmp=at.time.matrix[at.sample,]
      at.temp=as.numeric(at.tmp[!is.na(at.tmp[, j]),j])
      L2.matrix[i,j]=mean(failCompare::fc_pred(mod_obj=model.in,times=at.temp))
    }
  }

  return(list(L.matrix=L.matrix,L2.matrix=L2.matrix))
}
