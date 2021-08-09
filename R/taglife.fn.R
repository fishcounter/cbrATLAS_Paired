
#' @title Fitting tag-life function
#'
#' @param tags.in vector of observed time to failure (days)
#' @param model.in name of model to use.  Current option is "weibull"
#' @param tag.se compute SEs of estimates
#' @param plot.taglife create figure to see how model fits
#'
#' @return List with MLE estimates, mean time to failure based on selected model, andvector of tags used in the estimation of the tag-life curve.
#' @export
#'
taglife.fn=function(tags.in,model.in="weibull",tag.se=T,plot.taglife){
  lik.in=NULL
  ######################### start loop for fitting tag-life curve after changing model or censoring tags from tail-end of tag-life study

  #### currently set to default to weibull, need to add prompt to choose model here
  if(model.in=="weibull"){lik.in=logweib.lik; par.in=c(min(tags.in),.1,max(tags.in)) }
  stopifnot('Tag life model name not valid' =!is.null(lik.in)) #if no model picked

  # likelihood
  answer <- optim(par=par.in,fn=lik.in,hessian=F,tags.in=tags.in,control = list(reltol = 1e-10,maxit=500,ndeps=rep(.01,3)))
  if(tag.se){
    answer <- optim(par=answer$par,fn=lik.in,method="BFGS",hessian=T,tags.in=tags.in,control = list(parscale = abs(1/answer$par), reltol = 1e-20,maxit=500,ndeps=abs(1/answer$par)^2)) # estimate Hessian if se(tag life params wanted)
  }
  if(model.in=="weibull"){
    params.out=data.frame(matrix(NA,nrow=3,ncol=3))

    params.out[,1]=c("Beta:", "Gamma:","Eta")
    params.out[,2]=round(abs(answer$par),4)
    names(params.out)=c("Parameter","Estimate","s.e.")
  }

  if(tag.se){
    # calculate standard errors on parameter estimates, based on Hessian
    temp = answer$hessian
    singular = colSums(temp) == 0
    var.vec = rep(0, 3)
    var.vec[!singular] = diag(solve(temp[!singular, !singular]))
    params.out[,3]=apply(cbind("(",round(sqrt(var.vec),4),")"),1,paste,collapse="")
  }else{params.out[,3]=NA}


  if(plot.taglife){
    left=1-(1:length(tags.in))/length(tags.in)

    if(dev.cur()>1){dev.off()};dev.new(width=6, height=4.5)
    par(mar=c(4,3,2,1)+.3,lab=c(15,10,1),las=1,tck=-.01,mgp=c(2,.4,0))
    plot(tags.in,left,xlab="Time (days) to failure",
         ylab="Percent Tags Surviving",lwd=2,cex=1.1,pch=3,
         xlim=c(0,max(tags.in)+6),ylim=c(0,1))

    tag.pred=0:(round(max(tags.in))+5)
    fail.pred=fail.fn(tag.pred,params.out[1,2],params.out[2,2],params.out[3,2])
    tag.dead=(fail.pred<.0001) # cut-off for max predicted tag-life
    tag.pred=tag.pred[!tag.dead]
    fail.pred=fail.pred[!tag.dead]

    lines(tag.pred,fail.fn(tag.pred,params.out[1,2],params.out[2,2],params.out[3,2]),lwd=3,col=4)

    mtext(paste(c("Mean Time to Failure = ",mttf(params.out[,2])),collapse=""),outer=F,line=.5)
  }
  ##### end of fitting loop.  need to add prompt that fit is accepted.  If not, go back to model choice
  ##### fn will return last model fit
  out=list(model.out=model.in, lik.out=-answer$value,params=params.out,meantime2fail=mttf(params.out$Estimate),tags.out=tags.in)

  return(out)
}

#' Title
#'
#' @param x Log-likelihood of 3-parameter weibull distribution (beta, gamma, eta)
#' @param tags.in observed time to failure
#'
#' @return negative log-likelihood
#' @export
#'
logweib.lik=function(x,tags.in){
  x=abs(x)
  b.in=x[1] # beta
  g.in=x[2] # gamma
  n.in=x[3] # eta
  if (sum(g.in>tags.in)>0){g.in=min(tags.in)-.0000001}

  l.out=log(b.in)-(b.in*log(n.in))+((b.in-1)*log(tags.in-g.in))-(((tags.in-g.in)/n.in)^b.in)
  l.out=-sum(l.out)
  return(l.out)
}


#' @title Predicted survival fraction for 3-parameter weibull distribution
#'
#' @param x time
#' @param b.in beta
#' @param g.in gamma
#' @param n.in eta
#'
#' @return vector of proportion active at time x
#' @export
#'
fail.fn=function(x,b.in,g.in,n.in){
  out=exp(-((x-g.in)/n.in)^b.in)
  out[is.na(out)]=1
  return(out)
}
