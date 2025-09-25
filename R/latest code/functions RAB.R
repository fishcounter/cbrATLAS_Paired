
RABs edits to functions
Here are the functions that I edited. They are all in a single file, so the line numbers I use here all refer to this same file. All of these functions have "_rab" appended to the end of their names.
 
Changes I made:

# AdjSurv.fn_rab:  I added in arguments for the optim options.
# Note that I did not edit the code for the single release at all.

# adj.cjs.params: 
	I set the bootstrap argument to num.bootstrap in implementing cjs.paired.fn_rab().
DONE

# cjs.paired.fn_rab():
	I added in arguments for the optim options.
	I commented out the line: if(!is.null(L.in)){browser}
	Moved the lines that define the list "d" outside the brackets for if(is.null(seeds.in)).
	I added a second to remove parameters from seeds.in if they are fixed to 1. Currently applies only to p parameters and assumes that all site indices are < 10. This is lines 420-442.
	added argument "fix.1" to implementation of paired.cjs.lik_rab() to deal with the fixed-to-1 parameters - both times paired.cjs.like_rab() is run.
	Added in lines to add back the fixed-to-1 parameters to the output - see lines 569-582.

# paired.cjs.like_rab():
	Added argument "fix.1" to deal with parameters that are fixed to 1.
	Added code for fixed-to-1 parameters in lines 604-635.

# cjs.taglife.corr_rab():
	updated it to use boot.L_rab() instead of boot.L().

# boot.L_rab():
	Changed the "for" statement to a "while" statement (including increasing i at the end of the "while" code).
	Defined x using tryCatch() on failCompare::fc_fit.
	Changed the implementation of failCompare::fc_pred to use mod_obj=x in defining L.matrix and L2.matrix (lines 855 and 881).
 
 
 
# functions RAB.R
AdjSurv.fn_rab=function (taglife.file=NULL, taghist.file,taghist.format="atlas", taglife.model=NULL, num.release=1,
                     num.bootstrap = 1000, adjust.cjs = T, rounding = 4, plot.taglife = T, study.type="single",
                     parscale.in = 1e-04, reltol.in = 1e-30, maxit.in =2e+06) 
{
  # helper fns
  summaryf = function(x){cbind(summary(factor(x)))}
  sum.paste=function(x,sep=""){summaryf(apply(x,1,paste,collapse=sep))}
  
  data.file=NULL
  if(is.data.frame(taghist.file)){data.file=taghist.file}else{
    data.file = data.frame(read.csv(taghist.file, header = F, colClasses = c("character")))}
  if(!is.null(taglife.file)){
    if(is.data.frame(taglife.file)){tag.life=taglife.file}else{
      tag.life = read.csv(taglife.file, header = T)}
    tag.life = sort(tag.life[,1])
  }
  
  # reformat detetion file if in Program ATLAS format from ATLAS to flat
  if(taghist.format == "atlas"){data.file=atlas2flat.fn(data.file)}
  
  rel.groups=unique(data.file[,1])
  stopifnot('Number of release groups in file does not match specified number releases' = (length(rel.groups)==num.release))
  
  
  result.list=list(NULL)
  out=NULL
  detects=NULL
  capt.hist=list(NULL)
  #  if(num.release>1){result.list=NULL} # why only if num.releases>1?
  
  # single release study
  if(study.type=="single"){
    for(rel.i in 1:length(rel.groups)){
      print(paste0("Release Group ",rel.i,", ID Name: ",rel.groups[rel.i], collapse=""))
      
      data.in=data.file[data.file[,1] %in% rel.groups[rel.i],-c(1:2)] # assumes rel.group and bin columns in 1st 2 columns of data table
      
      # strip out detection history into "detects" matrix
      num.period = (dim(data.in)[2] - 3)/2
      site.names = names(data.in)[4:(3 + num.period)]
      detects = data.in[, 4:(3 + num.period)] # capture histories only
      for (i in 1:dim(detects)[2]) {
        detects[detects[, i] == "", i] = 0
        detects[, i] = factor(detects[, i])
      }
      detects=apply(detects,2,as.numeric) # convert dataframe to numeric matrix to speed up code
      capt.hist[[rel.i]]=sum.paste(detects)
      # Estimate unadjusted CJS
      unadj.cjs.params=cjs.fn(detect.in=detects,se.out=T)
      
      if(!is.null(taglife.file)){
        ################################################### Estimate tag-life curve
        require(failCompare)
        if(is.null(taglife.model)){
          mod_ls=failCompare::fc_fit(time=tag.life,model=c("weibull","weibull3", "gompertz", "gamma", "lognormal", "llogis", "gengamma","vitality.ku","vitality.4p"))
          mod_ls_ranked=fc_rank(mod_ls)
          
          # currently, select best-fitting model automatically
          taglife.fit=fc_select(mod_ls_ranked,model = as.character(mod_ls_ranked$"GOF_tab"[1,1]))
          
        }else{taglife.fit=failCompare::fc_fit(time=tag.life,model=taglife.model)
        }
        
        if(plot.taglife&(rel.i==1)){plot(taglife.fit)}
        
        ########### Estimate mean travel time to each site, use to calculate mean failure at each site
        tt.rel2site=mean.tt2site.fn(data.in,num.period,site.names)
        
        mean.tag.p=cjs.taglife.corr_rab(activetime.matrix=tt.rel2site$activetime.matrix,site.names=site.names,
                                    num.period=num.period,taglife.fit=taglife.fit,num.boots=0)
        ### once tag-life model has been selected, rerun to get bootstapped se's on P(Li)
        
        mean.tag.p=cjs.taglife.corr_rab(activetime.matrix=tt.rel2site$activetime.matrix,site.names=site.names,
                                    num.period=num.period,taglife.fit=taglife.fit,num.boots=num.bootstrap,cjs.est=unadj.cjs.params$cjs.param)
        
        # Estimate adjusted CJS
        if(adjust.cjs){
          adj.cjs.params=cjs.fn(detect.in=detects,L.in=mean.tag.p$L,seeds.in=unadj.cjs.params$cjs.param[,1],se.out=T,d=unadj.cjs.params$d)
          
          adj.cjs.params$cjs.param[c(1:(num.period-1),num.period*2-1),2]=mean.tag.p$adj.Si.se
        }
      }
      
      
      ### create output list from analysis
      if(is.data.frame(taghist.file)){file.out=deparse(substitute(taghist.file))}else{ # recheck this code
        file.out=taghist.file
      }
      
      out=list(release=rel.groups[rel.i])  #maybe try initializing / [[NULL]]?
      out[[rel.i]]$taghist=file.out
      
      out[[rel.i]]$capture.history=capt.hist
      out[[rel.i]]$unadjusted.cjs=unadj.cjs.params
      
      if(!is.null(taglife.file)){
        out[[rel.i]]$tagfile=taglife.file
        out[[rel.i]]$taglife.model=taglife.fit
        #	  out[[rel.i]]$taglife.gof=taglife.gof # pre-selecting a taglife model does not provide a GOF
        #	  out[[rel.i]]$
        out[[rel.i]]$mean.tag.pLive=mean.tag.p[2:3]
        
        if(adjust.cjs){out[[rel.i]]$adjusted.cjs=adj.cjs.params[1]}
      }
    }
    if(!is.null(taglife.file)){result.list=append(result.list,out)}
  }
  
  # paired release study
  if(study.type=="paired"){
    
    # define releases
    {
      rel.upper=0; rel.lower=0
      print(paste0(c(1,rel.groups[1]),collapse = ": "))
      print(paste0(c(2,rel.groups[2]),collapse = ": "))
      print("",quote=F)
      while(!((rel.upper %in% c(1,2))&(rel.lower %in% c(1,2))&(!(rel.upper == rel.lower)))){
        rel.upper=readline(prompt="select upper release (number): ")
        rel.lower=readline(prompt="select lower release (number): ")}
      
      
      rel.upper=as.numeric(rel.upper,n=1)
      rel.lower=as.numeric(rel.lower,n=1)
      rel.groups=rel.groups[c(rel.upper,rel.lower)]
    }
    
    # process detection data
    {
      detects.tmp=list()
      
      for(rel.i in 1:2){
        data.in=data.file[data.file[,1] %in% rel.groups[rel.i],-c(1:2)] # assumes rel.group and bin columns in 1st 2 columns of data table
        
        # strip out detection history into "detects" 
        num.period = (dim(data.in)[2] - 3)/2
        site.names = names(data.in)[4:(3 + num.period)]
        detects.tmp = data.in[, 4:(3 + num.period)] # capture histories only
        for (i in 1:dim(detects.tmp)[2]) {
          detects.tmp[detects.tmp[, i] == "", i] = 0
          detects.tmp[, i] = factor(detects.tmp[, i])
        }
        # each pooled rel group detections 1 element in list, in order of location
        detects[[rel.i]]=apply(detects.tmp,2,as.numeric) # convert dataframe to numeric matrix to speed up code
        capt.hist[[rel.i]]=sum.paste(detects[[rel.i]]) # output detection summaries
      }
    }
    
    # select point when common survival/detection probs start
    {
      print("Select first parameter to begin common survival & detection probabilities between release groups",quote=F)
      print("1 = first detection site common to both releases",quote=F)
      print("2 = 2nd reach (survival), 3 is 2nd detection site, etc.",quote=F)
      print("Last number represents P(survival to & detection at last site)",quote=F)
      common.start=(-1)
      while((common.start<0)|(common.start>((num.period-1)*2))){ # don't continue until legitimate common parameter entered
        print("Enter '0' for full CJS model (all params unique)",quote=F)
        common.start=readline(prompt=noquote(paste0(c("Enter number from 0 to ",(num.period-1)*2,": "),collapse="")))
        #		print(paste0(c("you have selected: ",common.start)))
      }
    }
    
    common.start=as.numeric(common.start)
    # Estimate unadjusted paired release
    
#print(paste0("AdjSurv.fn_rab: common.start = ",common.start))  

cat("\n AdjSurv.fn_rab: Computing the unadjusted CJS estimates for paired release study design \n")
    
    unadj.cjs.paired.params=cjs.paired.fn_rab(detects.in=detects,L.in=NULL,seeds.in=NULL,se.out=T,d=NULL,common.start.in=as.numeric(common.start))
    
print("AdjSurv.fn_rab: Have the unadjusted cjs params")    
#print(paste0("AdjSurv.fn_rab: unadj.cjs.paired.params = ",unadj.cjs.paired.params))    
    
#return(unadj.cjs.paired.params)

    #	unadj S ratio
    length.cjs.est=dim(unadj.cjs.paired.params$paired.cjs.param)[1]
    unadj.S.top=   unadj.cjs.paired.params$paired.cjs.param[1,]
    unadj.S.bottom=unadj.cjs.paired.params$paired.cjs.param[(length.cjs.est/2)+1,]
    
    unadj.S.ratio=vs2.fn(unadj.S.top[1],unadj.S.bottom[1],unadj.S.top[2],unadj.S.bottom[2])
    
    
    if(!is.null(taglife.file))
      {
      ################################################### Estimate tag-life curve     Need to make r code match code in D:/ATLAS pkg/cbrATLAS/R/AdjSurv.fn.R
      require(failCompare)
      if(is.null(taglife.model)){
        mod_ls=failCompare::fc_fit(time=tag.life,model=c("weibull","weibull3", "gompertz", "gamma", "lognormal", "llogis", "gengamma","vitality.ku","vitality.4p"))
        mod_ls_ranked=fc_rank(mod_ls)
        
        # currently, select best-fitting model automatically
        print(paste0(c("The taglife curve that best fits uses the ",as.character(mod_ls_ranked$"GOF_tab"[1,1])),collapse=""))
        
        taglife.fit=fc_select(mod_ls_ranked,model = as.character(mod_ls_ranked$"GOF_tab"[1,1]))
        
      }else{taglife.fit=failCompare::fc_fit(time=tag.life,model=taglife.model)
      }
      # taglife.fit=taglife.model
      
      if(plot.taglife){plot(taglife.fit)}
      
      mean.tag.p = list()
      
      for(rel.i in 1:2){
        data.in=data.file[data.file[,1] %in% rel.groups[rel.i],-c(1:2)] # assumes rel.group and bin columns in 1st 2 columns of data table
        
        ########### Estimate mean travel time to each site, use to calculate mean failure at each site
        
        tt.rel2site=mean.tt2site.fn(data.in,num.period,site.names)
        
cat("\n AdjSurv.fn_rab: ",paste0("rel ",rel.i,": "), "calculate mean prob tag survival to each site no_bootstrap \n")

        # mean.tag.p no boots
        mean.tag.p[[rel.i]]=cjs.taglife.corr_rab(activetime.matrix=tt.rel2site$activetime.matrix,
                                             site.names=site.names,
                                             num.period=num.period,
                                             taglife.fit=taglife.fit,
                                             num.boots=0) 
        
        ### once tag-life model has been selected, rerun to get bootstapped se's on P(Li)
        # need to adjust estimates by release
        unadj.cjs.length=dim(unadj.cjs.paired.params[[1]])[1]
        
        if(rel.i==1) {unadj.cjs.in=c(1:(unadj.cjs.length/2))}else{	  
          
          unadj.cjs.in=c(((unadj.cjs.length/2)+1):max(unadj.cjs.length))
          
        }
        
cat("\n AdjSurv.fn_rab: ",paste0("rel ",rel.i,": "), "calculate mean prob tag survival to each site w_bootstrap \n")
        
        # mean.tag.p w/boots
        mean.tag.p[[rel.i]]=cjs.taglife.corr_rab(activetime.matrix=tt.rel2site$activetime.matrix,
                                             site.names=site.names,
                                             num.period=num.period,
                                             taglife.fit=taglife.fit,
                                             num.boots=num.bootstrap,	
                                             cjs.est=unadj.cjs.paired.params$paired.cjs.param[unadj.cjs.in,])
      }			   
      
#cat("\n AdjSurv.fn_rab: mean prob tag survival to each site w_bootstrap = \n")
#print(mean.tag.p)
      
      ############ apply tag-life correction #########################################
      
      if(adjust.cjs){
        print(common.start)
        
cat("\n AdjSurv.fn_rab: Computing adjusted CJS paired parameters w_SE \n")  

        adj.cjs.params=cjs.paired.fn_rab(detects.in=detects,
                                         L.in=list(mean.tag.p[[1]]$L,mean.tag.p[[2]]$L),
                                         seeds.in=unadj.cjs.paired.params$paired.cjs.param[,1],
                                         se.out=T,
                                         bootstrap=num.bootstrap,#1000,
                                         common.start.in=as.numeric(common.start),
                                         parscale.in = parscale.in, reltol.in = reltol.in, maxit.in = maxit.in)
        
        adj.cjs.params$cjs.param[c(1:(num.period-1),num.period*2-1),2]=mean.tag.p$adj.Si.se 
        
#print(mean.tag.p$adj.Si.se)
print("AdjSurv.fn_rab: adjusted CJS paired parameter estimates")
#print(adj.cjs.params)

      }
      
      #	adjusted S ratio                                                  
      length.cjs.est=dim(adj.cjs.params$paired.cjs.param)[1]
      adj.S.top=   adj.cjs.params$paired.cjs.param[1,]
      adj.S.bottom=adj.cjs.params$paired.cjs.param[(length.cjs.est/2)+1,]
      adj.S.ratio=vs2.fn(adj.S.top[1],adj.S.bottom[1],adj.S.top[2],adj.S.bottom[2])
    }
    
    ### create output list from analysis
    if(is.data.frame(taghist.file)){file.out=deparse(substitute(taghist.file))}else{
      file.out=taghist.file
    }
    out=append(out,list(release=rel.groups))
    
    out=list(release=rel.groups,
             taghist=file.out, 
             capture.history=capt.hist, 
             unadjusted.cjs=unadj.cjs.paired.params,
             unadj.S.paired=unadj.S.ratio,
             adjusted.cjs=NULL,
             adj.S.paired=NULL)
    
    if(!is.null(taglife.file)){
      if(is.data.frame(taglife.file)){out$tagfile=deparse(substitute(taglife.file))}else{
        out$tagfile=taglife.file}
      out$taglife.model=taglife.fit
      out$mean.tag.pLive=mean.tag.p
      out$adj.ratio=adj.S.ratio
      
      if(adjust.cjs){
        #				 out$adjusted.cjs=adj.cjs.params
        adj.S.paired=adj.S.ratio}
    }
    
    result.list=append(result.list,out)
  }  
  #print(result.list) # for immediate answers
  return(result.list) # saves list to an assigned object
}

#####################################################################################################################################################


cjs.paired.fn_rab=function(detects.in,L.in=NULL,seeds.in=NULL,se.out=F,d=NULL,bootstrap=F,common.start.in=0, parscale.in = 1e-04, reltol.in = 1e-30, maxit.in =2e+06)
{
  
#cat("\n HERE WE ARE IN cjs.paired.fn_rab \n\n")
  
  if(is.null(seeds.in)){seeds.in=list(NULL,NULL)}else{
    seeds.tmp=list(NULL,NULL)
    seeds.tmp[[1]]=seeds.in[1:(length(seeds.in)/2)]
    seeds.tmp[[2]]=seeds.in[((length(seeds.in)/2)+1):length(seeds.in)]
    seeds.in=seeds.tmp}
  
#L.in.defined<-!(is.null(L.in)) 
#print(L.in.defined)
  
  counts.in=list(NULL,NULL)
  obs.hist=list(NULL,NULL)
  d=list(NULL,NULL)
  #if(!is.null(L.in)){browser}
  if(is.null(L.in)){L.in=list(NULL,NULL)}
  
  for(g in 1:2){
    # call each release separately (1=upper release; 2=lower release)
    # common lower release likelihood inputs will be set to upper release corresponding upper release likelihood inputs 
    
    detect.in=detects.in[[g]]
    num.period=dim(detect.in)[2]
    if(is.null(L.in[[g]])){L.in[[g]]=rep(1,num.period)} # if no tag failure, L = 100% for all sites
    
    # summary of counts.in/history and unique detection matrix in same order
    counts.in[[g]] = thist0(detect.in)$count
    obs.hist[[g]] = thist0(detect.in)$hist.matrix
    
    if(is.null(seeds.in[[g]])){
      

#if(L.in.defined) print("L.in IS DEFINED - WE ARE ADJUSTING BY TAGLIFE NOW for default seeds \n") else print("L.in IS NULL and seeds.in is NULL ALSO")

#cat("\n cjs.paired.fn_rab: seeds.in is null so define the default values based on detection data before adjusting by L.in for L.in", ifelse(L.in.defined,"","NOT") ,"< 1 \n \n")
            
      # seed estimates do not estimate survival/detection with censoring present.
      seeds.in[[g]] = rep(0.5, num.period * 2 - 1)
      names(seeds.in[[g]])=c(paste0(rep(c("S","p"),rep(num.period-1,2)),rep(1:(num.period-1),2),g),paste0("l",g))
      
      
      if (num.period > 2) {
        seeds.in[[g]][num.period] = sum(detect.in[apply(detect.in[,-1]==1, 1, sum) > 0, 1] == 1)/sum(apply(detect.in[, -1] == 1, 1, sum) > 0)
        for (i in (2:(num.period - 2))) {
          seeds.in[[g]][num.period + i - 1] = sum(detect.in[apply(detect.in[,-c(1:i)] == 1, 1, sum) > 0, i] == 1)/
            sum(apply(detect.in[,-c(1:i)] == 1, 1, sum) > 0)
        }
        seeds.in[[g]][num.period * 2 - 2] = sum(detect.in[detect.in[,num.period] == 1, num.period - 1] == 1)/sum(detect.in[,num.period] == 1)
      }else{seeds.in[[g]][num.period] = sum(detect.in[detect.in[, 2] == 1,1] == 1)/sum(detect.in[, 2] == 1)}
      
      seeds.in[[g]][(2 * num.period) - 1] = sum(detect.in[, num.period] == 1)/(sum(detect.in[, num.period - 1] == 1)/seeds.in[[g]][2 * (num.period - 1)])
      seeds.in[[g]][1] = (sum(detect.in[, 1] == 1)/dim(detect.in)[1])/seeds.in[[g]][num.period]
      if (num.period > 2) {
        for (i in 2:(num.period - 1)) {
          seeds.in[[g]][i] = (sum(detect.in[, i] == 1)/seeds.in[[g]][i + num.period - 1])/(sum(detect.in[, i - 1] == 1)/seeds.in[[g]][i + num.period - 2])                
        }
      }
      
      
      
#cat("\n cjs.paired.fn_rab: default seeds.in before adjusting by L.in for L.in", ifelse(L.in.defined,"","NOT") ,"< 1:", unlist(seeds.in), "\n\n")
      
    }
    
    # print("cjs.paired.fn_rab/p(censored)")
    d[[g]]=rep(0,num.period-1)
    for(i in (1:(num.period-1))){
      d[[g]][i]= sum(detect.in[,i]==2)/sum(detect.in[,i]>0)
    }
    #			print(d[[g]])
    names(d[[g]])=c(paste0(rep(paste0(rep("d",rep(num.period-1)))),rep(g,rep(num.period-1)),c(1:(num.period-1))))
    
    
    # set common parameters
    #		print("default seeds: ");print(seeds.in)
    if(common.start.in>0){
      if(g==2){seeds.in[[2]][!(c(1:length(seeds.in[[2]]))<(common.start.in+1))]=seeds.in[[1]][!(c(1:length(seeds.in[[2]]))<(common.start.in+1))]
      }}
    
#cat("\n cjs.paired.fn_rab: seeds.in for rel", g," before adjusting by taglife = ", unlist(seeds.in), "\n\n")
    
    # adjust S,lambda by expected tag-life.  L.in
    seeds.in[[g]][c(1:(num.period - 1), 2 * num.period - 1)] = seeds.in[[g]][c(1:(num.period - 1), 2 * num.period - 1)]/L.in[[g]]
    
#cat("\n cjs.paired.fn_rab: seeds.in for rel", g," after adjusting by taglife = ", unlist(seeds.in), "\n\n")
    
  }
  
#print(paste0("cjs.paired.fn_rab: seeds.in = ",seeds.in))  
#print("cjs.paired.fn_rab: seeds.in after all adjustment by taglife")
#print(seeds.in)
  
  ## maximum likelihood estimate of CJS parameters
  #parscale.in=1e-04
  #reltol.in=1e-30
  #maxit.in=2e+06

  
  ### print inputs 	
  # print("entering paired.cjs.lik 1st time");print(paste0(c("common.start.in= ",as.numeric(common.start.in))),collapse="")
  # if common.start.in is > 0, set lower release to equal upper release starting at first common site to end
  # print("seeds.in");print(seeds.in);print("counts.in");print(counts.in);print("common.start.in");print(common.start.in)
  # print(paste0("common.start.in 2: ",common.start.in))
  #browser()


# Remove from seeds.in the parameters to fix to 1 (if any)                                   STOPPED HERE!!
{
  cat("num.period: ", num.period, "\n")
  par.names<-c(paste0("S",1:(num.period-1)),paste0("p",1:(num.period-1)),"l")
  
  #cat("par.names: ", par.names, "\n")
  
  if(num.period>2) par.names<-c(sapply(c("S","p"), \(x) {sapply(1:(num.period-1),\(xx) paste0(x,xx))}),"l")

  #cat("par.names: ", par.names, "\n")
  
  pn.1<-paste0(par.names,"1") #gsub("S","S1",gsub("p","p1",gsub("l","l1",par.names)))  
  pn.2<-paste0(par.names,"2") #gsub("S","S2",gsub("p","p2",gsub("l","l2",par.names)))  
  par.names<-c(pn.1,pn.2)
  cat("par.names: ", par.names, "\n")
  
  seeds.in <- unlist(seeds.in)
  
  #cat("length of seeds.in: ", length(seeds.in), "\n")
  if(sum(seeds.in==1)>0) {fx.1<-par.names[seeds.in==1]; seeds.in<-seeds.in[seeds.in!=1]} else fx.1<-NULL
  cat("fx.1: ", fx.1, "\n")
  cat("seeds.in after fx.1: ", seeds.in, "\n")
}
browser()

#cat("\n cjs.paired.fn_rab: computing CJS parameters for paired release no_SE\n")
#cat("\n cjs.paired.fn_rab: computing CJS parameters for paired release no_SE this is d.in: ", unlist(d), "\n")

#init_fn_val<-paired.cjs.lik_rab(params.in = unlist(seeds.in),
#                                counts.in = counts.in,
#                                num.period = num.period,
#                                use.hist = obs.hist,
#                                L.in = L.in,
#                                d.in = d,
#                                common.start = as.numeric(common.start.in),
#                                like.check = F,
#                                hessian = F,
#                                fix.1 = fx.1)
#cat("\n cjs.paired.fn_rab: initial function value with L.in = ",init_fn_val, "\n\n")
  
  paired.cjs.out = optim(par=unlist(seeds.in), fn=paired.cjs.lik_rab, 
                         counts.in = counts.in, 
                         num.period=num.period, 
                         use.hist=obs.hist,
                         L.in=L.in, # was L=L.in
                         d.in=d,
                         common.start=as.numeric(common.start.in),
                         like.check=F, 
                         method = "BFGS", 
                         hessian = F, 
                         control = list(trace=F,
                                        parscale = rep(parscale.in, length(unlist(seeds.in))), 
                                        reltol = reltol.in,
                                        maxit=maxit.in,
                                        ndeps=rep(parscale.in, length(unlist(seeds.in)))),
                         fix.1=fx.1)
  
#print(paste0("cjs.paired.fn_rab: paired.cjs.out = ", paired.cjs.out))  
  
  # browser()
  # print("after paired.cjs.out first time")
  # print(paired.cjs.out)
  # print(paste0("common.param: ",common.param))
  #	print("back into cjs.paired.fn_rab")
  #	print(paired.cjs.out$par)		
  #	print("set common estimates")
  #	print(paired.cjs.out$par[c(1:(length(paired.cjs.out$par)/2))][!(common.param<(common.start.in+1))])
  #
  #	paired.cjs.out$par[length(paired.cjs.out$par)/2+c(1:(length(paired.cjs.out$par)/2))][!(common.param<(common.start.in+1))]=
  #			paired.cjs.out$par[c(1:(length(paired.cjs.out$par)/2))][!(common.param<(common.start.in+1))] # set 2nd release common params  = 1st release common params
  #		}  
  #print("Done with paired cjs w equal params")		
  paired.cjs.out$par=correct.fn(paired.cjs.out$par)
  
#print(paste0("cjs.paired.fn_rab: paired.cjs.out$par no_SE = ", paired.cjs.out$par))
  
  ##################
  
  #print(cjs.out$par)
  if(se.out){
    parscale.in=1e-04
    reltol.in=1e-10 #1e-30
    maxit.in=2000 #2e+06
    
    # rerun to get hessian
cat("\n cjs.paired.fn_rab: estimate Hessian for CJS estimates for paired release \n")    
    paired.cjs.out = optim(paired.cjs.out$par, 
                           fn=paired.cjs.lik_rab, 
                           counts.in = counts.in, 
                           num.period = num.period, 
                           use.hist=obs.hist,
                           L=L.in,
                           d.in=d,
                           common.start=common.start.in,
                           like.check=F,
                           method = "BFGS", 
                           hessian = T, 
                           control = list(parscale = rep(parscale.in, length(unlist(seeds.in))), 
                                          reltol = reltol.in,
                                          maxit=maxit.in,
                                          trace=1,
                                          ndeps=rep(parscale.in, length(unlist(seeds.in)))),
                           fix.1 = fx.1)
    
    paired.cjs.out$par=correct.fn(paired.cjs.out$par)
    
#print(paste0("cjs.paired.fn_rab: paired.cjs.out$par pt_est_w_Hessian = ", paired.cjs.out$par))    
    
    # calculate standard errors on CJS estimates, based on Hessian
    temp = paired.cjs.out$hessian
    singular = colSums(temp) == 0
    
#cat("cjs.paired.fn_rab: singular = ",singular,"\n")    
    
    #var.vec = rep(0, 2*(num.period * 2 - 1))
    var.vec <- rep(0, length(seeds.in))
    var.vec[!singular] = diag(solve(temp[!singular, !singular]))
    
#print("temp:")
#print(temp)
#print("var.vec:")
#print(var.vec)
    
    
    
    # format estimates into a matrix, 2 columns. 1 row per parameter
    paired.cjs.param = matrix(c(paired.cjs.out$par, sqrt(var.vec)), ncol = 2)

    #if(common.start.in>0){  ##########   in theory, the likelihood fn will estimate common params
    #			common.param=c(seq(1,(length(params[[1]])-1),2),seq(2,(length(params[[1]])-1),2),length(params[[1]]))
    #			print(common.param)
    #		paired.cjs.param[dim(paired.cjs.param)[1]/2+c(1:(dim(paired.cjs.param)[1]/2)),][!(common.param<(common.start.in+1)),]=
    #				paired.cjs.param[c(1:(dim(paired.cjs.param)[1]/2)),][!(common.param<(common.start.in+1)),] # set 2nd release common params  = 1st release common params
    #	}  
    
    # note: after paired.cjs.param	(line 112), estimates match ATLAS, but missing last site SE (rel2 lambda SE)
    
    row.names(paired.cjs.param)=setdiff(names(seeds.in),fx.1)
    dimnames(paired.cjs.param)[[2]] = c("Estimate", "s.e.")
  }else{ print("cjs.paired.fn_rab: why the heck is it popping up?")
    # though not used in current version, this option can be used to get point estimates only
    # format estimates into a matrix, 1 columns. 1 row per parameter
    paired.cjs.param = matrix(c(paired.cjs.out$par), ncol = 1)
    row.names(paired.cjs.param)=setdiff(names(unlist(seeds.in)),fx.1)
  }

#print(paired.cjs.param)    
#print(is.matrix(paired.cjs.param))

# add in parameters fixed to 1 if necessary
if(!is.null(fx.1))
{
  npars.est<-nrow(paired.cjs.param)
  paired.cjs.param<-rbind(paired.cjs.param,matrix(rep(c(1,0),each=length(fx.1)),ncol=2))
  row.names(paired.cjs.param)[(npars.est+1):(npars.est+length(fx.1))]<-fx.1
  
#print(paired.cjs.param)
  
#print(par.names)
  paired.cjs.param<-paired.cjs.param[par.names,]
  
#print(paired.cjs.param)
}


  return(list(paired.cjs.param=paired.cjs.param, d=d))
}







#####################################################################################################################################################



paired.cjs.lik_rab=function(params.in,counts.in,num.period,use.hist,L.in=NULL,d.in=NULL,
                            common.start=NULL,like.check=F,hessian=F,fix.1=NULL){
  
  # fix.1 = character vector naming the parameters to be fixed to 1 (e.g., fix.1=c("p11"))
  
  # helper fn
  bounding.fn=function(x){
    if(x<0){out=(-x%%1)}
    if(x>1){out=(1-x%%1)}
    if((!(x<0)) & (!(x>1))){out=x}
    return(out)}
  
  #cat("\n paired.cjs.lik_rab: params.in = ", params.in,"\n")
  
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
      for(i in fr1.sites) params.in<-c(params.in[1:(2*(i-1)+1)],1,params.in[(2*(i-1)+2):length(params.in)])
      names(params.in)[2*fr1.sites]<-fix1.rel1
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


#print(paste0("paired.cjs.lik_rab: params = ",params))

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
    
#cat("\n paired.cjs.lik_rab: final out after not detected in last period and not censored for rel", g, ": ", out, "\n")  
    
#IT IS IN THE NEXT PART THAT IT GETS INF FOR SECOND DATA SET
    
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

#####################################################################################################################################################


cjs.taglife.corr_rab<-function(activetime.matrix,site.names=NULL,num.period=num.period,taglife.fit=taglife.fit,num.boots=0,cjs.est=NULL){
  
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
  
  if(num.boots>0){
    stopifnot('Estimation of var(Li) requires CJS estimates' = !is.null(cjs.est))
    
#cat("\n cjs.taglife.corr_rab with num.boots>0: cjs.est = ", cjs.est, "\n")

    est.p=cjs.est[-c(1:(num.period-1)),1] # only need pi and lambda
    est.s=cjs.est[c(1:(num.period-1),dim(cjs.est)[1]),]
    
    boot.Ls=boot.L_rab(at.time.matrix=activetime.matrix,model.in=taglife.fit,num.boots=num.boots)
    #	print(c(mean(activetime.matrix[,3],na.rm=T),mean(activetime.matrix[,4],na.rm=T)))
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
    #print("var.Si");print(var.Si)
    adj.Si.se=sqrt(var.Si)
    #print("leaving cjs.taglife.corr!)")
    
  }else{L.se=NULL; adj.Si.se=NULL}
  
  return(list(model.out=model.used,L=L.out,L.se=L.se,adj.Si.se=adj.Si.se))
}


#####################################################################################################################################################

boot.L_rab<-function(at.time.matrix,model.in,num.boots=1000){
  num.period=dim(at.time.matrix)[2]
  
  L.matrix=matrix(0,nrow=num.boots,ncol=num.period)
  L2.matrix=matrix(0,nrow=num.boots,ncol=num.period)
  
  i<-1
  # estimate s^2[Ti|d] based on resampled taglife study and observed active time
  while(i <= num.boots) {
    
    #cat("i=",i," ")
    
    t.in=sort(sample(model.in$times[,1],replace=T))
    #x=failCompare::fc_fit(time=t.in,model=model.in$mod_choice)
    
    ii<-NULL
    x<-tryCatch(failCompare::fc_fit(time=t.in, model=model.in$mod_choice),
                 error = function(e) {ii <<- i-1},
                 warning = function(w) {ii <<- i-1}
                 )
    if(!is.null(ii)) i<-(i-1) else {
    
    for (j in 1:dim(at.time.matrix)[2]){
      at.temp=as.numeric(at.time.matrix[!is.na(at.time.matrix[, j]),j])
      L.matrix[i,j]=mean(failCompare::fc_pred(mod_obj=x,times=at.temp))
    }
    }
    
    i<-i+1
  }
  
  rm(i,ii,t.in,x,at.temp)
  
  i<-1
  # estimate s^2[Ti] based on resampled taglife study and resampled observed active time
  while(i <= num.boots) {
    t.in=sort(sample(model.in$times[,1],replace=T))
    #x=failCompare::fc_fit(time=t.in,model=model.in$mod_choice)
    ii<-NULL
    x<-tryCatch(failCompare::fc_fit(time = t.in,
                                    model = model.in$mod_choice),
                error = function(e) {ii <<- i-1},
                warning = function(w) {ii <<- i-1}
                )
    if(!is.null(ii)) i<-(i-1) else {
    
    for (j in 1:dim(at.time.matrix)[2]){
      at.sample=sample(1:dim(at.time.matrix)[1],replace=T)
      at.tmp=at.time.matrix[at.sample,]
      at.temp=as.numeric(at.tmp[!is.na(at.tmp[, j]),j])
      L2.matrix[i,j]=mean(failCompare::fc_pred(mod_obj=x,times=at.temp))
    }
    }
    
    i<-i+1
  }
  
  return(list(L.matrix=L.matrix,L2.matrix=L2.matrix))
}


