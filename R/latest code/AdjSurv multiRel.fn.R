
#' @title Active Tag-Life-Adjusted Survival Modeling FittingFG
#'
#' @description This function analyzes release-recapture data from a study that uses active-tag technology.
#' It estimates survival and detection probabilities and, if provided with tag-life data, the average probability
#' a tag is active at each detection site. These estimates are then used to adjust estimated survival for
#' potential tag-failure that may otherwise be interpreted as mortality. Methods are documented in
#' \href{http://www.cbr.washington.edu/sites/default/files/manuals/ATLAS_1_4_Manual_0.pdf}{Program ATLAS 1.4:Active Tag Life Adjusted Survival}.
#'
#' @param taghist.file  Required. Dataframe or a filename of .csv file with detection histories
#' @param taghist.format Required. Format of tag detection histories
#' \describe{
#'   \item{"atlas"}{format based on Program ATLAS input files. Eight columns. Each tag has one line per possible detection site.}
#'   \item{"flat"}{format has one line per tag: tag.code, activation date, release date, and
#' 	          one column per detection site filled with first detection times at that site}
#' }
#' @param taglife.file Optional. Name of .csv file with taglife times in 1st column.  other columns ignored.  Header expected.
#' @param taglife.model fc_obj. Failure time model object. (default = NULL).  If no fc_obj is provided, function will try models available in failCompare and select best fit
#' @param num.release (default = 1) If more than 1 release group, 1st column will be added to flat format file to  denote group name
#' @param num.bootstrap (default = 1000) Number of resample iterations to estimate additional variance on survival estimates
#' @param adjust.cjs (T|F) (default = T) adjust CJS estimates be adjusted for estimated tag-life
#' @param rounding  (default = 4)  Number of decimal places on estimate
#' @param plot.taglife (T|F) (default = T) Plot the estimated tag-life curve
#' @param study.type analyze ("single"|"paired") releases.
#'
#' @importFrom failCompare fc_fit fc_rank fc_select
#' @importFrom utils installed.packages packageVersion
#' @importFrom devtools install_github
#' @importFrom remotes install_remote
#'
#' @return Returns a list "out" with (if provided):
#' \describe{
#'   \item{taghist}{tagfile filename}
#'   \item{unadjusted.cjs}{table of unadjusted survival estimates}
#'   \item{tagfile}{tag life filename}
#'   \item{taglife.model}{taglife model name}
#'   \item{mean.tag.pLive}{expect proportion tags active to each site}
#'   \item{adjust.cjs}{adjusted survival estimates}
#'   }
#'
#' @export
#' 
AdjSurv.fn=function(taglife.file=NULL, taghist.file, taghist.format="atlas", taglife.model=NULL, num.release=1,
                     num.bootstrap = 1000, adjust.cjs = T, rounding = 4, plot.taglife = T, study.type="single") {
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
	  print("num.bootstrap")
	  print(num.bootstrap)
	  
      tt.rel2site=mean.tt2site.fn(data.in,num.period,site.names)

      mean.tag.p=cjs.taglife.corr(activetime.matrix=tt.rel2site$activetime.matrix,site.names=site.names,
                                num.period=num.period,taglife.fit=taglife.fit,num.bootstrap=0)
      ### once tag-life model has been selected, rerun to get bootstapped se's on P(Li)
      mean.tag.p=cjs.taglife.corr(activetime.matrix=tt.rel2site$activetime.matrix,site.names=site.names,
                                  num.period=num.period,taglife.fit=taglife.fit,num.bootstrap=num.bootstrap,cjs.est=unadj.cjs.params$cjs.param)

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

    out=list(release=rel.groups[rel.i])  
	out[[rel.i]]$taghist=file.out

	out[[rel.i]]$capture.history=capt.hist
	out[[rel.i]]$unadjusted.cjs=unadj.cjs.params

    if(!is.null(taglife.file)){
      out[[rel.i]]$tagfile=taglife.file
      out[[rel.i]]$taglife.model=taglife.fit
#	  out[[rel.i]]$taglife.gof=taglife.gof # pre-selecting a taglife model does not provide a GOF
#	  out[[rel.i]]$
      out[[rel.i]]$mean.tag.pLive=mean.tag.p[2:3]

      if(adjust.cjs){out[[rel.i]]$adjusted.cjs=adj.cjs.paired.params[1]}
		}
	}
	if(!is.null(taglife.file)){result.list=append(result.list,out)}
  }
  
  ################################################################################################################## start of paired release
  # paired release study
  if(study.type=="paired"){

	# define releases
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
		}
	}

	common.start=as.numeric(common.start)
    # Estimate unadjusted paired release

   unadj.cjs.paired.params=cjs.paired.fn(detects.in=detects,L.in=NULL,seeds.in=NULL,se.out=T,d=NULL,common.start.in=as.numeric(common.start))
	if(common.start>0){	
		# correct lower release estimates to match common parameters used in maximum likelihood in paired.cjs.lik
		num.params=dim(unadj.cjs.paired.params$paired.cjs.param)[1]/2
		unadj.cjs.paired.params$paired.cjs.param[(num.params+1):(2*num.params),][!((1:num.params)<(common.start+1)),1:2]=
			unadj.cjs.paired.params$paired.cjs.param[1:num.params,][!((1:num.params)<(common.start+1)),1:2]
		} 	 

#	unadj S ratio
	length.cjs.est=dim(unadj.cjs.paired.params$paired.cjs.param)[1]
	unadj.S.top=   unadj.cjs.paired.params$paired.cjs.param[1,]
	unadj.S.bottom=unadj.cjs.paired.params$paired.cjs.param[(length.cjs.est/2)+1,]

	unadj.S.ratio=vs2.fn(unadj.S.top[1],unadj.S.bottom[1],unadj.S.top[2],unadj.S.bottom[2])

    if(!is.null(taglife.file)){
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
			  
	  # mean.tag.p no boots
	  mean.tag.p[[rel.i]]=cjs.taglife.corr(activetime.matrix=tt.rel2site$activetime.matrix,
											site.names=site.names,
                                            num.period=num.period,
											taglife.fit=taglife.fit,
											num.bootstrap=0) 

	  ### once tag-life model has been selected, rerun to get bootstapped se's on P(Li)
	  # need to adjust estimates by release
	  unadj.cjs.length=dim(unadj.cjs.paired.params[[1]])[1]
	  if(rel.i==1){unadj.cjs.in=c(1:(unadj.cjs.length/2))}else{	   
				   unadj.cjs.in=c(((unadj.cjs.length/2)+1):max(unadj.cjs.length))
				   }
	  # mean.tag.p w/boots
	  mean.tag.p[[rel.i]]=cjs.taglife.corr(activetime.matrix=tt.rel2site$activetime.matrix,
											site.names=site.names,
											num.period=num.period,
											taglife.fit=taglife.fit,
											num.bootstrap=num.bootstrap,	
											cjs.est=unadj.cjs.paired.params$paired.cjs.param[unadj.cjs.in,])
	  }			   

############ apply tag-life correction #########################################
cat("ADJUSTING FOR TAGLIFE!")
      if(adjust.cjs){

		adj.cjs.paired.params=cjs.paired.fn(detects.in=detects,
										L.in=list(mean.tag.p[[1]]$L,mean.tag.p[[2]]$L),
										seeds.in=unadj.cjs.paired.params$paired.cjs.param[,1],
										se.out=T,
										num.bootstrap=1000,
										common.start.in=as.numeric(common.start))

		if(common.start>0){	
			# correct lower release estimates to match common parameters used in maximum likelihood in paired.cjs.lik
			num.params=dim(adj.cjs.paired.params$paired.cjs.param)[1]/2
			adj.cjs.paired.params$paired.cjs.param[(num.params+1):(2*num.params),][!((1:num.params)<(common.start+1)),1:2]=
				adj.cjs.paired.params$paired.cjs.param[1:num.params,][!((1:num.params)<(common.start+1)),1:2]
			} 	 

        adj.cjs.paired.params$cjs.param[c(1:(num.period-1),num.period*2-1),2]=mean.tag.p$adj.Si.se 
      }

#	adjusted S ratio                                                  
	length.cjs.est=dim(adj.cjs.paired.params$paired.cjs.param)[1]
	adj.S.top=   adj.cjs.paired.params$paired.cjs.param[1,]
	adj.S.bottom=adj.cjs.paired.params$paired.cjs.param[(length.cjs.est/2)+1,]
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
					 out$adjusted.cjs=adj.cjs.paired.params
					 adj.S.paired=adj.S.ratio}
	}

    result.list=append(result.list,out)
    }  
  ################################################################################################################## end of paired release

  return(result.list) # saves list to an assigned object
}


