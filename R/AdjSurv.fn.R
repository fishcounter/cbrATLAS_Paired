
#' @title Active Tag-Life-Adjusted Survival Modeling Fitting
#'
#' @description This function analyzes release-recapture data from a study that uses active-tag technology.
#' It estimates survival and detection probabilities and, if provided with tag-life data, the average probability
#' a tag is active at each detection site. These estimates are then used to adjust estimated survival for
#' potential tag-failure that may otherwise be interpreted as mortality. Methods are documented in
#' \href{http://www.cbr.washington.edu/sites/default/files/manuals/ATLAS_1_4_Manual_0.pdf}{Program ATLAS 1.4:Active Tag Life Adjusted Survival}.
#'
#' @param taglife.file Optional. Name of .csv file with taglife times in first column. Other columns ignored. Header expected
#' @param taghist.file  Required. Dataframe or a filename of .csv file with detection histories
#' @param taghist.format Required. Format of tag detection histories
#' \describe{
#'   \item{"atlas"}{format based on Program ATLAS input files. Eight columns. Each tag has one line per possible detection site.}
#'   \item{"flat"}{format has one line per tag: tag.code, activation date, release date, and
#' 	          one column per detection site filled with first detection times at that site}
#' }
#' @param taglife.model fc_obj. Failure time model object. (default = NULL).  If no fc_obj is provided, function will try models available in failCompare and select best fit
#' @param num.release (default = 1) If more than one release group, first column will be added to flat format file to denote group name (not implemented in this version)
#' @param num.bootstrap (default = 1000) Number of resample iterations to estimate additional variance on survival estimates
#' @param adjust.cjs (T|F) (default = T) adjust CJS estimates for estimated tag-life
#' @param rounding  (default = 4)  Number of decimal places on estimate
#' @param plot.taglife (T|F) (default = T) Plot the estimated tag-life curve
#'
#' @importFrom failCompare fc_fit fc_rank fc_select
#' @importFrom utils installed.packages packageVersion
#'
#' @return Returns a list "Out" with (if provided):
#' \describe{
#'   \item{taghist}{tagfile filename}
#'   \item{unadjusted.cjs}{table of unadjusted survival estimates}
#'   \item{tagfile}{tag life filename}
#'   \item{taglife.model}{taglife model name}
#'   \item{mean.tag.pLIve}{expect proportion tags active to each site}
#'   \item{adjust.cjs}{adjusted survival estimates}
#'   }
#'
#'
#' @details All date-times are assumed to be character vector with format: "%m/%d/%Y %H:%M"
#'
#'
#' @export
#'
AdjSurv.fn=function (taglife.file=NULL, taghist.file,taghist.format="atlas", taglife.model=NULL, num.release=1,
                     num.bootstrap = 250, adjust.cjs = T, rounding = 4, plot.taglife = T) {
  if(is.data.frame(taghist.file)){data.file=taghist.file}else{
    data.file = data.frame(utils::read.csv(taghist.file, header = F, colClasses = c("character")))}
  if(!is.null(taglife.file)){
    if(is.data.frame(taglife.file)){tag.life=taglife.file}else{
    tag.life = utils::read.csv(taglife.file, header = T)}
    tag.life = sort(tag.life[, 1])
  }

    # reformat detetion file if in Program ATLAS format from ATLAS to flat
  if(taghist.format == "atlas"){data.file=atlas2flat.fn(data.file)}

  rel.groups=unique(data.file[,1])
  stopifnot('Number of release groups in file does not match specified number releases' = (length(rel.groups)==num.release))

    data.in=data.file[,-c(1:2)] # assumes rel.group and bin columns in 1st 2 columns of data table

    # strip out detection history into "detects" matrix
    num.period = (dim(data.in)[2] - 3)/2
    site.names = names(data.in)[4:(3 + num.period)]
    detects = data.in[, 4:(3 + num.period)] # capture histories only
    for (i in 1:dim(detects)[2]) {
      detects[detects[, i] == "", i] = 0
      detects[, i] = factor(detects[, i])
      }
    detects=apply(detects,2,as.numeric) # convert dataframe to numeric matrix to speed up code

    # Estimate unadjusted CJS
    unadj.cjs.params=cjs.fn(detect.in=detects,se.out=T)

    if(!is.null(taglife.file)){
		  fc.loaded=("failCompare" %in% installed.packages())
		  fc.up2date=!(packageVersion("failCompare")<"1.0.0") #check if fc pkg up to date
		if((!fc.up2date)|(!fc.loaded)){print("Please download the latest version of the failCompare pkg at http://www.cbr.washington.edu/analysis/apps/failcompare/")}
     ################################################### Estimate tag-life curve
		if(is.null(taglife.model)){
	    mod_ls=failCompare::fc_fit(time=taglife.file$tag_life_days,model=c('weibull','weibull3','gompertz','gamma','lognormal','llogis','gengamma'))
			mod_ls_ranked=failCompare::fc_rank(mod_ls)

			# currently, select best-fitting model automatically
			taglife.fit=failCompare::fc_select(mod_ls_ranked,model = as.character(mod_ls_ranked$"GOF_tab"[1,1]))
			}else{taglife.fit=taglife.model}

		if(plot.taglife){plot(taglife.fit)}

      ########### Estimate mean travel time to each site, use to calculate mean failure at each site
      tt.rel2site=mean_tt2site.fn(data.in,num.period,site.names)

      mean.tag.p=cjs.taglife.corr(activetime.matrix=tt.rel2site$activetime.matrix,site.names=site.names,
                                num.period=num.period,taglife.fit=taglife.fit,num.boots=0)

      ### once tag-life model has been selected, rerun to get bootstapped se's on P(Li)
      mean.tag.p=cjs.taglife.corr(activetime.matrix=tt.rel2site$activetime.matrix,site.names=site.names,
                                num.period=num.period,taglife.fit=taglife.fit,num.boots=num.bootstrap,cjs.est=unadj.cjs.params$cjs.param)


      # Estimate adjusted CJS
      if(adjust.cjs){
        adj.cjs.params=cjs.fn(detect.in=detects,L.in=mean.tag.p$L,seeds.in=unadj.cjs.params$cjs.param[,1],se.out=T,d=unadj.cjs.params$d)
        adj.cjs.params$cjs.param[c(1:(num.period-1),num.period*2-1),2]=mean.tag.p$adj.Si.se
      }
    }


    ### create output list from analysis
    if(is.data.frame(taghist.file)){file.out=deparse(substitute(taghist.file))}else{
      file.out=taghist.file
      }
   out=list(taghist=file.out, unadjusted.cjs=unadj.cjs.params)
    if(!is.null(taglife.file)){
      if(is.data.frame(taglife.file)){out$tagfile=deparse(substitute(taglife.file))}else{
      out$tagfile=taglife.file}
      out$taglife.model=taglife.fit
      out$mean.tag.pLive=mean.tag.p[2:3]
      if(adjust.cjs){out$adjusted.cjs=adj.cjs.params[1]}
    }

  return(out)
}
