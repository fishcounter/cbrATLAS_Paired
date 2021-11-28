
#' @title Active Tag-Life-Adjusted Survival Modeling Fitting
#'
#' @param taglife.file Optional. Name of .csv file with taglife times in 1st column.  other columns ignored.  Header expected.
#' @param taghist.file  Required. dataframe or a filename of .csv with detection histories. 2 formats currently supported.
#' @param taghist.format Format of tag detection histories: "atlas" or"flat"
#' "atlas" = format based on Program ATLAS input files. 8 columns.  Each tag has 1 line per possible detection site
#' 	"flat" = format has 1 line per tag. tag.code, activation date, release date, and
#' 	          1 column per detection site filled with 1st detection times at that site (current format: "%m/%d/%Y %H:%M")
#' @param taglife.model Failure time model obs (fc_obj)
#' @param num.release if more than 1 release group 1st column will be added to flat format file that will denote group name
#' @param num.release if more than 1 release group 1st column will be added to flat format file that will denote group name (currently not working)
#' @param num.bootstrap number of resamples to estimate additional variance on survival estimates, default is 1000
#' @param adjust.cjs Should CJS estimates be adjusted for estimated tag-life? Logical.
#' @param rounding Number of decimal places on estimates
#' @param plot.taglife plot the estimated tag-life curve. Logical.
#'
#' @return Tables of estimates for unadjusted survival model and, if the taglife adjustment is used, taglife model estimates and adjusted survival estimates
#'
#' @export
#'
AdjSurv.fn=function (taglife.file=NULL, taghist.file,taghist.format="atlas", taglife.model=NULL, num.release=1,
                     num.bootstrap = 1000, adjust.cjs = T, rounding = 4, plot.taglife = T) {
                     num.bootstrap = 250, adjust.cjs = T, rounding = 4, plot.taglife = T) {
  if(is.data.frame(taghist.file)){data.file=taghist.file}else{
    data.file = data.frame(read.csv(taghist.file, header = F, colClasses = c("character")))}
  if(!is.null(taglife.file)){
    if(is.data.frame(taglife.file)){tag.life=taglife.file}else{
    tag.life = read.csv(taglife.file, header = T)}
    tag.life = sort(tag.life[, 1])
  }
  library(failCompare)

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
    ################################################### Estimate tag-life curve
      if(is.null(taglife.model)){
        mod_ls=fc_fit(time=taglife.data$tag_life_days,model=c('weibull','weibull3','gompertz','gamma','lognormal','llogis','gengamma'))
        mod_ls_ranked=fc_rank(mod_ls)

        # currently, select best-fitting model automatically
        taglife.fit=fc_select(mod_ls_ranked,model = as.character(mod_ls_ranked$"GOF_tab"[1,1]))
      }else{taglife.fit=taglife.model}

      if(plot.taglife){plot(taglife.fit)}

      ########### Estimate mean travel time to each site, use to calculate mean failure at each site
      tt.rel2site=mean.tt2site.fn(data.in,num.period,site.names)

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
