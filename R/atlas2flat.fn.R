
#' @title Converts ATLAS data file to a flat file format
#'
#' @description cbrATLAS allows the mark-recapture data to be loaded in one of two file formats. The ATLAS file format 
#' 	(detailed in the \href{https://www.cbr.washington.edu/analysis/apps/atlas}{ATLAS 1.9 manual}, Appendix A) is a vertical file
#' 	with one line per detection site. The default input format for the cbrATLAS package uses a flat file format, which details
#' 	the history of each tag on one line, including release and detection times at each site.  
#' 	This function will converts an ATLAS formatted file to a flat file.
#'
#' @param data.in Table: columns = release group name, bin number, tag id, tag activation date/time, tag release date/time, site name, 
#' 	detection (1 = yes, 0 = no), detection date/time
#'
#' @return Returns a table with columns: release group name, bin number, tag id, activation date/time, tag release date/time, one column 
#' 	per site name with detection time (blank if no detection)
#' @export
#'
atlas2flat.fn=function(data.in){
  names(data.in)[7:8]=c("detect","time")
  data.out = stats::reshape(data.in, v.names = c("detect","time"), idvar = c("V1","V2","V3","V4","V5"),
                     timevar = "V6", direction = "wide")
  data.out=data.out[,c(1:5,grep("detect",names(data.out)),grep("time",names(data.out)))]
  names(data.out)=sub("detect.","",names(data.out))
  return(data.out)
}
