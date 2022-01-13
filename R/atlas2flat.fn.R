
#' @title Converts ATLAS data file to a flat file format.
#'
#' @description The ATLAS file format (detailed in the ATLAS 1.4 manual, page 17) is a vertical file
#' with one line per detection site.  The input format for the cbrATLAS uses a flat file format, which details
#' the history of eacg tag on one line, including release and dection times at each site.
#'
#' @param data.in Table: columns = release group name, bin number, tag id, tag activation date/time, tag release date/time, site name, detection (1 = yes, 0 = no), detection date/time
#'
#' @return Creates a table with columns: release group name, bin number, tag id, activation date/time, tag release date/time, 1 column per site name with detection time (blank if no detection)
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
