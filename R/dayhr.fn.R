#' @title Standardize date-times input
#'
#' @description Convert date-time into yyyy-mm-dd hh:mm:ss format
#'
#' @param x character vector of date-time values
#'		current accepted input formats: "m/%d/%Y %I:%M:%S %p"                                (AM/PM)
#'									    "%Y-%m-%d %H:%M:%S" or "%m/%d/%Y %H:%M:%S"    (with seconds)
#'									    "Y-%m-%d %H:%M" or "%m/%d/%Y %H:%M"        (without seconds)
#' @param secs are seconds included in time values. Logical
#' @param AMPM is time format in AM, PM format. Logical
#'
#' @return striptime object in "%Y-%m-%d %H:%M:%S" or "%Y-%m-%d %H:%M" format
#'
#' @examples
#' \dontrun{
#' dayhr.fn("2021-05-10 14:22:34", secs =T , AMPM = F)
#' }
#' @export
#'
dayhr.fn=function(x,secs=T,AMPM=F){
	# return Y-M-D H:M:S or Y-M-D H:M
	if(AMPM){x=format(strptime(x, "%m/%d/%Y %I:%M:%S %p"), "%m/%d/%Y %H:%M:%S")}

	if(secs){ #if seconds included
		if(sum(is.na(strptime(x,format="%Y-%m-%d %H:%M:%S")))==length(x)){out=strptime(x,format="%m/%d/%Y %H:%M:%S")}else{
		out=strptime(x,format="%Y-%m-%d %H:%M:%S")}
		}else{ # no seconds
			if(sum(is.na(strptime(x,format="%Y-%m-%d %H:%M")))==length(x)){out=strptime(x,format="%m/%d/%Y %H:%M")}else{
			out=strptime(x,format="%Y-%m-%d %H:%M")}
		}
#  stopifnot('Time format in detection history is not correct'=sum(is.na(out))==length(x))
  return(out)
	}

