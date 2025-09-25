#' @title Standardizes date-time input
#'
#' @description This function converts date and time into a "yyyy-mm-dd hh:mm:ss" format.
#'
#' @param x Character vector of date-time values currently accepted as input formats
#' \describe{
#'    \item{"m/d/Y I:M:S p"}{AM/PM}
#'	  \item{"Y-m-d H:M:S" or "m/d/Y H:M:S"}{with seconds}
#'	  \item{"Y-m-d H:M" or "m/d/Y H:M"}{without seconds}
#' }
#' @param secs (T|F) Seconds are included in time values
#' @param AMPM (T|F) Time format is in AM/PM format
#'
#' @return Returns a strptime object in "%Y-%m-%d %H:%M:%S" or "%Y-%m-%d %H:%M" format
#'
#' @examples
#' \dontrun{dayhr.fn("2021-05-10 14:22:34", secs=T , AMPM = F)}
#'
#' @export
#'
dayhr.fn=function(x,secs=T,AMPM=F){
	# return Y-M-D H:M:S or Y-M-D H:M
	if(AMPM){x=format(strptime(x, "%m/%d/%Y %I:%M:%S %p"), "%m/%d/%Y %H:%M:%S")}
	if((!AMPM)&(sum(is.na(strptime(x,format="%Y-%m-%d %H:%M:%S")))==length(x))){secs=F} # check if seconds included

	if(secs){ #if seconds included
		if(sum(is.na(strptime(x,format="%Y-%m-%d %H:%M:%S")))==length(x)){out=strptime(x,format="%m/%d/%Y %H:%M:%S")}else{
		out=strptime(x,format="%Y-%m-%d %H:%M:%S")}
		}else{ # no seconds
			if(sum(is.na(strptime(x,format="%Y-%m-%d %H:%M")))==length(x)){out=strptime(x,format="%m/%d/%Y %H:%M")}else{
			out=strptime(x,format="%Y-%m-%d %H:%M")}
		}
    return(out)
	}

