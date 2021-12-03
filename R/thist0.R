#' @title Create a list with counts per history and vertical matrix of observed histories only
#'
#' @param x detection history in wide ATlAS format, without rel.group or bin columns
#'
#' @return summary of counts.in per history and unique detection matrix in same order
#' @export
#'
thist0=function(x) {
  count=summary(factor(apply(x,1,paste,collapse="")))
  tmp=names(count)
  hist.matrix=as.data.frame(matrix(unlist(strsplit(tmp,split="")),nrow=length(tmp),byrow=T))
  return(list(count=count,hist.matrix=hist.matrix))
}
