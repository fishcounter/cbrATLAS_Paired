#' @title Creates a list with counts per history and a vertical matrix of observed histories only.
#'
#' @param x Detection history in wide ATLAS format without rel.group or bin columns.
#'
#' @return This function returns a summary of counts.in per history and unique detection matrix in the same order.
#' @export
#'
thist0=function(x) {
  count=summary(factor(apply(x,1,paste,collapse="")))
  tmp=names(count)
  hist.matrix=as.data.frame(matrix(unlist(strsplit(tmp,split="")),nrow=length(tmp),byrow=T))
  return(list(count=count,hist.matrix=hist.matrix))
}
