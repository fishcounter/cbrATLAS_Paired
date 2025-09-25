#' @title Creates a list with counts per history and corresponding matrix of unique histories
#'
#' @description This function creates a list with counts per history and a corresponding matrix of unique histories.
#'
#' @param x Detection history in wide ATLAS format without rel.group or bin columns
#'
#' @return This function returns a list with:
#' \describe{
#'    \item{count}{total counts in the same order as unique histories in hist.matrix}
#'    \item{hist.matrix}{unique detection matrix, columns indicate detection site. (0=not detected/1=detected)}
#' }
#'
#' @export
#'
thist0=function(x) {
  count=summary(factor(apply(x,1,paste,collapse="")))
  tmp=names(count)
  hist.matrix=as.data.frame(matrix(unlist(strsplit(tmp,split="")),nrow=length(tmp),byrow=T))
  return(list(count=count,hist.matrix=hist.matrix))
}
