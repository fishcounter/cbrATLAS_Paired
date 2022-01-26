#' @title Estimates travel time harmonic mean
#'
#' @description This function estimates the travel time harmonic mean.
#'
#' @param data.in Detection history and detection times (flat format)
#' @param num.period Number of detection sites
#' @param site.names Vector of site names
#'
#' @return  This function returns a list with:
#' \describe{
#'		\item{harmonic.tt}{harmonic mean travel time}
#'		\item{tt.matrix}{individual travel time, release to first detection time for each site (days)}
#'		\item{activetime.matrix}{individual time tag active, time of tag activation to first detection for each site (days)}
#'}
#' @export
#'
mean_tt2site.fn=function(data.in,num.period,site.names){
  mean.tt = data.frame(matrix(0, nrow = num.period, ncol = 2))
  mean.tt[, 1] = site.names
  names(mean.tt) = c("Release to", "Travel Time (s.e.) num.tags")
  tt.matrix = data.frame(matrix(0, nrow = dim(data.in)[1], ncol = num.period))
  names(tt.matrix)=site.names
  active.matrix=tt.matrix

  # estimate travel time, release to 1st detection time at each site
  for (i in 1:num.period) {
    tt.matrix[,i] = difftime(dayhr.fn(data.in[, 3 + num.period + i]),
                             dayhr.fn(data.in[, 3]),units="days")
  }

    for (i in 1:num.period) {
    temp.time = as.numeric(tt.matrix[!is.na(tt.matrix[, i]), i])
    mean.tt[i, 2] = harmonic(temp.time)
  }

  # estimate travel time, tag activation to 1st detection time at each site
  for (i in 1:num.period) {
    active.matrix[,i] = difftime(dayhr.fn(data.in[, 3 + num.period + i]),
                             dayhr.fn(data.in[, 2]),units="days")
    }

  return(list(harmonic.tt=mean.tt,tt.matrix=tt.matrix,activetime.matrix=active.matrix))
}
