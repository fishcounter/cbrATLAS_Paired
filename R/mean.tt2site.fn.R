#' @title Estimate the travel time harmonic mean, release to each detection site
#'
#' @param data.in detection history and detection times (flat format)
#' @param num.period number of detection sites
#' @param site.names vector of site names
#'
#' @return harmonic.tt: harmonic mean of travel time (days); tt.matrix: time of release to 1st detection time for each site (days); activetime.matrix: time of tag activation to 1st detection for each site (days)
#'
#' @export
#'
mean.tt2site.fn=function(data.in,num.period,site.names){
  #
  # INPUT
  # data.in: detection history and detection times (flat format)
  # num.period: number of detection sites
  # site.names: vector of site names
  # OUTPUT
  # harmonic.tt: harmonic mean of travel time (days)
  # tt.matrix: time of release to 1st detection time for each site (days)
  # activetime.matrix: time of tag activation to 1st detection for each site (days)

  mean.tt = data.frame(matrix(0, nrow = num.period, ncol = 2))
  mean.tt[, 1] = site.names
  names(mean.tt) = c("Release to", "Travel Time (s.e.) num.tags")
  tt.matrix = data.frame(matrix(0, nrow = dim(data.in)[1], ncol = num.period))
  names(tt.matrix)=site.names
  active.matrix=tt.matrix

  # estimate travel time, release to 1st detection time at each site
  for (i in 1:num.period) {
    tt.matrix[, i] = strptime(data.in[, 3 + num.period + i],
                              "%m/%d/%Y %H:%M") - strptime(data.in[, 3], "%m/%d/%Y %H:%M")
    units(tt.matrix[, i]) = "days"
  }
  for (i in 1:num.period) {
    temp.time = as.numeric(tt.matrix[!is.na(tt.matrix[, i]),
                                     i])
    mean.tt[i, 2] = harmonic(temp.time)
  }
  # estimate travel time, tag activation to 1st detection time at each site
  for (i in 1:num.period) {
    active.matrix[, i] = strptime(data.in[, 3 + num.period + i],
                                  "%m/%d/%Y %H:%M") - strptime(data.in[, 2], "%m/%d/%Y %H:%M")
    units(active.matrix[, i]) = "days"
  }

  return(list(harmonic.tt=mean.tt,tt.matrix=tt.matrix,activetime.matrix=active.matrix))
}

harmonic=function(t.in){
  # estimate mean travel times between 2 sites.  Harmonic tt is less influence by large outliers.
  # INPUT
  # vector of times (days)
  # OUTPUT
  # vector of harmonic travel time, standard error, and number of times used in estimate
  w=1/t.in
  t.out=length(w)/sum(w)
  t.var=(var(w)/length(w))*(t.out^4)
  t.se=sqrt(t.var)
  return(paste(c(sprintf("%.4f",t.out),"  (",sprintf("%.4f",t.se),")      ",length(t.in)),collapse=""))
}
