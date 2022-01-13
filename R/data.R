#' Single Release ATLAS Practice File.
#'
#' @format A data frame with 8 columns and 9996 rows in the ATLAS file format:
#' \describe{
#' \item{V1}{Release group (ie. R1)}
#' \item{V2}{tag lot number (ie. 1,2,...)}
#' \item{V3}{tag id code}
#' \item{V4}{tag activation date/time  (yyyy-mm-dd hh:mm:ss)}
#' \item{V5}{release date/time  (yyyy-mm-dd hh:mm:ss)}
#' \item{V6}{detecton site designation}
#' \item{V7}{detection indicator (0: not detected; 1: detected; 2: detected & censored)}
#' \item{V8}{tag detection date/time  (yyyy-mm-dd hh:mm:ss) at that site}
#' ...
#' }
"single.rel"

#' Taglife ATLAS Practice File.
#'
#' @format A data frame with 1 columns and 100 rows:
#' \describe{
#' \item{tag_life_days}{length of time from activation to tag failure in days}
#' ...
#' }
"taglife.data"
