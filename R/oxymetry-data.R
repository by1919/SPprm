#' Oxymetry comparison data for 20 subjects
#'
#' A representative data from a single-center, non-randomized, prospective study conducted by
#' a medical device company to demonstrate equivalency in performance between a FDA-cleared oxymetry system and an investigational system.
#' Thirty-two healthy, non-smoking (or has refrained from smoking for 2 days) adults
#' and adolescent volunteers, ages 18 to 50 years olds, have been recruited.
#' Multiple measures of paired difference of tissue oxygen saturation levels have been obtained for each subject.
#' Due to the proprietary nature of the data, we simulated a representative data with 20 subjects
#' to mimic the original data set based on estimated model from the original data.
#'
#' @format A dataframe with two components
#' \describe{
#'   \item{Y}{ differences of repeatedly measured oxygen levels between a FDA-clearted oxymetry and an investigational oxymetry system.  }
#'   \item{A}{ individual ids (coded as 1 to 20) }
#' }
#' 
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2018) Statistical test with sample size and power calculation for paired repeated measures designs of method comparison studies.
## ' @source \url{http://www.github.com/baolinwu/SPprm}
#' @examples
#' data(poData)
#' summary(poData$Y); table(poData$A)
#' RMSt(poData$Y,poData$A)
"poData"


