#' The 'tburden' package.
#'
#' @docType package
#' @name    tburden-package
#' @aliases tburden
#' @useDynLib tburden, .registration = TRUE
#'
#' @import methods
#' @import stats
#' @import ggplot2
#'
#' @importFrom grDevices colors
#' @importFrom graphics axis box legend lines par plot points text arrows grid rect
#' @importFrom parallel detectCores
#' @importFrom utils as.roman
#' @importFrom dplyr %>% group_by_ group_by summarize mutate count mutate_if
#'     rename filter select arrange ungroup n distinct left_join if_else rowwise
#' @importFrom tidyr gather
#' @importFrom data.table rbindlist
#' @importFrom flexsurv flexsurvreg
#' @importFrom survival Surv survfit coxph survdiff
#' @importFrom survminer ggsurvplot
#' @importFrom nlme lme
#' @importFrom JM jointModel
#' @importFrom lme4 lmer
#'
#' @description Tumor burden based treatment effect evaluation
#'
#' @references
#'
NULL
