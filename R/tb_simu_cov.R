## ------------------------------------------------------------------
##
##                SIMULATE COVARIATES
##
## ------------------------------------------------------------------

#' Simulate covariates
#'
#' Simulate covariates based on existing real studies
#'
#' @param dta_es           Existing study data
#' @param v_covs           Vector of covariates to be sampled
#' @param n                Sample size for two arms
#' @param rnd_ratio_trt    Randomization ratio vs. control (=1)
#' @param label_arm        Column name for Arm in the result
#' @param label_id         Column name for patient id in the existing data
#' @param replace          Sample with or without replacement
#'
#' @export
#'
tb_simu_cov_es <- function(dta_es, v_covs, n = 500,
                           label_arm     = "ARM",
                           label_id      = "SUBJID",
                           label_randt   = "RANDT",
                           rnd_ratio_trt = 1, replace = TRUE,
                           ...,
                           seed = NULL) {

    if (!is.null(seed)) {
        message(paste("tb_simu_cov_es: Random seed set to ", seed))
        old_seed <- set.seed(seed)
    }


    ## remove duplications
    dta_es <- dta_es %>%
        select(c(label_id, label_randt, v_covs)) %>%
        na.omit() %>%
        distinct()

    ## covariates
    n_pt <- nrow(dta_es)
    smps <- sample(1:n_pt, n, replace = replace)
    rst  <- dta_es[smps, v_covs]

    ## arm assignment
    smp_arm          <- rbinom(n, 1, prob = rnd_ratio_trt / (rnd_ratio_trt + 1))
    rst[[label_arm]] <- smp_arm

    ## randomization time
    if (!is.null(label_randt)) {
        stopifnot(label_randt %in% names(dta_es))
        smps_dt <- sample(dta_es[[label_randt]],
                          n,
                          replace = replace)

        rst[[label_randt]] <- smps_dt
    }

    ## id
    rst$SUBJID <- as.character(seq_len(nrow(rst)))

    ## reset random seed
    if (!is.null(seed)) {
        set.seed(old_seed)
    }

    ## return
    rownames(rst) <- NULL
    rst
}
