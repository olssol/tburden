#' Get survival data for multi-state model
#'
#'
#' @export
#'
tb_surv_imp <- function(fit_rst, imp_m, ..., seed = NULL) {

    if (!is.null(seed)) {
        message(paste("tb_surv_imp: Random seed set to ", seed))
        old_seed <- set.seed(seed)
    }

    ## subjects
    vec_covs <- c("SUBJID", "ARM", "RANDT",
                  "PFS_DAYS", "OS_DAYS",
                  "PFS_CNSR", "OS_CNSR",
                  "time",     "status",
                  "T_Event",  "T_Premean",
                  all.vars(fit_rst$fml_surv)[-(1:2)])

    d_subs <- fit_rst$dat_imp_surv %>%
        select(any_of(vec_covs)) %>%
        distinct()

    ## impute
    n_sub <- nrow(d_subs)
    rst   <- NULL
    for (i in seq_len(n_sub)) {
        d     <- d_subs[i, ]
        f_imp <- switch(fit_rst$method,
                        msm     = tb_msm_imp_single,
                        weibull = tb_weibull_imp_single)

        cur_rst        <- f_imp(d, fit_rst, imp_m, ...)
        cur_rst$SUBJID <- d$SUBJID
        cur_rst$RANDT  <- d$RANDT
        rst            <- rbind(rst, cur_rst)
    }

    ## reset seed
    if (!is.null(seed)) {
        set.seed(old_seed)
    }

    ## note: weibull distribution imputes only one type of event note:
    ## the event is marked as death to be compatible with the msm program
    ## setting
    rst <- rst %>%
        mutate(IT_EVENT = factor(IT_Event,
                                 1:2,
                                 c("Progression", "Death")))

    ## return
    rst
}

#' Summarize imputed survival data
#'
#'
#' @export
#'
tb_summary_imp <- function(imp_surv, dat_surv, inx_imp = NULL,
                           by_var = c("ARM")) {
    dat_surv <- imp_surv %>%
        left_join(dat_surv, by = "SUBJID")

    if (!is.null(inx_imp)) {
        dat_surv <- dat_surv %>%
            filter(Imp == inx_imp)
    }

    dat_surv %>%
        group_by(!!as.name(by_var)) %>%
        summarize(Progression_Rate  = mean(is.na(IT_PFS)),
                  Progession_Mean   = mean(IT_PFS,   na.rm = T),
                  Progession_Median = median(IT_PFS, na.rm = T),
                  OS_Mean           = mean(IT_OS,    na.rm = T),
                  OS_Median         = median(IT_OS,  na.rm = T)
                  )
}
