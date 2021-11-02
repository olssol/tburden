#'  Plot survival cures
#'
#'
#'
#'
#' @export
#'
plot_km <- function(dat_surv, var_time, var_status, lab_y = "PFS",
                    by_var = c("ARM"), event = 0,
                    pval  = FALSE, lim_x = 0, lim_y = c(0, 1)) {

    dat_surv$time   <- dat_surv[[var_time]]
    dat_surv$status <- dat_surv[[var_status]]

    dat_surv <- dat_surv %>%
        mutate(status = if_else(event == status, 1, 0))

    s_fml <- paste("Surv(time, status) ~",
                   paste(by_var, collapse = "+"))

    s_fml            <- as.formula(s_fml)
    fit              <- survfit(as.formula(s_fml), data = dat_surv)
    fit$call$formula <- s_fml

    rst <- ggsurvplot(fit,
                      data = dat_surv,
                      pval = pval)$plot
    rst <- rst +
        labs(y = lab_y) +
        ylim(lim_y[1], lim_y[2])

    if (lim_x > 0) {
        rst <- rst + coord_cartesian(xlim = c(0, lim_x))
    }

    rst
}

#' Draw bootstrap samples
#'
#'
#' @export
#'
tb_bs_draw <- function(dat_tb, dat_surv, seed = NULL) {

    d_subjid <- dat_tb %>%
        select(SUBJID) %>%
        distinct()

    d_subjid <- d_subjid[sample(nrow(d_subjid), replace = TRUE), ,
                         drop = FALSE]

    dat_tb   <- d_subjid %>%
        left_join(dat_tb)
    dat_surv <- d_subjid %>%
        left_join(dat_surv)

    list(dat_tb   = dat_tb,
         dat_surv = dat_surv)
}
