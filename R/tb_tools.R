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


#' Get survival curves cut by a time point
#'
#'
#' @export
#'
tb_surv_cut <- function(surv_f, t_dur = NULL) {
    old_f  <- rbind(c(0, 1), surv_f)

    if (!is.null(t_dur))  {
        surv_f <- rbind(old_f, c(Inf, 0))
        inx    <- which(surv_f[, 1] <= t_dur)
        inx    <- max(inx)
        surv_f_dur <- rbind(surv_f[1:inx, ],
                        cbind(t_dur, surv_f[inx, 2]))
    } else {
        surv_f_dur <- old_f
    }

    list(surv_f     = old_f,
         surv_f_dur = surv_f_dur)
}

#' Get data from All results
#'
#'
#' @export
#'
tb_extract_rst <- function(rst_all) {
    list(imp_surv      = rst_all$rst_orig$imp_surv,
         fit_msm       = rst_all$rst_orig$msm_fit$mdl_fit,
         dat_tb        = rst_all$rst_orig$params$dat_tb,
         dat_surv      = rst_all$rst_orig$params$dat_surv,
         formula_surv  = rst_all$rst_orig$params$fml_surv,
         formulat_tb   = rst_all$rst_orig$params$fml_tb,
         uti_gamma     = rst_all$rst_orig$params$uti_gamma,
         date_dbl      = rst_all$rst_orig$params$date_dbl,
         reg_tb        = rst_all$rst_orig$reg_tb,
         estimate      = rst_all$rst_orig$estimate_sub,
         params        = rst_all$rst_orig$params,
         results       = rst_all$summary)
}

#' Get AIC
#'
#'
#' @export
#'
tb_extract_aic <- function(rst_orig) {
    reg_tb <- rst_orig$reg_tb
    rst    <- NULL

    for (i in 1:length(reg_tb)) {
        for (j in names(reg_tb[[i]])) {
            rst <- rbind(rst,
                         data.frame(Imp = i,
                                    ARM = j,
                                    AIC = reg_tb[[i]][[j]]$fit_aic))
        }
    }

    rst
}


#'  Cox Regression
#'
#'
#'
#'
#' @export
#'
tb_coxph <- function(dat_surv, var_time, var_status, fml = "ARM", event = 0) {

    dat_surv$time   <- dat_surv[[var_time]]
    dat_surv$status <- dat_surv[[var_status]]

    dat_surv <- dat_surv %>%
        mutate(status = if_else(event == status, 1, 0))

    s_fml <- paste("Surv(time, status) ~", fml)
    s_fml <- as.formula(s_fml)
    fit   <- coxph(s_fml, data = dat_surv)

    fit
}
