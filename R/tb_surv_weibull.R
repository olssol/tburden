#' Fit PFS Weibull
#'
#'
#' @export
#'
tb_weibull_fit <- function(dat_surv,
                           fml_surv,
                           type_surv = c("PFS", "OS"),
                           ind_event = 0,
                           by_arm    = TRUE,
                           ...) {

    type_surv       <- match.arg(type_surv)
    var_status      <- paste(type_surv, "_", "CNSR", sep = "")
    var_time        <- paste(type_surv, "_", "DAYS", sep = "")
    dat_surv$time   <- dat_surv[[var_time]]
    dat_surv$status <- dat_surv[[var_status]] == ind_event

    arms  <- unique(dat_surv$ARM)
    fml   <- as.formula(paste("Surv(time, status)",
                              fml_surv))

    if (by_arm) {
        rst   <- list()
        for (a in arms) {
            cur_d <- dat_surv %>%
                filter(ARM == a)

            cur_rst  <- survreg(fml, data = cur_d, dist = 'weibull')
            rst[[a]] <- cur_rst
        }
    } else {
        rst <- survreg(fml, data = dat_surv, dist = 'weibull')
    }

    list(dat_imp_surv = dat_surv,
         mdl_fit      = rst,
         fml_surv     = fml,
         type_surv    = type_surv,
         by_arm       = by_arm,
         method       = "weibull")
}

#' Impute survival for single patient based on weibull regression
#'
#' @export
#'
tb_weibull_imp_single <- function(d, fit_rst, imp_m, ...) {

    mdl_fit <- fit_rst$mdl_fit

    ## event time observed
    if (1 == d$status) {
        ## 3: pfs
        cur_rst <- cbind(seq_len(imp_m), d$time, 2)
    } else {
        ## 3: pfs
        if (fit_rst$by_arm) {
            mf <- mdl_fit[[d$ARM]]
        } else {
            mf <- mdl_fit
        }

        cur_imp <- tb_weibull_imp(d, mf, imp_m, d$time)
        cur_rst <- cbind(seq_len(imp_m), cur_imp, 2)
    }

    colnames(cur_rst) <- c("Imp", "IT_Time", "IT_Event")
    rownames(cur_rst) <- NULL

    ## return
    rst <- data.frame(cur_rst)

    ## placeholder for plots etc
    rst$IT_PFS <- NA
    rst$IT_OS  <- rst$IT_Time

    rst
}

#' Impute weibull survial
#'
#'
#' @export
#'
tb_weibull_imp <- function(d, fit_rst, imp_m, t_censor = 0, offset = 0) {
    pred_mean <- predict(fit_rst, newdata = d) * exp(offset)
    cur_imp   <- NULL
    while (length(cur_imp) < imp_m) {
        cur_weibull <- rweibull(10, 1, scale = 1)
        cur_t       <- pred_mean * cur_weibull^fit_rst$scale

        ## longer than censored time
        inx <- which(cur_t > t_censor)
        if (0 < length(inx))
            cur_imp <- c(cur_imp, cur_t[inx])
    }

    cur_imp <- cur_imp[1 : imp_m]
}
