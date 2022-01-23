#' Get survival censoring status
#'
#'
#'@export
#'
tb_set_surv <- function(dat_surv) {

    nd  <- nrow(dat_surv)
    rst <- NULL
    for (i in seq_len(nd)) {
        d       <- dat_surv[i, , drop = TRUE]
        cur_rst <- NULL
        if (0 == d$OS_CNSR &
            0 == d$PFS_CNSR) {
            delta    <- 1
            t_censor <- NA
            t_death  <- d$OS_DAYS

            if ("DEATH" == d$PFS_EVENT) {
                t_pfs <- NA
            } else {
                t_pfs <- d$PFS_DAYS
            }
        } else if (0 != d$OS_CNSR &
                   0 != d$PFS_CNSR) {
            delta    <- 0
            t_censor <- d$OS_DAYS
            t_pfs    <- NA
            t_death  <- NA
        } else if (0 == d$PFS_CNSR &
                   0 != d$OS_CNSR) {
            delta    <- 2
            t_censor <- d$OS_DAYS
            t_pfs    <- d$PFS_DAYS
            t_death  <- NA
        }

        cur_rst <- c(delta, t_censor, t_pfs, t_death)
        rst     <- rbind(rst, cur_rst)
    }

    colnames(rst) <- c("Delta", "T_Censor", "T_PFS", "T_Death")
    rownames(rst) <- NULL
    cbind(dat_surv, data.frame(rst))
}



#' Impute survival
#'
#'
#'@export
#'
tb_imp_surv <- function(dat_surv, seed = NULL, ext = 1.2) {

    f_imp0 <- function(d) {
        cur_d <- dat_delta1 %>%
            filter(ARM == d$ARM &
                   min_days > d$T_Censor) %>%
            data.frame()

        ncur <- nrow(cur_d)
        if (0 == ncur) {
            print(d$T_Censor)
            rst <- c(NA,
                     round(d$T_Censor * ext))
        } else {
            smp <- sample(1:ncur, 1)
            rst <- c(cur_d[smp, "T_PFS"],
                     cur_d[smp, "T_Death"])
        }

        rst
    }

    f_imp2 <- function(d) {
        cur_d <- dat_lambda %>%
            filter(ARM == d$ARM) %>%
            data.frame()

        lambda <- cur_d$lambda[1]
        smp    <- rexp(1, lambda)

        c(d$T_PFS, d$T_Censor + round(smp))
    }

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## delta_1
    dat_delta1 <- dat_surv %>%
        filter(Delta == 1) %>%
        rowwise() %>%
        mutate(min_days = min(T_Death, T_PFS, na.rm = T))

    ## delta_1 lambda
    dat_lambda <- dat_delta1 %>%
        mutate(days = T_Death - T_PFS) %>%
        group_by(ARM) %>%
        summarize(lambda = 1 / mean(days, na.rm = TRUE))

    nd  <- nrow(dat_surv)
    rst <- NULL
    for (i in seq_len(nd)) {
        d       <- dat_surv[i, , drop = TRUE]
        cur_rst <- NULL
        if (1 == d$Delta) {
            cur_rst <- c(d$T_PFS, d$T_Death)
        } else if (0 == d$Delta) {
            cur_rst <- f_imp0(d)
        } else {
            cur_rst <- f_imp2(d)
        }
        rst <- rbind(rst, cur_rst)
    }

    if (!is.null(seed))
        set.seed(old_seed)

    colnames(rst) <- c("IT_PFS", "IT_OS")
    rownames(rst) <- NULL
    cbind(dat_surv, data.frame(rst))

}


#' Get survival data for multi-state model
#'
#'
#' @export
#'
tb_msm_set <- function(dat_surv) {

    ## no progression, death
    f_1 <- function(d) {
        rbind(c(1, 1, 2, 0, d$PFS_DAYS, 0),
              c(2, 1, 3, 0, d$OS_DAYS,  1))
    }

    ## progression, death
    f_2 <- function(d) {
        rbind(c(1, 1, 2, 0,          d$PFS_DAYS, 1),
              c(2, 1, 3, 0,          d$PFS_DAYS, 0),
              c(3, 2, 3, d$PFS_DAYS, d$OS_DAYS,  1))
    }

    ## no progression, no death
    f_3 <- function(d) {
        rbind(c(1, 1, 2, 0, d$PFS_DAYS, 0),
              c(2, 1, 3, 0, d$OS_DAYS,  0))
    }

    ## progression, no death
    f_4 <- function(d) {
        rbind(c(1, 1, 2, 0,          d$PFS_DAYS, 1),
              c(2, 1, 3, 0,          d$PFS_DAYS, 0),
              c(3, 2, 3, d$PFS_DAYS, d$OS_DAYS,  0))
    }


    nd  <- nrow(dat_surv)
    rst <- NULL
    for (i in seq_len(nd)) {
        d       <- dat_surv[i, , drop = TRUE]
        cur_rst <- NULL
        if (0 == d$OS_CNSR &
            0 == d$PFS_CNSR) {
            if (d$PFS_DAYS == d$OS_DAYS) {
                cur_rst <- f_1(d)
            } else {
                cur_rst <- f_2(d)
            }
        } else if (0 != d$OS_CNSR &
                   0 != d$PFS_CNSR) {
            cur_rst <- f_3(d)
        } else if (0 == d$PFS_CNSR &
                   0 != d$OS_CNSR) {
            cur_rst <- f_4(d)
        } else if (0 != d$PFS_CNSR &
                   0 == d$OS_CNSR) {
            cur_rst <- f_1(d)
        }

        rst <- rbind(rst, cbind(i, cur_rst))
    }

    ## observed status = 1 or censored status = 0
    colnames(rst) <- c("inx", "trans", "from", "to",
                       "t_start", "t_stop", "status")

    dat_surv %>%
        mutate(inx = 1:n()) %>%
        left_join(data.frame(rst)) %>%
        select(-inx) %>%
        mutate(time = t_stop - t_start)
}

#' Fit Survival MSM
#'
#'
#' @export
#'
tb_msm_fit <- function(msm_surv, fml_surv, v_censor = 0) {

    fml   <- as.formula(paste("Surv(time, status)",
                              fml_surv))
    arms  <- unique(msm_surv$ARM)
    trans <- unique(msm_surv$trans)

    rst   <- list()
    for (a in arms) {
        lst_a <- rep(list(NULL), length(trans))

        for (i in trans) {
            cur_d <- msm_surv %>%
                filter(ARM   == a &
                       trans == i) %>%
                mutate(status = if_else(status == v_censor, 0, 1))

            if (0 == nrow(cur_d)) {
                warning(paste("Not enough survival outcomes for ",
                              a, " trans ", i))
                next
            }

            cur_rst    <- flexsurvreg(fml, data = cur_d, dist = "exp")
            lst_a[[i]] <- cur_rst
        }

        rst[[a]] <- lst_a
    }

    list(msm_surv = msm_surv,
         mdl_fit  = rst,
         fml_surv = fml)
}

#' Impute survival for single patient
#'
#' @export
#'
tb_msm_imp_single <- function(d, mdl_fit, imp_m) {

    ## progression, death
    f_1 <- function() {
        if (d$PFS_DAYS == d$OS_DAYS) {
            pfs <- NA
        } else {
            pfs <- d$PFS_DAYS
        }

        cbind(rep(pfs,        imp_m),
              rep(d$OS_DAYS,  imp_m))
    }

    ## no progression, no death
    f_2 <- function() {
        imp_pfs    <- d$PFS_DAYS + rexp(imp_m, pred_mean[1])
        imp_os     <- d$OS_DAYS  + rexp(imp_m, pred_mean[2])
        imp_pfs_os <- imp_pfs    + rexp(imp_m, pred_mean[3])

        ## censoring
        imp_pfs[which(imp_pfs >= imp_os)] <- NA
        imp_os <- apply(cbind(imp_os, imp_pfs_os), 1, min)

        cbind(imp_pfs, imp_os)
    }

    ## progression, no death
    f_3 <- function() {
        imp_os  <- d$OS_DAYS + rexp(imp_m, pred_mean[3])

        cbind(rep(d$PFS_DAYS, imp_m),
              imp_os)
    }

    ## no progression, death
    f_4 <- function() {
        imp_pfs <- d$PFS_DAYS + rexp(imp_m, pred_mean[1])

        ## censored by death
        imp_pfs[which(imp_pfs >= d$OS_DAYS)] <- NA

        cbind(imp_pfs,
              rep(d$OS_DAYS, imp_m))
    }

    ## mdl_fit result
    fit_msm   <- mdl_fit[[d$ARM]]
    pred_mean <- NULL
    for (trans in 1:3) {
        if (is.null(fit_msm[[trans]]))
            stop("Fitting multistate model failed")

        cur_mean  <- predict(fit_msm[[trans]],
                             newdata = d,
                             type    = "response")
        pred_mean <- c(pred_mean, 1 / as.numeric(cur_mean))
    }


    ## imputation
    if (0 == d$OS_CNSR &
        0 == d$PFS_CNSR) {
        cur_rst <- f_1()
    } else if (0 != d$OS_CNSR &
               0 != d$PFS_CNSR) {
        cur_rst <- f_2()
    } else if (0 == d$PFS_CNSR &
               0 != d$OS_CNSR) {
        cur_rst <- f_3()
    } else if (0 != d$PFS_CNSR &
               0 == d$OS_CNSR) {
        cur_rst <- f_4()
    }

    stopifnot(all(!is.nan(cur_rst)))

    ## Earliest event
    cur_event <- cur_rst

    cur_event[is.na(cur_event)] <- Inf
    cur_event <- apply(cur_event, 1, function(x) {
        time  <- min(x)
        event <- if_else(time == x[1], 1, 2)
        c(time, event)
    })

    ## combine results
    cur_rst <- cbind(seq_len(imp_m), cur_rst, t(cur_event))

    colnames(cur_rst) <- c("Imp", "IT_PFS", "IT_OS", "IT_Time", "IT_Event")
    rownames(cur_rst) <- NULL

    ## return
    rst <- data.frame(cur_rst)
    rst$SUBJID <- d$SUBJID
    rst$RANDT  <- d$RANDT
    rst$IT_Event <- factor(rst$IT_Event,
                           1:2,
                           c("Progression", "Death"))

    rst
}


#' Get survival data for multi-state model
#'
#'
#' @export
#'
tb_msm_imp <- function(fit_rst, ..., seed = 10000) {
    if (!is.null(seed)) {
        old_seed <- set.seed(seed)
    }

    msm_surv <- fit_rst$msm_surv

    ## subjects
    vec_covs <- c("SUBJID", "ARM", "RANDT",
                  "PFS_DAYS", "OS_DAYS",
                  "PFS_CNSR", "OS_CNSR",
                  all.vars(fit_rst$fml_surv)[-(1:2)])

    d_subs <- msm_surv %>%
        select(any_of(vec_covs)) %>%
        distinct()

    ## impute
    n_sub <- nrow(d_subs)
    rst   <- NULL
    for (i in seq_len(n_sub)) {
        cur_rst <- tb_msm_imp_single(d_subs[i, ], fit_rst$mdl_fit, ...)
        rst     <- rbind(rst, cur_rst)
    }

    ## reset seed
    if (!is.null(seed)) {
        set.seed(old_seed)
    }

    ## return
    rst
}


#' Fit KM survival curve
#'
#'
#'  @export
#'
tb_km_surv <- function(dta_surv, formula = "Surv(time, event) ~ arm") {

    surv_fit  <- survfit(as.formula(formula), data = dta_surv)
    surv_test <- survdiff(as.formula(formula), data = dta_surv)
    pval      <- 1 - pchisq(surv_test$chisq,
                            length(surv_test$n) - 1)

    list(dta_surv  = dta_surv,
         surv_fit  = surv_fit,
         surv_test = surv_test,
         pval      = pval)
}

#' Summarize imputed survival data
#'
#'
#' @export
#'
tb_summary_imp <- function(imp_surv, dat_surv, inx_imp = NULL,
                           by_var = c("ARM")) {
    dat_surv <- imp_surv %>%
        left_join(dat_surv)

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


## -----------------------------------------------------------------
##
##                  PRESENTATION
##
## -----------------------------------------------------------------

#' Plot survival curve with area under the curve
#'
#' @export
#'
tb_plt_surv <- function(surv_f, t_dur = NULL,
                        type = c("rmf", "rmst", "none"),
                        y_lim = c(0, 1), x_lim = NULL) {

    type     <- match.arg(type)
    surv_dur <- tb_surv_cut(surv_f, t_dur)$surv_f_dur
    surv_km  <- tb_surv_cut(surv_f, x_lim)$surv_f_dur

    ## survival curves
    rst    <- ggplot(data = data.frame(Time = surv_km[, 1],
                                       Y    = surv_km[, 2]),
                     aes(x = Time, y = Y)) +
        labs(x = "Time", y = "Survival Probability") +
        ylim(y_lim) +
        theme_bw() +
        geom_step()

    if (is.null(x_lim)) {
        rst <- rst + xlim(c(0, x_lim))
    }

    if (type == "none")
        return(rst)

    ## ploygon
    if (!is.null(t_dur)) {
        y_dur    <- switch(type,
                           rmst = rbind(c(surv_dur[nrow(surv_dur), 1],
                                          0),
                                        c(0, 0)),
                           rmf  = c(surv_dur[nrow(surv_dur), 1], 1))

        surv_poly <- NULL
        for (i in 1:(nrow(surv_dur) - 1)) {
            surv_poly <- rbind(surv_poly,
                               surv_dur[i, ],
                               c(surv_dur[i + 1, 1], surv_dur[i, 2]))
        }

        surv_poly <- rbind(surv_poly, y_dur)
        rst <- rst + geom_polygon(data = data.frame(x = surv_poly[, 1],
                                                    y = surv_poly[, 2]),
                                  aes(x = x, y = y),
                                  alpha = 0.2)
    }

    rst
}

#' Plot patients by enrollment
#'
#' @export
#'
tb_plt_onstudy <- function(t_enroll, t_time, event, t_dur,
                           add_auc = FALSE, auc_k = 1.2,
                           add_lab = TRUE, size_lab = 8, hjust_lab = -1,
                           h = 0.4) {

    lab_e <- c("Censored", "Event", "Enrolled")
    dat <- data.frame(t_enroll = t_enroll,
                      time     = t_time,
                      event    = event) %>%
        arrange(t_enroll) %>%
        mutate(y     = row_number(),
               event = factor(event, 0:2, lab_e))

    rst <- ggplot(data = dat, aes(x = time, y = y)) +
        geom_point(aes(pch = event, color = event)) +
        geom_vline(xintercept = t_dur, lty = 2) +
        geom_vline(xintercept = 0, lty = 2) +
        geom_text(aes(x = 0, y = 5, label = "Study Started"),
                  angle = 90, vjust = -0.5) +
        geom_text(aes(x = t_dur, y = 5, label = "Study Finished"),
                  angle = 90, vjust = 1) +
        labs(xlim = c(-0.1, t_dur * 1.05), lty = 2) +
        theme_bw() +
        theme(
            axis.line    = element_blank(),
            axis.text    = element_blank(),
            axis.ticks   = element_blank(),
            axis.title   = element_blank(),
            panel.grid   = element_blank(),
            legend.title = element_blank(),
            legend.position = "bottom")

    ## geom_point(data = data.frame(x     = dat$t_enroll,
    ##                              y     = dat$y,
    ##                              event = factor(2, 0:2, lab_e)),
    ##            aes(x     = x,
    ##                y     = y,
    ##                pch   = event,
    ##                color = event)) +

    ## add line
    for (i in 1:nrow(dat)) {
        pt    <- dat[i, "time"]
        t_enr <- dat[i, "t_enroll"]
        pe    <- "Event" == dat[i, "event"]
        pd    <- data.frame(x = c(t_enr, pt),
                            y = c(i, i))
        rst <- rst +
            geom_line(data = pd, aes(x = x, y = y))

        if (add_auc) {
            if (pe) {
                rst <- rst +
                    geom_line(data = data.frame(x = c(pt, pt, t_dur),
                                               y = c(i, i + h, i + h)),
                              aes(x = x, y = y))
                pt_imp <- pt
            } else {
                pt_imp <- runif(1, pt, pt + auc_k * (t_dur - pt))
                pt_imp <- min(pt_imp, t_dur)
                rst    <- rst +
                    geom_line(data = data.frame(
                                  x = c(pt, pt_imp, pt_imp, t_dur),
                                  y = c(i,  i, i + h, i + h)),
                              aes(x = x, y = y), lty = 2)
            }

            rst <- rst +
                geom_polygon(data = data.frame(
                                 x = c(pt_imp, pt_imp, t_dur, t_dur),
                                 y = c(i, i + h, i + h, i)),
                             aes(x = x, y = y),
                             fill = "gray30",
                             alpha = 0.2)
        }
    }

    ## add pt label
    if (add_lab) {
        dat_lab <- data.frame(x    = dat$t_enroll,
                              y    = dat$y,
                              labs = paste("P", dat$y, sep = ""))

        rst <- rst +
            geom_label(data = dat_lab,
                       aes(x = x, y = y, label = labs),
                       hjust = hjust_lab,
                       size  = size_lab)
    }

    ## return
    rst
}
