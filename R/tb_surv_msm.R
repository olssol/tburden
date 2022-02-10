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
        left_join(data.frame(rst), by = "inx") %>%
        select(-inx) %>%
        mutate(time = t_stop - t_start)
}

#' Fit Survival MSM
#'
#'
#' @export
#'
tb_msm_fit <- function(msm_surv, fml_surv,
                       ind_event = 0,
                       by_arm    = TRUE,
                       ...) {

    f_trans <- function(md) {
        lst_a <- rep(list(NULL), length(trans))
        for (i in trans) {
            cur_d <- md %>%
                filter(trans == i)

            if (0 == nrow(cur_d)) {
                warning(paste("Not enough survival outcomes for ",
                              a, " trans ", i))
                next
            }

            cur_rst    <- flexsurvreg(fml, data = cur_d, dist = "exp")
            lst_a[[i]] <- cur_rst
        }

        lst_a
    }

    fml   <- as.formula(paste("Surv(time, status)",
                              fml_surv))

    arms  <- unique(msm_surv$ARM)
    trans <- unique(msm_surv$trans)
    msm_surv <- msm_surv %>%
        mutate(status = if_else(status == ind_event, 1, 0))

    if (by_arm) {
        rst  <- list()
        for (a in arms) {
            md <- msm_surv %>%
                filter(ARM == a)
            rst[[a]] <- f_trans(md)
        }
    } else {
        rst <- f_trans(msm_surv)
    }

    list(dat_imp_surv = msm_surv,
         mdl_fit      = rst,
         fml_surv     = fml,
         by_arm       = by_arm,
         method       = "msm")
}

#' Impute survival for single patient
#'
#' @export
#'
tb_msm_imp_single <- function(d, fit_rst, imp_m, ind_event = 0, ...) {

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
    mdl_fit   <- fit_rst$mdl_fit
    if (fit_rst$by_arm) {
        fit_msm   <- mdl_fit[[d$ARM]]
    } else {
        fit_msm   <- mdl_fit
    }

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
    if (ind_event == d$OS_CNSR &
        ind_event == d$PFS_CNSR) {
        cur_rst <- f_1()
    } else if (ind_event != d$OS_CNSR &
               ind_event != d$PFS_CNSR) {
        cur_rst <- f_2()
    } else if (ind_event == d$PFS_CNSR &
               ind_event != d$OS_CNSR) {
        cur_rst <- f_3()
    } else if (ind_event != d$PFS_CNSR &
               ind_event == d$OS_CNSR) {
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
    rst
}
