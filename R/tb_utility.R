#' Utility based on observed survival only
#'
#'
#'
#' @export
#'
tb_uti_surv_obs_single <- function(t_event, t_dur, uti_gamma,
                                   type = c("jump", "connect")) {

    type <- match.arg(type)

    if (t_event > t_dur) {
        rst <- switch(type,
                      jump    = 0,
                      connect = uti_gamma * t_dur^2 / t_event / 2)
    } else {
        rst <- switch(type,
                      jump    = uti_gamma * (t_dur - t_event),
                      connect = uti_gamma * (t_dur - t_event / 2))
    }

    rst
}

#' Utility based on survival function only
#'
#'
#' @export
#'
tb_uti_surv_prob_single <- function(surv_f, t_dur, uti_gamma,
                                    type = c("jump", "connect")) {
    type   <- match.arg(type)
    surv_f <- tb_surv_cut(surv_f, t_dur)$surv_f_dur
    uti    <- c_uti_surv(surv_f) * uti_gamma

    uti
}


#' Generate history of a patient
#'
#'
#' @export
#'
tb_pt_hist <- function(time_tb, tb, time_event, t_ana,
                       uti_gamma = c(0.2, 0.5),
                       locf = FALSE, ...) {

    ## remove NA pfs from event
    if (is.na(time_event[1])) {
        time_event <- time_event[-1]
        uti_gamma  <- uti_gamma[-1]
        ## os
        ind_event  <- 2
    } else {
        ## pfs os
        ind_event  <- 1:2
    }

    ## remove tb after events
    inx <- which(time_tb > min(time_event, na.rm = TRUE))
    if (0 < length(inx)) {
        time_tb <- time_tb[-inx]
        tb      <- tb[-inx]
    }

    ## tumor burden
    rst_x   <- c(time_tb, time_event)
    rst_y   <- c(tb, uti_gamma)
    rst_ind <- c(rep(0, length(tb)), ind_event)
    ntot    <- length(rst_x)

    ## locf setting
    if (locf) {
        last_t  <- time_tb[length(time_tb)]
        last_tb <- tb[length(tb)]

        rst_y[which(rst_x > last_t)] <- last_tb
    }

    ## adjust duration
    npost <- sum(rst_x >= t_ana)

    if (0 == npost) {
        pri_ana <- cbind(x = c(rst_x, t_ana),
                         y = c(rst_y, rst_y[ntot]),
                         z = c(rst_ind, 3))

        pos_ana <- NULL
    } else {
        inx <- ntot - npost
        y0  <- rst_y[inx]
        y1  <- rst_y[inx + 1]
        x0  <- rst_x[inx]
        x1  <- rst_x[inx + 1]
        y_ana <- y0 + (y1 - y0) * (t_ana - x0) / (x1 - x0)

        pri_ana <- cbind(x = c(rst_x[1:inx],   t_ana),
                         y = c(rst_y[1:inx],   y_ana),
                         z = c(rst_ind[1:inx], 3))

        pos_ana <- cbind(x = c(t_ana, rst_x[- (1:inx)]),
                         y = c(y_ana, rst_y[- (1:inx)]),
                         z = c(3,     rst_ind[- (1:inx)]))
    }

    ## utility
    uti_tb    <- 0
    uti_event <- 0
    for (i in 2:nrow(pri_ana)) {
        y0  <- pri_ana[i - 1, 2] + 1
        y1  <- pri_ana[i,     2] + 1
        x0  <- pri_ana[i - 1, 1]
        x1  <- pri_ana[i,     1]

        cur_u   <- (y1 + y0) * (x1 - x0) / 2

        if (0 == pri_ana[i - 1, 3]) {
            uti_tb  <- uti_tb + cur_u
        } else {
            uti_event <- uti_event + cur_u
        }
    }
    rst_uti <- uti_tb + uti_event

    if (min(time_event) < t_ana) {
        dur_event <- t_ana - min(time_event)
    } else {
        dur_event <- 0
    }
    dur_tb <- t_ana - dur_event

    ## return
    rst <- list(pri_ana     = pri_ana,
                pos_ana     = pos_ana,
                t_ana       = t_ana,
                dur_tb      = dur_tb,
                dur_event   = dur_event,
                utility     = rst_uti - t_ana,
                adj_utility = (rst_uti - t_ana) / t_ana,
                uti_tb      = uti_tb    - dur_tb,
                uti_event   = uti_event - dur_event,
                uti_gamma   = uti_gamma)
}
