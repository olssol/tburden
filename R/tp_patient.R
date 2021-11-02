#' Generate history of a patient
#'
#'
#' @export
#'
tb_pt_hist <- function(time_tb, tb, time_event, t_ana,
                       gamma = c(0.2, 0.5), locf = FALSE) {

    ## remove NA pfs from event
    if (is.na(time_event[1])) {
        time_event <- time_event[-1]
        gamma      <- gamma[-1]
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
    rst_y   <- c(tb, gamma)
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
    rst_uti <- 0
    for (i in 2:nrow(pri_ana)) {
        y0  <- pri_ana[i - 1, 2] + 1
        y1  <- pri_ana[i,     2] + 1
        x0  <- pri_ana[i - 1, 1]
        x1  <- pri_ana[i,     1]

        cur_u   <- (y1 + y0) * (x1 - x0) / 2
        rst_uti <- rst_uti + cur_u
    }

    ## return
    rst <- list(pri_ana = pri_ana,
                pos_ana = pos_ana,
                t_ana   = t_ana,
                utility = rst_uti - t_ana,
                gamma   = gamma)
}


#' Get patient's data
#'
#'
#' @export
#'
tb_get_pt <- function(id, imp_surv, dat_tb, date_dbl,
                      t_ana = NULL, imp_inx = 1, ...) {

    d_surv <- imp_surv %>%
        filter(SUBJID == id &
               Imp    == imp_inx) %>%
        data.frame()

    d_tb <- dat_tb %>%
        filter(SUBJID == id &
               !is.na(PCHG)) %>%
        arrange(DAY)

    time_event <- as.numeric(d_surv[1, c("IT_PFS", "IT_OS")])
    time_tb    <- d_tb$DAY
    tb         <- d_tb$PCHG

    if (is.null(t_ana)) {
        t_ana <- as.Date(date_dbl) - as.Date(d_tb$RANDT[1])
        t_ana <- as.numeric(t_ana)
    }

    rst_hist    <- tb_pt_hist(time_tb, tb, time_event, t_ana = t_ana, ...)
    rst_hist$id <- id

    ## return
    rst_hist
}
