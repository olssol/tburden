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

#' Utility from the original results
#'
#'
#' @export
#'
tb_uti_rst <- function(rst_orig, ...) {

    dat_tb  <- rst_orig$params$dat_tb
    dat_sub <- dat_tb %>%
        select(SUBJID, ARM) %>%
        distinct()

    rst <- NULL
    for (w in seq(0.05, 0.95, by = 0.05)) {
        cur_g   <- w / (1 - w)
        cur_est <- tb_estimate(dat_sub   = dat_sub,
                               imp_surv  = rst_orig$imp_surv,
                               dat_tb    = dat_tb,
                               reg_tb    = rst_orig$reg_tb,
                               date_dbl  = rst_orig$params$date_dbl,
                               uti_gamma = c(cur_g, cur_g),
                               ...)

        cur_summary <- cur_est$estimate %>%
            group_by(SUBJID, ARM) %>%
            summarise(uti_tb    = mean(uti_tb),
                      uti_event = mean(uti_event)) %>%
            group_by(ARM) %>%
            summarise(n             = n(),
                      var_uti_tb    = var(uti_tb) / n,
                      var_uti_event = var(uti_event) / n) %>%
            data.frame()

        cur_var <- sum(cur_summary[1, 3:4]) +
            sum(cur_summary[2, 3:4])

        rst <- rbind(rst, c(w, cur_g, cur_var))
    }

    rst
}

#' Draw bootstrap samples
#'
#'
#' @export
#'
tb_draw_bs <- function(dat_tb, dat_surv, inx_bs = 0) {
    if (0 != inx_bs) {
        d_subjid <- dat_tb %>%
            select(SUBJID) %>%
            distinct()

        d_subjid <- d_subjid[sample(nrow(d_subjid),
                                    replace = TRUE), ,
                             drop = FALSE] %>%
            mutate(bs_id = as.character(1:n()))

        dat_tb   <- d_subjid %>%
            left_join(dat_tb, by = "SUBJID") %>%
            select(- SUBJID) %>%
            rename(SUBJID = bs_id)
        dat_surv <- d_subjid %>%
            left_join(dat_surv, by = "SUBJID") %>%
            select(- SUBJID) %>%
            rename(SUBJID = bs_id)
    }

    list(dat_tb   = dat_tb,
         dat_surv = dat_surv)
}
