#' Get patient data
#'
#' Get turmor burden, survival and utility data of a specific patient
#'
#' @param date_dbl Date of database lock
#'
#' @export
#'
tb_get_pt <- function(id, imp_surv, dat_tb,
                      date_dbl = NULL, t_ana = NULL,
                      imp_inx = 1, reg_tb = NULL,
                      ...) {

    d_surv <- imp_surv %>%
        filter(SUBJID == id &
               Imp    == imp_inx) %>%
        data.frame()

    d_tb <- dat_tb %>%
        filter(SUBJID == id &
               !is.na(PCHG)) %>%
        arrange(DAY)

    time_tb <- d_tb$DAY
    tb      <- d_tb$PCHG
    tb_obs  <- cbind(time_tb, tb)

    ## analysis time
    if (is.null(t_ana)) {
        stopifnot(!is.null(date_dbl))
        t_ana <- as.Date(date_dbl) - as.Date(d_tb$RANDT[1])
        t_ana <- as.numeric(t_ana)
    }


    ## event time
    time_event <- as.numeric(d_surv[1, c("IT_PFS", "IT_OS")])

    ## tumor burden
    if (!is.null(reg_tb)) {
        cur_pred <- tb_predict(id,
                               reg_tb[[imp_inx]][[d_tb[1, "ARM"]]],
                               ...)

        time_tb  <- cur_pred$time_tb
        tb       <- cur_pred$tb
    }

    ## get patient history
    rst_hist        <- tb_pt_hist(time_tb, tb, time_event, t_ana = t_ana, ...)
    rst_hist$id     <- id
    rst_hist$tb_obs <- tb_obs

    ## return
    rst_hist
}
