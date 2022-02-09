## Simulation

#' Simulate studies
#'
#'
#' @export
#'
tb_simu_all <- function(dat_tb,
                        n,
                        surv_fit,
                        tb_fit,
                        trt_effect_surv  = 0,
                        trt_effect_tb    = 0,
                        rand_effect_surv = 0,
                        rand_effect_tb   = c("(Intercept)"     = 0,
                                             "poly(reg_t, 2, raw = TRUE)1" = 0,
                                             "poly(reg_t, 2, raw = TRUE)2" = 0),

                        ..., seed = NULL) {

    ## random seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    simu_pt_cov  <- tb_simu_cov_es(dta_es = dat_tb, n = n, ...)
    simu_pt_surv <- tb_simu_surv_wb(simu_pt_cov,
                                    surv_fit,
                                    trt_effect  = trt_effect_surv,
                                    rand_effect = rand_effect_surv,
                                    ...)

    simu_pt_tb <- tb_simu_tb(simu_pt_surv$pt,
                             tb_fit,
                             trt_effect  = trt_effect_tb,
                             rand_effect = rand_effect_tb,
                             ...)

    ## reset random seed
    if (!is.null(seed)) {
        set.seed(old_seed)
    }

    ## result
    dat_tb   <- simu_pt_tb %>%
        mutate(ARM = as.character(ARM)) %>%
        select(- reg_t)

    dat_surv <- simu_pt_surv$pt %>%
        mutate(ARM = as.character(ARM))

    list(dat_tb   = dat_tb,
         dat_surv = dat_surv,
         date_dbl = simu_pt_surv$date_dbl)

}

#' Summarize simulated patients
#'
#'
#' @export
#'
tb_simu_present <- function(simu_pt) {
    plt_surv <- plot_km(simu_pt, "PFS_DAYS", "PFS_CNSR", event = 1)
    plt_tb   <- ggplot(data = simu_pt %>%
                           filter(TB_mis == 0) %>%
                           mutate(ARM = factor(ARM)),
                       aes(x = DAY, y = PCHG)) +
        geom_line(aes(group = SUBJID, col = ARM)) +
        theme_bw()

    tb_mis <- simu_pt %>%
        group_by(DAY) %>%
        summarize(m = mean(TB_mis))

    list(plt_surv = plt_surv,
         plt_tb   = plt_tb,
         tb_mis   = tb_mis)
}
