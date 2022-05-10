#' Simulate Survival time for all
#'
#'
#'@export
#'
tb_simu_surv <- function(n, hd_prog = 0.5, hd_death = 1, hd_cens = 1,
                         dur_enroll = 1, dur_study = 3,
                         seed = NULL) {

    ## random seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)


    t_enroll <- runif(n, min = 0, max = dur_enroll)
    t_prog   <- sapply(1:n, function(x) {
        rst <- tkt_simu_pwexp(hd_prog)
        rst[1]
    })

    t_death  <- sapply(1:n, function(x) {
        rst <- tkt_simu_pwexp(hd_death)
        rst[1]
    })

    t_cens   <- sapply(1:n, function(x) {
        rst <- tkt_simu_pwexp(hd_cens)
        rst[1]
    })

    ## summarize
    rst <- cbind(t_enroll, t_cens, t_prog, t_death)
    rst <- apply(rst, 1, function(x) {
        o_dur <- dur_study - x[1]

        if (x[3] < min(x[2], x[4], o_dur)) {
            o_prog <- x[3]
        } else {
            o_prog <- NA
        }

        if (x[4] < min(x[2], o_dur)) {
            o_death <- x[4]
        } else {
            o_death <- NA
        }

        if (is.na(o_death)) {
            o_cens <- min(x[2], o_dur)
        } else {
            o_cens <- NA
        }

        c(x, o_prog, o_death, o_cens, o_dur)
    })

    ## set seed
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst           <- t(rst)
    colnames(rst) <- c("T_Enroll", "T_Cens",  "T_Prog", "T_Death",
                       "O_Prog",   "O_Death", "O_Cens", "O_Dur")
    rownames(rst) <- NULL

    data.frame(rst)
}


#' Simulate Survival time for all based on weibull
#'
#' @param min_fu_day minimum follow up days
#'
#'@export
#'
tb_simu_surv_wb <- function(pt_cov, mdl_fit,
                            trt_effect     = 0,
                            rand_effect    = 0.1,
                            annual_dropout = 0.1,
                            min_fu_days    = 0.6 * 365.25,
                            ...,
                            seed = NULL) {
    ## random seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    mdl_fit$coefficients["ARM"] <- trt_effect
    n_sub   <- nrow(pt_cov)
    rand_ef <- rnorm(n_sub)

    ## survival time
    all_st_mean <- NULL
    all_st      <- NULL
    for (i in seq_len(n_sub)) {
        cur_d  <- pt_cov[i, ]
        cur_t  <- tb_weibull_imp(cur_d, mdl_fit,
                                 imp_m    = 1,
                                 t_censor = 0,
                                 offset   = rand_ef[i] * rand_effect)
        all_st      <- c(all_st, cur_t$cur_imp)
        all_st_mean <- c(all_st_mean, cur_t$pred_mean)
    }

    ## censoring
    lambda  <- tkt_lambda(annual_dropout, tp  = 365.25)
    if (0 == lambda) {
        all_ct <- rep(Inf, n_sub)
    } else {
        all_ct <- rexp(n_sub, lambda)
    }

    ## observed
    date_dbl   <- max(pt_cov$RANDT) + min_fu_days
    all_obs    <- apply(cbind(all_st,
                              all_ct,
                              pt_cov$RANDT,
                              date_dbl),
                        1,
                        function(x) {
                            y    <- x[1:3]
                            y[3] <- x[4] - x[3]
                            c(min(y), y[1] <= y[2] & y[1] <= y[3])
                        })

    ## append survival results
    pt_cov$T_Event   <- all_st
    pt_cov$T_Censor  <- all_ct
    pt_cov$T_Premean <- all_st_mean
    pt_cov$PFS_DAYS  <- all_obs[1, ]
    pt_cov$PFS_CNSR  <- all_obs[2, ]

    ## random effect
    pt_cov$rand_ef  <- rand_ef

    ## reset seed
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    list(pt       = pt_cov,
         date_dbl = date_dbl)
}
