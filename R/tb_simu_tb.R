#' Simulate Tumor Burden
#'
#' @param miss_rate missing probability at scaled time 1 (i.e., event time)
#'
#'@export
#'
tb_simu_tb <- function(pt_surv, mdl_fit,
                       visit_interval = 63,
                       trt_effect     = 0,
                       tb_sig         = 2,
                       rand_effect    = c("(Intercept)"     = 0,
                                            "poly(reg_t, 2, raw = TRUE)1" = 0,
                                          "poly(reg_t, 2, raw = TRUE)2" = 0),
                       miss_rate      = 0.5,
                       ...,
                       seed = NULL) {

    f_mis <- function(rt) {
        miss_rate * rt
    }

    ## random seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    n_sub <- nrow(pt_surv)
    rst   <- NULL

    ## visit time
    for (i in seq_len(n_sub)) {
        cur_d <- pt_surv[i, ]
        days  <- seq(0, cur_d$T_Event, by = visit_interval)
        reg_t <- days
        if (mdl_fit$scale_time_by_surv) {
            reg_t <- reg_t / cur_d$T_Event
        }

        cur_d       <- cur_d[rep(1, length(reg_t)), ]
        cur_d$DAY   <- days
        cur_d$reg_t <- reg_t
        rst  <- rbind(rst, cur_d)
    }

    ## turmor burden
    fml     <- as.formula(paste(mdl_fit$fml_tb, "+", mdl_fit$fml_tb_t))
    des_mat <- model.matrix(fml, data = rst)

    coef_fix        <- fixef(mdl_fit$fit_reg)
    coef_fix["ARM"] <- trt_effect
    coef_l          <- names(rand_effect)

    simu_tb <- apply(cbind(rst$rand_ef,
                           des_mat), 1,
                     function(x) {
                         cur_coef <- coef_fix
                         cur_coef[coef_l] <- x[1] * rand_effect +
                             cur_coef[coef_l]
                         cur_mean <- sum(x[-1] * cur_coef)
                         rnorm(1, cur_mean, sd = tb_sig)
                     })

    rst$TB <- simu_tb

    rst <- rst %>%
        arrange(SUBJID, DAY)

    ## percent change
    cur_id   <- -1
    all_pchg <- NULL
    for (i in seq_len(nrow(rst))) {
        if (rst[i, "SUBJID"] != cur_id) {
            cur_base <- rst[i, "TB"]
            cur_id   <- rst[i, "SUBJID"]
        }

        all_pchg <- c(all_pchg,
                      rst[i, "TB"] / cur_base - 1)
    }
    rst$PCHG <- all_pchg

    ## missingness
    all_mis <- sapply(rst$reg_t, function(x) {
        pmis <- f_mis(x)
        rbinom(1, 1, pmis)
    })
    rst$TB_mis <- all_mis

    ## set seed
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}
