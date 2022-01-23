## -----------------------------------------------------------------
##
##     FUNCTIONS RELATED TO EVALUATE THE IMPACT OF GAMMA
##           SIMULATION
##           PRESENTATION
##           OPTIMAZATION
##
## -----------------------------------------------------------------

#' Simulate survival outcomes for RCT
#'
#' @export
#'
tb_gamma_simu_surv <- function(n = 200, hd_death = c(1, 1.5), ...) {

    rst <- NULL
    for (j in 1:2) {
        simu_pt <- tb_simu_surv(n          = n,
                                hd_prog    = 0.5,
                                hd_death   = hd_death[j],
                                ...)

        simu_pt <- simu_pt %>%
            mutate(arm   = j - 1,
                   time  = if_else(!is.na(O_Death), O_Death, O_Cens),
                   event = if_else(!is.na(O_Death), 1,       0))

        rst <- rbind(rst, simu_pt)
    }

    rst

}

#' Evaluate treatment effect by AUC_SURV
#'
#' @export
#'
tb_gamma_simu_uti <- function(simu_fit, uti_gamma = 1) {
    dta_surv  <- simu_fit$dta_surv
    lst_fit   <- get_surv_f(simu_fit)
    uti_true  <- get_uti_true(dta_surv,  uti_gamma = uti_gamma)
    uti_fsurv <- get_uti_fsurv(dta_surv, lst_fit, uti_gamma = uti_gamma)

    dta_surv$Uti_True <- uti_true
    dta_surv$Uti_Est  <- uti_fsurv

    ## summary results
    dta_summary <- dta_surv %>%
        group_by(arm) %>%
        summarize(O_Dur    = mean(O_Dur),
                  Uti_True = mean(Uti_True),
                  Uti_Est  = mean(Uti_Est))

    uti_edur <- NULL
    for (j in 1:2) {
        cur_uti <- tb_uti_surv_prob_single(lst_fit[[j]],
                                           dta_summary$O_Dur[j],
                                           uti_gamma)
        uti_edur <- c(uti_edur, cur_uti)
    }

    dta_summary$Uti_Edur <- uti_edur

    ## effect
    est_uti      <- dta_summary$Uti_Est[2]  - dta_summary$Uti_Est[1]
    est_uti_true <- dta_summary$Uti_True[2] - dta_summary$Uti_True[1]


    ## return
    list(dta_surv    = dta_surv,
         dta_summary = dta_summary,
         est_effect  = c(est_uti, est_uti_true),
         pval        = simu_fit$pval)
}


#' Evaluate treatment effect by AUC_SURV by BS
#'
#' @export
#'
tb_gamma_simu_bs <- function(simu_pt, uti_gamma,
                             n_bs = 200, seed = 10000,
                             n_cores = 5, ...) {

    ## random seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    simu_fit <- tb_km_surv(simu_pt)
    simu_rst <- tb_gamma_simu_uti(simu_fit, uti_gamma = uti_gamma)

    n      <- nrow(simu_pt)
    rst_bs <- parallel::mclapply(seq_len(n_bs),
                                 function(x) {
                                     bs_pt <- simu_pt[sample(1 : n,
                                                             n,
                                                             replace = TRUE), ]
                                     bs_fit <- tb_km_surv(bs_pt)
                                     bs_rst <- tb_gamma_simu_uti(bs_fit, uti_gamma)
                                     bs_rst$est_effect
                                 }, mc.cores = n_cores)

    rst_bs <- simplify2array(rst_bs)
    rst_bs <- t(rst_bs)
    est_sd <- apply(rst_bs, 2, sd)
    est    <- simu_rst$est_effect

    pvals  <- c(tkt_pval(est[1], est_sd[1]),
                tkt_pval(est[2], est_sd[2]))

    ## random seed
    if (!is.null(seed))
        old_seed <- set.seed(old_seed)

    rst <- c("KM"    = simu_rst$pval,
             "EST"   = pvals[1],
             "TRUTH" = pvals[2])
}

#' Plot Gamma related utilities
#'
#' @export
#'
tb_gamma_plt_uti <- function(dta_surv) {
    ggplot(data = dta_surv, aes(x = Uti_True, y = Uti_Est)) +
        geom_point() +
        theme_bw() +
        facet_wrap(~arm)
}


## -----------------------------------------------------------------------------
##
##                   PRIVATE FUNCTIONS
##
## -----------------------------------------------------------------------------

#' AUC_Surv based on true data
#'
get_uti_true <- function(dta_surv, uti_gamma) {
    ## true utility
    uti_true <- apply(dta_surv[, c("T_Death", "O_Dur")], 1,
                      function(x) {
                          tb_uti_surv_obs_single(x[1],
                                                 x[2],
                                                 uti_gamma)
                      })
}

#' AUC_Surv based on true data
#'
get_uti_fsurv <- function(dta_surv, lst_fit, uti_gamma) {
    apply(dta_surv[, c("arm", "O_Dur")], 1,
          function(x) {
              tb_uti_surv_prob_single(lst_fit[[x[1] + 1]],
                                      x[2],
                                      uti_gamma)
          })
}

## get survival functions for arms 0 and 1
get_surv_f <- function(simu_fit) {
    surv_fit <- summary(simu_fit$surv_fit)
    ## arms 0:1
    lst_fit <- list()
    for (j in 1:2) {
        cur_strata   <- paste("arm=", j - 1, sep = "")
        cur_inx      <- which(cur_strata == surv_fit$strata)
        lst_fit[[j]] <- cbind(surv_fit$time[cur_inx],
                              surv_fit$surv[cur_inx])
    }

    lst_fit
}
