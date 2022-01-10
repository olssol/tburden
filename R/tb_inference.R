#' Get average utility
#'
#'
#'@export
#'
tb_estimate <- function(dat_sub, dat_tb, imp_surv, reg_tb = NULL,
                        ...) {
    imp_m <- max(imp_surv$Imp)
    nsub  <- nrow(dat_sub)

    rst   <- NULL
    for (imp in seq_len(imp_m)) {
        cur_imp    <- NULL
        tot_tana   <- 0
        for (i in seq_len(nsub)) {
            cur_uti <- tb_get_pt(id       = dat_sub[i, "SUBJID"],
                                 imp_surv = imp_surv,
                                 dat_tb   = dat_tb,
                                 imp_inx  = imp,
                                 reg_tb   = reg_tb,
                                 ...)

            tot_tana <- tot_tana + cur_uti$t_ana
            cur_rst  <- c(i,
                          cur_uti$utility,
                          cur_uti$adj_utility,
                          cur_uti$uti_tb,
                          cur_uti$uti_event,
                          cur_uti$t_ana)

            cur_imp <- rbind(cur_imp, cur_rst)
        }

        rst <- rbind(rst, cbind(imp, cur_imp))
    }

    colnames(rst) <- c("imp",     "inx",
                       "utility", "adj_utility",
                       "uti_tb",  "uti_event",
                       "t_ana")
    dat_sub %>%
        mutate(inx = 1:n()) %>%
        left_join(data.frame(rst))
}

#' Summarize estimation results
#'
#'
#'@export
#'
tb_estimate_summary <- function(rst_estimate, arm_control = "Chemotherapy") {

    ## check control arm label
    stopifnot(arm_control %in% unique(rst_estimate$ARM))

    rst_estimate %>%
        group_by(ARM, imp) %>%
        summarize(utility      = mean(utility),
                  adj_utility  = mean(adj_utility),
                  uti_tb       = mean(uti_tb),
                  uti_event    = mean(uti_event),
                  t_ana        = mean(t_ana)) %>%
        gather(Outcome, Value, utility, adj_utility,
               uti_tb, uti_event, t_ana) %>%
        mutate(Value = if_else(ARM == arm_control, -Value, Value)) %>%
        ungroup() %>%
        group_by(Outcome, imp) %>%
        summarize(Value = sum(Value)) %>%
        ungroup() %>%
        group_by(Outcome) %>%
        summarize(Value = mean(Value))
}

#' Overall results
#'
#' @export
#'
tb_get_all <- function(dat_tb, dat_surv,
                       inx_bs    = 0,
                       fml_surv  = "~BASE+AGE+SEX+STRATA1+P1TERTL",
                       fml_tb    = "~AGE+SEX+STRATA1+P1TERTL",
                       imp_m     = 5,
                       fit_tb    = TRUE,
                       date_dbl  = "2020-03-01",
                       uti_gamma = c(0.2, 0.5),
                       scenario  = "scenario",
                       ...) {

    params <- c(as.list(environment()),
                list(...))

    ## bootstrap samples
    if (0 != inx_bs) {
        d_subjid <- dat_tb %>%
            select(SUBJID) %>%
            distinct()

        d_subjid <- d_subjid[sample(nrow(d_subjid), replace = TRUE), ,
                             drop = FALSE]

        dat_tb   <- d_subjid %>%
            left_join(dat_tb)
        dat_surv <- d_subjid %>%
            left_join(dat_surv)
    }

    ## survival: multi-state survival data
    msm_surv <- tb_msm_set(dat_surv) %>%
        mutate(time = if_else(0 == time,
                              10,
                              time))

    ## survival: fit model
    msm_fit  <- tb_msm_fit(msm_surv, fml_surv)

    if (is.null(msm_fit)) {
        warning("Failure in survival model fitting")
        return(NULL)
    }

    ## survival imputation
    imp_surv <- tb_msm_imp(msm_fit, imp_m = imp_m)

    ## tb regression
    if (fit_tb) {
        reg_tb <- tb_regression(dat_tb, imp_surv,
                                fml_tb    = fml_tb,
                                uti_gamma = uti_gamma, ...)
    } else {
        reg_tb <- NULL
    }

    ## estimate
    dat_sub <- dat_tb %>%
        select(SUBJID, ARM) %>%
        distinct()

    rst_est <- tb_estimate(dat_sub   = dat_sub,
                           imp_surv  = imp_surv,
                           dat_tb    = dat_tb,
                           reg_tb    = reg_tb,
                           date_dbl  = date_dbl,
                           uti_gamma = uti_gamma,
                           ...)

    ## estimate summary
    rst <- tb_estimate_summary(rst_est) %>%
        mutate(inx_bs   = inx_bs,
               Scenario = scenario) %>%
        data.frame()

    ## bootstrap
    if (0 != inx_bs) {
        rst_est <- NULL
    }

    ## return
    list(scenario     = scenario,
         params       = params,
         msm_surv     = msm_surv,
         msm_fit      = msm_fit,
         imp_surv     = imp_surv,
         reg_tb       = reg_tb,
         f_estimate   = "tb_get_all",
         estimate_sub = rst_est,
         estimate     = rst)
}

#' Overall results with bootstrap results
#'
#' @export
#'
tb_get_all_bs <- function(rst_orig, nbs = 100, seed = 1234, n_cores = 5) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    if (0 == nbs) {
        rst <- list(rst_orig = rst_orig,
                    nbs      = nbs)
        return(rst)
    }

    f_estimate <- get(rst_orig$f_estimate)
    params     <- rst_orig$params

    rst <- parallel::mclapply(seq_len(nbs),
                              function(x) {
                                  cat("--Rep ", x, "\n")
                                  params$inx_bs <- x
                                  rst           <- do.call(f_estimate, params)

                                  rst$estimate
                              }, mc.cores = n_cores)

    if (!is.null(seed))
        set.seed(old_seed)

    ## summary
    rst     <- rbind(rst_orig$estimate,
                     rbindlist(rst))

    summary <- rst %>%
        filter(0 == inx_bs) %>%
        select(-inx_bs) %>%
        left_join(rst %>%
                  filter(0 != inx_bs) %>%
                  group_by(Scenario, Outcome) %>%
                  summarize(bs_sd = sd(Value))) %>%
        mutate(LB     = Value - 1.96 * bs_sd,
               UB     = Value + 1.96 * bs_sd,
               pvalue = 2 * (1 - pnorm(abs(Value / bs_sd))))

    ## return
    list(rst_orig = rst_orig,
         summary  = summary,
         bs_rst   = rst,
         nbs      = nbs)
}
