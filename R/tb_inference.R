#' Get average utility
#'
#'
#'@export
#'
tb_estimate <- function(dat_sub, dat_tb, imp_surv, reg_tb = NULL,
                        ...) {

    if (is.null(imp_surv))
        return(NULL)

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
        left_join(data.frame(rst), by = "inx")
}

#' Summarize estimation results
#'
#'
#'@export
#'
tb_estimate_summary <- function(rst_estimate,
                                arm_control = "Chemotherapy",
                                ...) {

    if (is.null(rst_estimate))
        return(NULL)

    ## check control arm label
    stopifnot(arm_control %in% unique(rst_estimate$ARM))

    rst_value <- rst_estimate %>%
        group_by(ARM, imp) %>%
        summarize(utility      = mean(utility),
                  adj_utility  = mean(adj_utility),
                  uti_tb       = mean(uti_tb),
                  uti_event    = mean(uti_event),
                  t_ana        = mean(t_ana)) %>%
        gather(Outcome, Value, utility, adj_utility,
               uti_tb, uti_event, t_ana)

    rst_arm <- rst_value %>%
        group_by(ARM, Outcome) %>%
        summarize(Value = mean(Value))

    rst_effect <- rst_value %>%
        mutate(Value = if_else(ARM == arm_control, -Value, Value)) %>%
        ungroup() %>%
        group_by(Outcome, imp) %>%
        summarize(Value = sum(Value)) %>%
        ungroup() %>%
        group_by(Outcome) %>%
        summarize(Value = mean(Value))

    list(rst_arm    = rst_arm,
         rst_effect = rst_effect)
}

#' Overall results
#'
#' @param dat_tb tumor burden dataset
#' @param dat_surv survival dataset
#' @param fml_surv formula for survival model
#' @param fml_tb formula for tumor burden model
#' @param imp_m number of imputations for each subject
#' @param fit_tb whether fit tumor burden curve or only use the observed tumor
#'     burden
#' @param date_dbl database lock date, i.e., analysis date
#' @param uti_gamm utility gamma for progression and death
#' @param scenario label of scenarios
#' @param mdl_surv regression method for survival. msm: multi-state model;
#'     weibull: weibull regression for PFS
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
                       mdl_surv  = c("msm", "weibull"),
                       ...,
                       seed      = NULL) {

    f_abs <- function(d) {
        d %>%
            mutate(inx_bs   = inx_bs,
                   Scenario = scenario) %>%
            data.frame()
    }

    if (!is.null(seed)) {
        message(paste("tb_get_all: Random seed set to ", seed))
        old_seed <- set.seed(seed)
    }

    ## collect all parameters
    params <- c(as.list(environment()),
                list(...))

    ## remove seed for bootstrap
    params$seed <- NULL

    ## survival regression model approach
    mdl_surv <- match.arg(mdl_surv)

    ## bootstrap samples
    bs_smp   <- tb_draw_bs(dat_tb, dat_surv, inx_bs)
    dat_tb   <- bs_smp$dat_tb
    dat_surv <- bs_smp$dat_surv

    ## survival and imputation
    if ("msm" == mdl_surv) {
        ##multi-state survival data
        dat_imp_surv <- tb_msm_set(dat_surv) %>%
            mutate(time = if_else(0 == time,
                                  10,
                                  time))
        f_surv_fit <- tb_msm_fit
    } else {
        dat_imp_surv <- dat_surv
        f_surv_fit   <- tb_weibull_fit
    }

    surv_fit <- NULL
    imp_surv <- NULL
    tryCatch({
        ## survival: fit model
        surv_fit  <- f_surv_fit(dat_imp_surv, fml_surv, ...)
        ## survival imputation
        imp_surv <- tb_surv_imp(surv_fit, imp_m = imp_m, ...)
    }, error = function(e) {
        message("Error in survival imputation")
    })

    ## tb regression
    if (fit_tb) {
        reg_tb <- tb_regression(dat_tb,
                                imp_surv,
                                fml_tb    = fml_tb,
                                uti_gamma = uti_gamma,
                                ...)
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
    rst <- tb_estimate_summary(rst_est, ...)
    if (!is.null(rst)) {
        rst$rst_arm    <- f_abs(rst$rst_arm)
        rst$rst_effect <- f_abs(rst$rst_effect)
    }

    ## estimate parameters
    est_par <- tb_get_para_est(surv_fit, reg_tb)

    ## reset random seed
    if (!is.null(seed)) {
        message("tb_get_all: Random seed reset")
        set.seed(old_seed)
    }

    ## return
    if (0 != inx_bs) {
        rst_est <- NULL
    }

    list(scenario     = scenario,
         params       = params,
         surv_fit     = surv_fit,
         imp_surv     = imp_surv,
         reg_tb       = reg_tb,
         f_estimate   = "tb_get_all",
         estimate_par = est_par,
         estimate_sub = rst_est,
         estimate_arm = rst$rst_arm,
         estimate     = rst$rst_effect)
}

#' Overall results with bootstrap results
#'
#' @export
#'
tb_get_all_bs <- function(rst_orig, nbs = 100, n_cores = 5, seed = NULL) {

    if (!is.null(seed)) {
        message(paste("tb_get_all_bs: Random seed set to", seed))
        old_seed <- set.seed(seed)
    }

    all_seeds <- ceiling(abs(rnorm(nbs)) * 100000)

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
                                  cur_rst <- do.call(f_estimate,
                                                     c(params,
                                                       list(seed = all_seeds[x])))
                                  cur_rst$estimate
                              },
                              mc.cores = n_cores)

    if (!is.null(seed)) {
        message(paste("tb_get_all_bs: Random seed reset"))
        set.seed(old_seed)
    }

    ## summary
    rst <- rbind(rst_orig$estimate,
                 rbindlist(rst))

    summary <- rst %>%
        filter(0 == inx_bs) %>%
        select(-inx_bs) %>%
        left_join(rst %>%
                  filter(0 != inx_bs) %>%
                  group_by(Scenario, Outcome) %>%
                  summarize(n     = n(),
                            bs_sd = sd(Value)),
                  by = c("Scenario", "Outcome")) %>%
        mutate(LB     = Value - 1.96 * bs_sd,
               UB     = Value + 1.96 * bs_sd,
               pvalue = 2 * (1 - pnorm(abs(Value / bs_sd))))

    ## return
    list(rst_orig = rst_orig,
         summary  = summary,
         bs_rst   = rst,
         nbs      = nbs)
}
