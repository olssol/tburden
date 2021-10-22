#' Get average utiltity
#'
#'
#'@export
#'
tb_estimate <- function(dat_sub, imp_surv, ...) {

    imp_m <- max(imp_surv$Imp)
    nsub  <- nrow(dat_sub)
    rst   <- NULL

    for (i in seq_len(nsub)) {
        for (imp in seq_len(imp_m)) {
            cur_uti <- tb_get_pt(id       = dat_sub[i, "SUBJID"],
                                 imp_surv = imp_surv,
                                 imp_inx  = imp, ...)
            cur_rst <- c(i,
                         imp,
                         cur_uti$utility,
                         cur_uti$auc,
                         cur_uti$auc / cur_uti$t_ana)
            rst <- rbind(rst, cur_rst)
        }
    }

    colnames(rst) <- c("inx", "imp", "utility", "auc", "adj_auc")
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
    rst_estimate %>%
        group_by(ARM, imp) %>%
        summarize(utility = mean(utility),
                  auc     = mean(auc),
                  adj_auc = mean(adj_auc)) %>%
        gather(Outcome, Value, utility, auc, adj_auc) %>%
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
                       inx_bs   = 0,
                       fml_surv = "~BASE+AGE+SEX+STRATA1+P1TERTL",
                       imp_m    = 5,
                       lst_par  = list(Calendar =
                                           list(date_dbl = "2020-03-01",
                                                gamma    = c(0.2, 0.5))),
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

    ## multistate survival data
    msm_surv <- tb_msm_set(dat_surv) %>%
        mutate(time = max(time, 10))

    ## fit imputation model
    msm_fit <- tb_msm_fit(msm_surv, fml_surv)

    ## imputation
    imp_surv <- tb_msm_imp(msm_fit, imp_m = imp_m)

    ## estimate
    dat_sub <- dat_tb %>%
        select(SUBJID, ARM) %>%
        distinct()

    rst <- NULL
    for (i in seq_len(length(lst_par))) {
        cur_par   <- lst_par[[i]]
        cur_label <- names(lst_par)[i]

        cur_rst   <- do.call(tb_estimate,
                             c(list(dat_sub  = dat_sub,
                                    imp_surv = imp_surv,
                                    dat_tb   = dat_tb),
                               cur_par,
                               list(...)))

        cur_rst   <- tb_estimate_summary(cur_rst) %>%
            mutate(inx_bs   = inx_bs,
                   Scenario = cur_label) %>%
            data.frame()

        rst <- rbind(rst, cur_rst)
    }

    ## return
    list(params   = params,
         msm_surv = msm_surv,
         msm_fit  = msm_fit,
         imp_surv = imp_surv,
         estimate = rst)
}

#' Overall results with bootstrap results
#'
#' @export
#'
tb_get_all_bs <- function(rst_orig, nbs = 100, seed = 1234, n_cores = 5) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    rst <- parallel::mclapply(seq_len(nbs),
                              function(x) {
                                  cat("--Rep ", x, "\n")
                                  params        <- rst_orig$params
                                  params$inx_bs <- x
                                  rst <- do.call(tb_get_all, params)

                                  rst$estimate
                              }, mc.cores = n_cores)

    if (!is.null(seed))
        set.seed(old_seed)

    ## summary
    rst     <- rbind(rst_orig$estimate, rbindlist(rst))
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
