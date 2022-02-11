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
    simu_pt  <- simu_pt$dat_tb %>%
        mutate(ARM = factor(ARM,
                            levels = 0:1,
                            labels = c("Control", "Treatment")))

    my_theme <- theme_bw() +
        theme(strip.background = element_blank())

    plt_surv <- plot_km(simu_pt, "PFS_DAYS", "PFS_CNSR", event = 1)
    plt_tb   <- ggplot(data = simu_pt %>%
                           filter(TB_mis == 0),
                       aes(x = DAY, y = PCHG)) +
        geom_line(aes(group = SUBJID)) +
        labs(y = "Percent Change from Baseline") +
        facet_wrap(~ARM, ncol = 1) +
        my_theme

    tb_mis <- simu_pt %>%
        group_by(ARM, DAY) %>%
        summarize(m = mean(TB_mis),
                  n = n()) %>%
        filter(DAY > 0)

    plt_mis <- ggplot(data = tb_mis,
                      aes(x = DAY, y = m)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = n),
                  vjust = 1.6,
                  size = 3) +
        labs(y = "Missing Rate") +
        facet_wrap(~ARM, ncol = 1) +
        my_theme

    list(plt_surv = plt_surv,
         plt_tb   = plt_tb,
         plt_mis  = plt_mis)
}

#' Summarize simulated results
#'
#'
#' @export
#'
tb_simu_summary <- function(nrep = 100, prefix = "./Results") {
    rst_value <- list()
    rst_par   <- list()
    rst_true  <- NULL

    for (i in 1:nrep) {
        fi  <- paste(prefix, "/simu_rst_",
                     i, ".Rdata", sep = "")

        if (file.exists(fi)) {
            load(fi)
        } else {
            print(paste(fi, "missing"))
            next
        }

        rst_value <- c(rst_value,
                       list(rst_all$summary %>%
                            select("Outcome", "Value", "pvalue")))

        rst_par   <- c(rst_par,
                       list(rst_all$rst_orig$estimate_par %>%
                            mutate(inx = 1:n())))

        if (is.null(rst_true)) {
            rst_true <- true_par %>%
                mutate(inx = 1:n())
        }
    }

    rst_value <- rbindlist(rst_value) %>%
        mutate(rej = pvalue < 0.05) %>%
        group_by(Outcome) %>%
        summarize(Mean      = mean(Value),
                  SD        = sd(Value),
                  Rejection = mean(rej))

    rst_par <- rbindlist(rst_par) %>%
        left_join(rst_true %>%
                  rename(Truth = Value) %>%
                  select(inx, Truth),
                  by = "inx") %>%
        mutate(bias = Value - Truth,
               mse  = (Value - Truth)^2) %>%
        group_by(Type, Para) %>%
        summarize(estimate = mean(Value),
                  truth    = mean(Truth),
                  bias     = mean(bias),
                  mse      = mean(mse))

    list(rst_value, rst_par)
}
