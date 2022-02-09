#' Get survival censoring status
#'
#'
#' @export
#'
tb_set_surv <- function(dat_surv) {

    nd  <- nrow(dat_surv)
    rst <- NULL
    for (i in seq_len(nd)) {
        d       <- dat_surv[i, , drop = TRUE]
        cur_rst <- NULL
        if (0 == d$OS_CNSR &
            0 == d$PFS_CNSR) {
            delta    <- 1
            t_censor <- NA
            t_death  <- d$OS_DAYS

            if ("DEATH" == d$PFS_EVENT) {
                t_pfs <- NA
            } else {
                t_pfs <- d$PFS_DAYS
            }
        } else if (0 != d$OS_CNSR &
                   0 != d$PFS_CNSR) {
            delta    <- 0
            t_censor <- d$OS_DAYS
            t_pfs    <- NA
            t_death  <- NA
        } else if (0 == d$PFS_CNSR &
                   0 != d$OS_CNSR) {
            delta    <- 2
            t_censor <- d$OS_DAYS
            t_pfs    <- d$PFS_DAYS
            t_death  <- NA
        }

        cur_rst <- c(delta, t_censor, t_pfs, t_death)
        rst     <- rbind(rst, cur_rst)
    }

    colnames(rst) <- c("Delta", "T_Censor", "T_PFS", "T_Death")
    rownames(rst) <- NULL
    cbind(dat_surv, data.frame(rst))
}



#' Impute survival
#'
#'
#'@export
#'
tb_imp_surv <- function(dat_surv, seed = NULL, ext = 1.2) {

    f_imp0 <- function(d) {
        cur_d <- dat_delta1 %>%
            filter(ARM == d$ARM &
                   min_days > d$T_Censor) %>%
            data.frame()

        ncur <- nrow(cur_d)
        if (0 == ncur) {
            print(d$T_Censor)
            rst <- c(NA,
                     round(d$T_Censor * ext))
        } else {
            smp <- sample(1:ncur, 1)
            rst <- c(cur_d[smp, "T_PFS"],
                     cur_d[smp, "T_Death"])
        }

        rst
    }

    f_imp2 <- function(d) {
        cur_d <- dat_lambda %>%
            filter(ARM == d$ARM) %>%
            data.frame()

        lambda <- cur_d$lambda[1]
        smp    <- rexp(1, lambda)

        c(d$T_PFS, d$T_Censor + round(smp))
    }

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## delta_1
    dat_delta1 <- dat_surv %>%
        filter(Delta == 1) %>%
        rowwise() %>%
        mutate(min_days = min(T_Death, T_PFS, na.rm = T))

    ## delta_1 lambda
    dat_lambda <- dat_delta1 %>%
        mutate(days = T_Death - T_PFS) %>%
        group_by(ARM) %>%
        summarize(lambda = 1 / mean(days, na.rm = TRUE))

    nd  <- nrow(dat_surv)
    rst <- NULL
    for (i in seq_len(nd)) {
        d       <- dat_surv[i, , drop = TRUE]
        cur_rst <- NULL
        if (1 == d$Delta) {
            cur_rst <- c(d$T_PFS, d$T_Death)
        } else if (0 == d$Delta) {
            cur_rst <- f_imp0(d)
        } else {
            cur_rst <- f_imp2(d)
        }
        rst <- rbind(rst, cur_rst)
    }

    if (!is.null(seed))
        set.seed(old_seed)

    colnames(rst) <- c("IT_PFS", "IT_OS")
    rownames(rst) <- NULL
    cbind(dat_surv, data.frame(rst))

}



#' Fit KM survival curve
#'
#'
#'  @export
#'
tb_km_surv <- function(dta_surv, formula = "Surv(time, event) ~ arm") {

    surv_fit  <- survfit(as.formula(formula), data = dta_surv)
    surv_test <- survdiff(as.formula(formula), data = dta_surv)
    pval      <- 1 - pchisq(surv_test$chisq,
                            length(surv_test$n) - 1)

    list(dta_surv  = dta_surv,
         surv_fit  = surv_fit,
         surv_test = surv_test,
         pval      = pval)
}
