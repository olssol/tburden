## ----------------------------------------------------------
##
##              FUNCTIONS FOR STUDY DESIGN: TUMOR BURDEN
##
## ----------------------------------------------------------

#' Polynomial tumor burden
#'
#'
#' @export
#'
tb_des_tpoly <- function(vec_t, par_t,
                         intercept = FALSE,
                         ft = NULL) {
    if (intercept) {
        rst   <- par_t[1]
        par_t <- par_t[-1]
    } else {
        rst <- 0
    }

    if (!is.null(ft))
        vec_t <- ft(vec_t)

    for (j in seq_len(length(par_t))) {
        rst <- rst + par_t[j] * vec_t^j
    }

    rst
}


#' Define Pseudo Overall Response
#'
#' @param thresh_cr Threshold for CR/PR: >=30% shrinkage from baseline
#' @param thresh_pd Threshold for PD: >= 20% increase than the minimum
#'
#' @export
#'
tb_des_pseu_or <- function(pchg,
                           vec_t     = NULL,
                           thresh_cr = -0.3,
                           thresh_pd = 0.2) {
    f_or <- function(cur_pchg) {
        inx <- which(!is.na(cur_pchg))
        res <- c_pseudo_response(cur_pchg[inx], thresh_cr, thresh_pd)
        if (-1 == res[1]) {
            rt <- NA
            rp <- NA
        } else {
            rt <- vec_t[inx][res[1]]
            rp <- cur_pchg[inx][res[1]]
        }

        c(res[2], rt, rp, res[3])
    }

    pchg <- rbind(pchg)
    if (is.null(vec_t)) {
        vec_t <- seq_len(ncol(pchg)) - 1
    }

    ## check response
    rst           <- apply(pchg, 1, f_or)
    rst           <- t(rst)
    colnames(rst) <- c("Response", "Time", "PCHG", "PCHG_Min")
    rownames(rst) <- NULL
    rst           <- data.frame(rst) %>%
        mutate(Response = factor(Response,
                                 levels = -1:1,
                                 labels = c("SD", "PD", "Response"))) %>%
        mutate(OR = if_else(Response == "Response",
                            "Yes", "No"))
    rst
}

#' Find tumor burden related design parameters
#'
#'
#' @export
#'
tb_des_tpar <- function(mean_pchg,
                        sd_pchg    = NULL,
                        vec_t      = NULL,
                        pred_vec_t = vec_t,
                        poly_o     = 3,
                        label      = "Control",
                        ...) {

    fr_1 <- function(par) {
        rst <- tb_des_tpoly(vec_t, par, ...) - mean_pchg
        sum(rst^2)
    }

    if (is.null(vec_t)) {
        vec_t <- seq_len(length(mean_pchg)) - 1
    }

    ## mean
    par_tb  <- optim(rep(0, poly_o), fr_1)$par
    pred_y  <- tb_des_tpoly(pred_vec_t, par_tb, ...)

    list(par_tb     = par_tb,
         vec_t      = c(0, vec_t),
         mean_pchg  = c(0, mean_pchg),
         pred_vec_t = c(0, pred_vec_t),
         pred_y     = c(0, pred_y),
         label      = label)
}

#' Find tumor burden related design parameters
#'
#'
#' @export
#'
tb_des_tpar_rand <- function(obs_surv, dta_surv, par_tb, par_tb_sig = 0.05,
                             par_tb_rand = c(0, 0), k = 100, ...) {

    frst <- function(pars, keep = FALSE) {
        dta_tb <- tb_des_simu_tb_giv_surv(dta_surv,
                                          par_tb,
                                          par_tb_rand = pars[-1],
                                          par_tb_sig  = exp(pars[1]))

        ## survival
        ds <- NULL
        for (tc in tb_cuts) {
            cur_ds <- tb_tb_surv(dta_tb, dta_surv,
                                 tb_cut      = tc,
                                 var_time    = "T_Event",
                                 var_status  = "event",
                                 surv_quants = NULL)

            ds <- rbind(ds, cur_ds)
        }

        ## std
        dt <- dta_tb %>%
            filter(DAY > 0) %>%
            group_by(DAY) %>%
            summarize(N         = n(),
                      PCHG_MEAN = mean(PCHG),
                      PCHG_SIG  = sd(PCHG)) %>%
            filter(N > 10)

        if (keep) {
            rst <- list(obs_surv = ds,
                        obs_tb   = dt,
                        dta_tb   = dta_tb)
        } else {
            rst <- c(ds$Days, k * mean(dt$PCHG_SIG))
        }

        rst
    }

    fval <- function(pars) {
        val <- frst(pars)
        val <- val - c(obs_surv$Days, k * par_tb_sig)
        val <- sum(val^2)
        val
    }

    dta_surv <- dta_surv %>% mutate(event = 0)
    tb_cuts  <- unique(obs_surv$TB_Cut)
    par_rst  <- optim(par = c(log(par_tb_sig), par_tb_rand),
                      fn  = fval,
                      ...)

    rst      <- frst(par_rst$par, keep = TRUE)

    list(par_tb_rand = par_rst$par[-1],
         par_tb_sig  = exp(par_rst$par[1]),
         dta_tb      = rst$dta_tb,
         obs_surv    = rst$obs_surv,
         obs_tb      = rst$obs_tb)
}

#' Simulate tumor burden for study design
#'
#' @param sub_rand Random characteristics of each subject
#' @param vec_t   Time without 0
#'
#'
#' @export
#'
tb_des_simu_tb <- function(vec_t, par_tb,
                           par_tb_rand = 0,
                           par_tb_sig  = 0.25,
                           n = 500, sub_rand = NULL,
                           label = "Control",
                           rand_cut = 0, ...) {

    f_tb <- function(w) {
        rst <- w * par_tb_rand[1]
        rst <- rst + tb_des_tpoly(vec_t, par_tb)
        rst <- rst + rnorm(nt, 0, par_tb_sig)

        ## set boundary
        inx <- which(rst < -1)
        if (length(inx) > 0)
            rst[inx] <- -1

        rst
    }

    if (is.null(sub_rand)) {
        sub_rand <- rnorm(n)
    }

    nt  <- length(vec_t)
    ni  <- length(sub_rand)
    sid <- seq_len(ni)

    rst <- sapply(sub_rand, f_tb)
    rst <- t(rst)

    ## pseudo response
    rst_or         <- tb_des_pseu_or(rst, vec_t, ...)
    rst_or$SUBJID  <- sid
    rst_or$rand_ef <- sub_rand

    ## tumor burden
    all_tb <- cbind(0, rst)
    rst    <- cbind(rep(sid, each = ncol(all_tb)),
                    rep(c(0, vec_t), ni),
                    c(t(all_tb)))
    colnames(rst) <- c("SUBJID", "DAY", "PCHG")

    rst %>%
        data.frame() %>%
        left_join(rst_or %>%
                  rename(PCHG_OR = PCHG),
                  by = "SUBJID") %>%
        mutate(ARM    = label,
               PosSub = if_else(rand_ef > rand_cut,
                                "Positive",
                                "Negative"))
}

#' Simulate tumor burden for study design
#'
#' @param sub_rand Random characteristics of each subject
#' @param vec_t   Time without 0
#'
#'
#' @export
#'
tb_des_simu_tb_giv_surv <- function(dta_surv,
                                    par_tb,
                                    par_tb_rand = c(0, 0, 0),
                                    par_tb_sig  = 0.25,
                                    label       = NULL) {

    fw <- function(w) {
        rst <- rep(w, nrnd)^seq_len(nrnd)
        sum(rst * par_tb_rand)
    }

    frand <- function(g_w) {
        rst <- rnorm(length(vec_t), g_w + mean_pchg, par_tb_sig)

        inx <- which(rst < -1)
        if (length(inx) > 0)
            rst[inx] <- -1

        c(0, rst)
    }

    if (is.null(label)) {
        label <- par_tb$label
    }

    vec_t              <- par_tb$vec_t[-1]
    mean_pchg          <- par_tb$mean_pchg[-1]
    nrnd               <- length(par_tb_rand)
    dta_surv$g_rand_ef <- sapply(dta_surv$rand_ef, fw)

    ## fix mean
    par_a0 <- NULL
    for (vt in vec_t) {
        cur_id <- dta_surv %>%
            filter(T_Event > vt)

        par_a0 <- c(par_a0, -mean(cur_id$g_rand_ef))
    }
    mean_pchg <- par_a0 + mean_pchg

    ## random tb
    rand_pchg <- sapply(dta_surv$g_rand_ef, frand)
    dta_tb    <- data.frame(SUBJID = rep(dta_surv$SUBJID,
                                         each = nrow(rand_pchg)),
                            DAY    = rep(c(0, vec_t), nrow(dta_surv)),
                            PCHG   = c(rand_pchg),
                            ARM    = label) %>%
        left_join(dta_surv %>% select(SUBJID, T_Event, rand_ef, g_rand_ef),
                  by = "SUBJID") %>%
        filter(DAY < T_Event)

    dta_tb
}

## ----------------------------------------------------------
##
##              FUNCTIONS FOR STUDY DESIGN: SURVIVAL
##
## ----------------------------------------------------------


#' Find survival related design parameters
#'
#' @param median_surv Median survival in days
#'
#' @export
#'
tb_des_surv_par <- function(median_surv, par_surv_rand = 0,
                            label = "Control", nlarge = 20000,
                            interval = c(-10, 10)) {

    fr_1 <- function(par) {
        lambda <- par + par_surv_rand * surv_rand
        smps   <- rexp(nlarge, exp(lambda))
        (median(smps) - median_surv)^2
    }


    ## mean
    l0        <- log(log(2) / median_surv)
    surv_rand <- rnorm(nlarge)
    par_surv  <- optimize(fr_1, interval = interval, tol = 1e-6)$minimum

    list(par_surv      = par_surv,
         par_surv_rand = par_surv_rand,
         median_surv   = median_surv,
         label         = label)
}

#' Find survival related design parameters
#'
#'
#' @export
#'
tb_des_surv_par_rand <- function(surv_days, surv_quants,
                                 label = "Control", nlarge = 20000) {

    fr_1 <- function(par) {
        lambda <- par[1] + par[2] * surv_rand
        smps   <- rexp(nlarge, exp(lambda))

        quants <- quantile(smps, surv_quants)
        sum((surv_days - quants)^2)
    }

    ## mean
    l0        <- log(log(1 / surv_quants[1]) / surv_days[1])
    surv_rand <- rnorm(nlarge)
    par_surv  <- optim(par = c(l0, 0), fr_1)$par

    list(par_surv      = par_surv[1],
         par_surv_rand = par_surv[2],
         surv_days     = surv_days,
         surv_quants   = surv_quants,
         label         = label)
}


#' Simulate PFS in Days
#'
#' @param sub_rand Random characteristics of each subject
#' @param vec_c    Time without 0
#'
#'
#' @export
#'
tb_des_simu_pfs <- function(par_surv, n = 500,
                            sub_rand    = NULL,
                            label       = NULL,
                            mth_to_days = 1) {

    if (is.null(sub_rand)) {
        sub_rand <- rnorm(n)
    }

    if (is.null(label)) {
        label <- par_surv$label
    }

    lambda <- par_surv$par_surv + par_surv$par_surv_rand  * sub_rand
    lambda <- exp(lambda)
    pfs    <- rexp(length(sub_rand), lambda) * mth_to_days

    data.frame(SUBJID      = seq_len(length(sub_rand)),
               rand_ef     = sub_rand,
               surv_lambda = lambda,
               T_Event     = pfs,
               ARM         = label)
}


## ----------------------------------------------------------
##
##              FUNCTIONS FOR STUDY DESIGN: SIMULATE ALL
##
## ----------------------------------------------------------

#' Simulate data for design
#'
#' @export
#'
tb_des_simu_all <- function(true_par,
                            n               = g_sample_n,
                            enroll_dur_mth  = 6,
                            annual_dropout  = 0.05,
                            min_fu_mth      = 12,
                            visit_days      = seq(9, 120, by = 9),
                            visit_max_days  = 300,
                            miss_rate       = g_tb_mis_rate,
                            arms            = c("Control", "Treatment"),
                            seed            = NULL) {


    ## set random seed
    if (!is.null(seed)) {
        old_seed <- set.seed(seed)
    }

    ## enrollment
    dta_enl <- NULL
    for (a in arms) {
        cur <- tkt_simu_enroll(ntot            = n,
                               enroll_duration = enroll_dur_mth,
                               min_fu          = min_fu_mth)

        cur$ARM <- a
        dta_enl <- rbind(dta_enl, cur)
    }

    dta_enl <- dta_enl %>%
        rename(SUBJID = sid)

    ## tumor burden
    sub_rand <- list(rnorm(n), rnorm(n))
    par_tb   <- true_par$par_tb
    dta_tb   <- NULL
    for (a in 1:2) {
        cur <- tb_des_simu_tb(vec_t       = visit_days,
                              par_tb      = par_tb[[a]]$par_tb,
                              sub_rand    = sub_rand[[a]],
                              par_tb_rand = par_tb$rand[a],
                              par_tb_sig  = par_tb$sd,
                              label       = arms[a])

        dta_tb <- rbind(dta_tb, cur)
    }

    ## survival
    par_surv <- true_par$par_surv
    dta_surv <- NULL
    for (a in 1:2) {
        cur   <- tb_des_simu_pfs(par_surv[[a]],
                                 sub_rand = sub_rand[[a]],
                                 label    = arms[1])

        dta_surv <- rbind(dta_surv, cur)
    }

    ## censoring data
    dta_cens <- NULL
    for (a in arms) {
        cur <- data.frame(ARM    = a,
                          SUBJID = seq_len(n),
                          T_Cens = tkt_rexp(n, annual_drop = annual_drop))

        dta_cens <- rbind(dta_cens, cur)
    }

    ## reset random seed
    if (!is.null(seed)) {
        set.seed(old_seed)
    }

    ## generate all data
    rst <- tb_des_simu_obs(dta_enl, dta_tb, dta_surv, dta_cens)

    rst
}
