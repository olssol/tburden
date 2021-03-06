## -----------------------------------------------------------------------
##
##               FUNCTIONS RELATED TO TUMOR BURDEN
##                                 REGRESSION
##                                 PREDICTION
##                                 BEST RESPONSE
##
## -----------------------------------------------------------------------

#' Combine survival and tumor burden data
#'
#' Prepare data set for longitudinal tumor burden regression
#'
#' @param add_surv Whether include uti_gamma at the event time as a tumor burden
#'     measurement
#' @param f_conv Function for converting tumor burden
#' @param scale_time_by_surv whether conditioning on survival time and scale
#'     DAYS to [0,1]
#'
#' @export
#'
tb_prepare_tb_dta <- function(dat_tb, imp_surv, imp_inx = 1,
                              uti_gamma = c("Progression" = 0.2,
                                            "Death"       = 0.5),
                              tb_reg_pow_conv    = 1,
                              tb_reg_add_surv    = FALSE,
                              scale_time_by_surv = TRUE,
                              ...) {

    imp_surv <- imp_surv %>%
        filter(Imp == imp_inx) %>%
        select(SUBJID, IT_Time, IT_Event)

    rst <- dat_tb %>%
        left_join(imp_surv, by = "SUBJID") %>%
        mutate(reg_tb    = (1 + PCHG) * BASE,
               reg_gamma = (uti_gamma[IT_Event] + 1) * BASE)

    if (scale_time_by_surv) {
        rst$reg_t <- rst$DAY / rst$IT_Time
    } else {
        rst$reg_t <- rst$DAY
    }

    ## last point corresponding to event
    if (tb_reg_add_surv) {
        last_point <- rst %>%
            group_by(SUBJID) %>%
            slice(n = 1) %>%
            mutate(reg_tb = reg_gamma)

        if (scale_time_by_surv) {
            last_point$reg_t <- 1
        } else {
            last_point$reg_t <- last_point$IT_Time
        }

        rst <- rbind(rst, last_point)
    }

    ## convert tumor burden data
    rst$reg_tb <- rst$reg_tb ^ tb_reg_pow_conv

    ## return
    rst %>%
        arrange(SUBJID, reg_t)

}

#' Tumor burden regression
#'
#' Conduct tumor burden regression
#'
#'
#' @param tb_reg_mixed Mixed model or not
#'
#' @export
#'
tb_regression_single <- function(tb_surv, fml_tb = NULL,
                                 tb_reg_poly_order = 2,
                                 tb_reg_pow_conv   = 1,
                                 tb_reg_poly_raw   = FALSE,
                                 tb_reg_mixed      = TRUE,
                                 ...) {

    ## time
    tb_reg_poly_coefs <- poly(tb_surv$reg_t,
                              tb_reg_poly_order,
                              raw = tb_reg_poly_raw)
    tb_reg_poly_coefs <- attr(tb_reg_poly_coefs, "coefs")

    ## event time
    id_time <- tb_surv %>%
        select(SUBJID, IT_Time) %>%
        distinct()

    ## design matrix
    fml_1   <- paste("poly(reg_t, ",
                     tb_reg_poly_order,
                     ", raw = ",
                     tb_reg_poly_raw,
                     ")")

    if (tb_reg_mixed) {
        fml_2 <- paste("+(", fml_1, "| SUBJID)")
    } else {
        fml_2 <- NULL
    }

    des_mat <- NULL
    if (!is.null(fml_tb)) {
        fml     <- paste("reg_tb", fml_tb, "+", fml_1, fml_2)

        vars    <- all.vars(as.formula(fml_tb))
        tmp_tb  <- tb_surv %>%
            select(c("SUBJID", vars)) %>%
            distinct()

        tmp_tb_comp <- tmp_tb %>%
            na.omit()

        if (nrow(tmp_tb_comp) < nrow(tmp_tb)) {
            warning(paste(nrow(tmp_tb) - nrow(tmp_tb_comp),
                          "patients with missing covariates were",
                          "excluded in tb_regression_single"))
        }

        des_mat <- model.matrix(as.formula(fml_tb),
                                data = tmp_tb_comp)
        row.names(des_mat) <- tmp_tb_comp$SUBJID
    } else {
        fml <- paste("reg_tb ~", fml_1, fml_2)
    }

    ## regression
    if (tb_reg_mixed) {
        fit_reg <- lmer(as.formula(fml), data = tb_surv)
        fit_rst <- coefficients(fit_reg)$SUBJID
    } else {
        fit_reg <- lm(as.formula(fml), data = tb_surv)
        fit_rst <- coefficients(fit_reg)
    }

    ## return
    list(tb_reg_poly_coefs  = tb_reg_poly_coefs,
         tb_reg_poly_order  = tb_reg_poly_order,
         tb_reg_poly_raw    = tb_reg_poly_raw,
         tb_reg_pow_conv    = tb_reg_pow_conv,
         tb_reg_mixed       = tb_reg_mixed,
         id_time            = id_time,
         des_mat            = des_mat,
         fml_tb             = fml_tb,
         fml_tb_t           = fml_1,
         fml_tb_rand        = fml_2,
         fit_rst            = fit_rst,
         fit_reg            = fit_reg,
         fit_aic            = AIC(fit_reg)
         )
}


#' Tumor burden regression
#'
#' Conduct tumor burden regression
#'
#'
#' @export
#'
tb_regression <- function(dat_tb, imp_surv,
                          scale_time_by_surv = TRUE,
                          by_arm             = TRUE,
                          ...) {

    imp_m <- max(imp_surv$Imp)
    arms  <- unique(dat_tb$ARM)

    rst   <- list()
    for (imp in seq_len(imp_m)) {
        cur_dta <- tb_prepare_tb_dta(dat_tb, imp_surv,
                                     imp_inx            = imp,
                                     scale_time_by_surv = scale_time_by_surv,
                                     ...)
        if (by_arm) {
            cur_rst <- list()
            for (a in arms) {
                cur_arm <- cur_dta %>%
                    filter(ARM == a)
                cur_rst[[a]]  <- tb_regression_single(cur_arm, ...)
                cur_rst[[a]]$scale_time_by_surv <- scale_time_by_surv
            }
        } else {
            cur_rst                    <- tb_regression_single(cur_dta, ...)
            cur_rst$scale_time_by_surv <- scale_time_by_surv
        }

        rst[[imp]] <- cur_rst
    }

    ## return rst
    rst
}

#' Tumor burden prediction
#'
#' Predict tumor burden curve
#'
#'
#' @export
#'
tb_predict <- function(id, tb_reg_fit, by_days = 7, ...) {
    pow_conv  <- tb_reg_fit$tb_reg_pow_conv
    reg_mixed <- tb_reg_fit$tb_reg_mixed
    des_mat   <- tb_reg_fit$des_mat[id, ]
    id_time   <- tb_reg_fit$id_time %>%
        filter(SUBJID == id)
    id_time   <- id_time$IT_Time[1]

    ## coefficients
    fit_rst <- tb_reg_fit$fit_rst
    if (reg_mixed) {
        fit_rst  <- fit_rst[id, ]
    }

    ## time points
    pred_t <- seq(0, id_time, by = by_days)
    if (tb_reg_fit$scale_time_by_surv) {
        pred_t_scale <- pred_t / id_time
    } else {
        pred_t_scale <- pred_t
    }

    if (tb_reg_fit$tb_reg_poly_raw) {
        poly_t  <- poly(pred_t_scale,
                        tb_reg_fit$tb_reg_poly_order,
                        raw = TRUE)
    } else {
        poly_t  <- poly(pred_t_scale,
                        tb_reg_fit$tb_reg_poly_order,
                        coefs = tb_reg_fit$tb_reg_poly_coefs)
    }

    ## predicted y
    pred_y <- sum(des_mat * fit_rst[1 : length(des_mat)])
    par_t   <- fit_rst[- (1 : length(des_mat))]
    pred_y  <- pred_y + apply(poly_t, 1, function(x) sum(x * par_t))
    pred_y  <- pred_y ^ (1 / pow_conv)
    pchg    <- pred_y / pred_y[1] - 1

    ## return
    list(time_tb = pred_t,
         pred_y  = pred_y,
         tb      = pchg)
}

#' Define Pseudo Overall Response
#'
#' @param thresh_cr Threshold for CR/PR: >=30% shrinkage from baseline
#' @param thresh_pd Threshold for PD: >= 20% increase than the minimum
#'
#' @export
#'
tb_pseu_or <- function(pchg,
                       vec_t     = NULL,
                       thresh_cr = -0.3,
                       thresh_pd = 0.2) {

    f_vt <- function(r, vi) {
        if (!is.na(r)) {
            rst <- vi[r]
        } else {
            rst <- NA
        }

        rst
    }

    f_or <- function(cur_pchg) {
        inx <- which(!is.na(cur_pchg))
        if (0 == length(inx))
            return(rep(NA, 8))

        rst <- c_pseudo_response(cur_pchg[inx], thresh_cr, thresh_pd)

        ## translate time
        vi     <- vec_t[inx]
        rst[1] <- f_vt(rst[1], vi)
        rst[3] <- f_vt(rst[3], vi)
        rst[5] <- f_vt(rst[5], vi)

        rst
    }

    pchg <- rbind(pchg)
    if (is.null(vec_t)) {
        vec_t <- seq_len(ncol(pchg)) - 1
    }

    ## check response
    rst           <- apply(pchg, 1, f_or)
    rst           <- t(rst)
    colnames(rst) <- c("T_OR",  "PCHG_OR",
                       "T_PD",  "PCHG_PD",
                       "T_Min", "PCHG_Min")

    rownames(rst) <- NULL
    rst <- data.frame(rst) %>%
        mutate(Response = if_else(is.na(T_OR), "No", "Yes"))
    rst
}


#' Get improved tumor burden
#'
#' @param perc_improve Percentage improvement
#'
#' @export
#'
tb_tb_impr <- function(vec, perc_improve = 0.1,
                       type = c("ratio", "absolute")) {

    type <- match.arg(type)

    if ("ratio" == type) {
        rst <- vec * (1 - perc_improve) - perc_improve
    } else {
        rst <- vec - perc_improve
    }

    inx <- which(rst < -1)
    if (length(inx) > 0)
        rst[inx] <- -1

    rst
}


#' Tumor burden response
#'
#' Best tumor burden response
#'
#'
#' @export
#'
tb_tb_resp <- function(dat_tb,
                        var_day  = "DAY",
                        var_tb   = "PCHG",
                        var_id   = "SUBJID",
                        ...) {
    f_tb <- function(tbl, id) {
        tb <- tbl %>%
            arrange(!!as.name(var_day))

        rst <- tb_pseu_or(tb[[var_tb]],
                          vec_t = tb[[var_day]],
                          ...)

        cbind(id, rst)
    }

    dat_tb %>%
        group_by(!! as.name(var_id)) %>%
        group_map(.f = ~ f_tb(.x, .y), .keep = TRUE) %>%
        bind_rows()
}

#' Median survival by best TB
#'
#' @param cut_tb Best response measured by PCHG
#'
#'  @export
#'
tb_tb_surv <- function(dta_tb, dta_surv, tb_cut = -0.3,
                       var_tb      = "PCHG",
                       var_id      = "SUBJID",
                       var_time    = "PFS_DAYS",
                       var_status  = "PFS_CNSR",
                       surv_quants = c(0.5, 0.7),
                       ...) {

    f_m <- function(d, bt) {
        fit <- tb_surv_fit(d, var_time, var_status,
                           fsurv = survfit, ...)
        rst            <- tkt_survfit_median(fit, surv_quants)
        rst$TB_Cut     <- tb_cut
        rst$BetterThan <- bt
        rst
    }

    dta <- dta_surv %>%
        left_join(dta_tb %>%
                  group_by(!!as.name(var_id)) %>%
                  summarize(PCHG_Min = min(!!as.name(var_tb),
                                           na.rm = TRUE)),
                  by = var_id)

    rst <- NULL
    inx <- which(dta$PCHG_Min <= tb_cut)
    if (length(inx) > 0) {
        rst <- rbind(rst, f_m(dta[inx, ], "Yes"))
    }

    if (length(inx) < nrow(dta)) {
        rst <- rbind(rst, f_m(dta[-inx, ], "No"))
    }

    rst
}


#' Observed TB mean, sd and missing rate
#'
#'
#'  @export
#'
tb_tb_obs <- function(dta_tb, dta_surv,
                      var_tb   = "PCHG",
                      var_day  = "DAY",
                      var_id   = "SUBJID",
                      var_time = "PFS_DAYS") {

    dta_surv$PFS <- dta_surv[[var_time]]
    dta_surv$ID  <- dta_surv[[var_id]]
    dta_tb$DAY   <- dta_tb[[var_day]]
    dta_tb$ID    <- dta_tb[[var_id]]

    rst <- dta_tb %>%
        left_join(dta_surv %>% select(ID, PFS), by = "ID") %>%
        filter(PFS >= DAY) %>%
        group_by(ARM, DAY) %>%
        summarize(N = n(),
                  PCHG_Mean = mean(PCHG, na.rm = TRUE),
                  PCHG_SD   = sd(PCHG, na.rm = TRUE))

    d_surv <- NULL
    vec_t  <- unique(rst$DAY)
    for (i in vec_t) {
        cur_d <- dta_surv %>%
            filter(PFS >= i) %>%
            group_by(ARM) %>%
            summarize(N_Surv = n())
        cur_d$DAY <- i

        d_surv <- rbind(d_surv, cur_d)
    }

    rst %>%
        left_join(d_surv,
                  by = c("DAY", "ARM")) %>%
        mutate(Miss_Rate = 1 - N / N_Surv)
}
