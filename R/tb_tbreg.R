#' Combine survival and tumor burden data
#'
#' Prepare data set for longitudinal tumor burden regression
#'
#' @param add_surv Whether include uti_gamma at the event time as a tumor burden
#'     measurement
#' @param f_conv Function for converting tumor burden
#'
#' @export
#'
tb_prepare_tb_dta <- function(dat_tb, imp_surv, imp_inx = 1,
                              uti_gamma = c("Progression" = 0.2,
                                            "Death"       = 0.5),
                              tb_reg_pow_conv = 1,
                              tb_reg_add_surv = FALSE,
                              ...) {

    imp_surv <- imp_surv %>%
        filter(Imp == imp_inx) %>%
        select(SUBJID, IT_Time, IT_Event)

    rst <- dat_tb %>%
        left_join(imp_surv) %>%
        mutate(reg_t     = DAY / IT_Time,
               reg_tb    = (1 + PCHG) * BASE,
               reg_gamma = (uti_gamma[IT_Event] + 1) * BASE)

    ## last point corresponding to event
    if (tb_reg_add_surv) {
        rst <- rst %>%
            group_by(SUBJID) %>%
            slice(n = 1) %>%
            mutate(reg_t  = 1,
                   reg_tb = reg_gamma) %>%
            rbind(rst)
    }

    ## convert tumor burden data
    rst$reg_tb <- rst$reg_tb ^ tb_reg_pow_conv

    ## return
    rst %>% arrange(SUBJID, reg_t)

}

#' Tumor burden regression
#'
#' Conduct tumor burden regression
#'
#'
#' @export
#'
tb_regression_single <- function(tb_surv, fml_tb = NULL,
                                 tb_reg_poly_order = 2,
                                 tb_reg_pow_conv = 1,
                                 ...) {

    ## time
    tb_reg_poly_coefs <- poly(tb_surv$reg_t, tb_reg_poly_order)
    tb_reg_poly_coefs <- attr(tb_reg_poly_coefs, "coefs")

    ## event time
    id_time <- tb_surv %>%
        select(SUBJID, IT_Time) %>%
        distinct()

    ## design matrix
    fml_1   <- paste("poly(reg_t, ", tb_reg_poly_order, ")")
    fml_2   <- paste("(", fml_1, "| SUBJID)")

    des_mat <- NULL
    if (!is.null(fml_tb)) {
        fml     <- paste("reg_tb", fml_tb, "+", fml_1, "+", fml_2)

        vars    <- all.vars(as.formula(fml_tb))
        tmp_tb  <- tb_surv %>%
            select(c("SUBJID", vars)) %>%
            distinct()

        des_mat <- model.matrix(as.formula(fml_tb),
                                data = tmp_tb)
        rownames(des_mat) <- tmp_tb$SUBJID
    } else {
        fml <- paste("reg_tb ~", fml_1, "+", fml_2)
    }

    ## regression
    fit_reg <- lmer(as.formula(fml), data = tb_surv)
    fit_rst <- coefficients(fit_reg)$SUBJID

    ## return
    list(tb_reg_poly_coefs  = tb_reg_poly_coefs,
         tb_reg_poly_order  = tb_reg_poly_order,
         tb_reg_pow_conv    = tb_reg_pow_conv,
         id_time            = id_time,
         des_mat            = des_mat,
         fit_rst            = fit_rst,
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
tb_regression <- function(dat_tb, imp_surv, ...) {
    imp_m <- max(imp_surv$Imp)
    arms  <- unique(dat_tb$ARM)

    rst   <- list()
    for (imp in seq_len(imp_m)) {
        cur_dta <- tb_prepare_tb_dta(dat_tb, imp_surv, imp, ...)
        cur_rst <- list()
        for (a in arms) {
            cur_arm <- cur_dta %>%
                filter(ARM == a)
            cur_rst[[a]] <- tb_regression_single(cur_arm, ...)
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

    pow_conv <- tb_reg_fit$tb_reg_pow_conv
    fit_rst  <- tb_reg_fit$fit_rst[id, ]
    des_mat  <- tb_reg_fit$des_mat[id, ]
    id_time  <- tb_reg_fit$id_time %>% filter(SUBJID == id)
    id_time  <- id_time$IT_Time[1]

    pred_t  <- seq(0, id_time, by = by_days)
    poly_t  <- poly(pred_t / id_time,
                    tb_reg_fit$tb_reg_poly_order,
                    coefs = tb_reg_fit$tb_reg_poly_coefs)

    par_t   <- fit_rst[- (1 : length(des_mat))]
    pred_y  <- sum(des_mat * fit_rst[1 : length(des_mat)])
    pred_y  <- pred_y + apply(poly_t, 1, function(x) sum(x * par_t))
    pred_y  <- pred_y ^ (1 / pow_conv)
    pred_y  <- pred_y / pred_y[1] - 1

    ## return
    list(time_tb = pred_t,
         tb      = pred_y)
}
