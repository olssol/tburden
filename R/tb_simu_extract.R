#' Extract model parameters
#'
#'
#' @export
#'
tb_simu_extract <- function(dat_surv, dat_tb, fml_surv, fml_tb, ...) {
    ## pfs
    es_surv_fit <- tb_weibull_fit(dat_surv,
                                  fml_surv  = fml_surv,
                                  by_arm    = FALSE,
                                  ...)

    ## imputed survival
    es_surv_imp <- tb_surv_imp(es_surv_fit, imp_m = 1, ...)

    ## tumor burden
    es_tb_fit <- tb_regression(dat_tb, es_surv_imp,
                               fml_tb          = fml_tb,
                               tb_reg_poly_raw = TRUE,
                               by_arm          = FALSE,
                               ...)

    ## return
    list(surv_fit = es_surv_fit,
         tb_fit   = es_tb_fit,
         fml_tb   = fml_tb,
         fml_surv = fml_surv)
}
