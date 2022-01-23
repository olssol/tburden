#' Spider plot of tumor burden
#'
#' @export
#'
tb_plt_tb <- function(dat_tb, sel_ids = NULL, by_var = c("ARM")) {

    s_fml <- paste("~", paste(by_var, collapse = "+"))

    rst <- ggplot(data = dat_tb, aes(x = DAY, y = PCHG)) +
        geom_line(aes(group = SUBJID), col = "brown") +
        facet_wrap(as.formula(s_fml)) +
        theme_bw() +
        theme(legend.position = "none")

    ## selected pt
    d_sel <- dat_tb %>%
        filter(SUBJID %in% sel_ids)

    if (nrow(d_sel) > 0)
        rst <- rst +
            geom_line(data = d_sel,
                      aes(x = DAY, y = PCHG, group = SUBJID),
                      col = "green",
                      lwd = 1.5)

    rst

}


#' Survival curves
#'
#' @export
#'
tb_plt_km <- function(dat_surv, type = c("PFS", "OS"), ...) {

    type       <- match.arg(type)
    var_status <- paste(type, "_", "CNSR", sep = "")
    var_time   <- paste(type, "_", "DAYS", sep = "")

    plot_km(dat_surv, var_time, var_status, lab_y = type, ..., )
}


#' Survival curves for imputed survival
#'
#' @export
#'
tb_plt_km_imp <- function(imp_surv, dat_surv, inx_imp = NULL,
                          type = c("PFS", "OS"), ...) {
    type     <- match.arg(type)
    dat_surv <- imp_surv %>%
        left_join(dat_surv) %>%
        mutate(status = 0)

    if (!is.null(inx_imp)) {
        dat_surv <- dat_surv %>%
            filter(Imp == inx_imp)
    }

    stopifnot(nrow(dat_surv) > 0)

    if ("PFS" == type) {
        dat_surv$time <- apply(dat_surv[, c("IT_PFS", "IT_OS")], 1,
                               function(x) min(x, na.rm = TRUE))
    } else {
        dat_surv$time <- dat_surv$IT_OS
    }

    plot_km(dat_surv, "time", "status", lab_y = type, ...)
}



#' Plot correlation of utilities
#'
#'
#'
#'
#' @export
#'
tb_plt_estimate <- function(rst_estimate, var1 = "uti_tb", var2 = "uti_event") {

    rst_estimate$x <- rst_estimate[[var1]]
    rst_estimate$y <- rst_estimate[[var2]]

    sum_lm <- rst_estimate %>%
        group_by(imp, ARM) %>%
        summarize(R2 =  cor(x, y)) %>%
        mutate(R2 = round(R2, 3))

    ggplot(data = rst_estimate, aes(x = x, y = y)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        theme_bw() +
        facet_grid(imp ~ ARM) +
        labs(x = var1, y = var2) +
        geom_label(data = sum_lm,
                   aes(x = -Inf, y = Inf,
                       label = paste("R2 = ", R2, sep = " ")),
                   hjust = 0, vjust = 1)
}
