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
        left_join(dat_surv, by = "SUBJID") %>%
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

## -----------------------------------------------------------------
##
##                SURVIVAL  PRESENTATION
##
## -----------------------------------------------------------------

#' Plot survival curve with area under the curve
#'
#' @export
#'
tb_plt_surv <- function(surv_f, t_dur = NULL,
                        type = c("rmf", "rmst", "none"),
                        y_lim = c(0, 1), x_lim = NULL) {

    type     <- match.arg(type)
    surv_dur <- tb_surv_cut(surv_f, t_dur)$surv_f_dur
    surv_km  <- tb_surv_cut(surv_f, x_lim)$surv_f_dur

    ## survival curves
    rst    <- ggplot(data = data.frame(Time = surv_km[, 1],
                                       Y    = surv_km[, 2]),
                     aes(x = Time, y = Y)) +
        labs(x = "Time", y = "Survival Probability") +
        ylim(y_lim) +
        theme_bw() +
        geom_step()

    if (is.null(x_lim)) {
        rst <- rst + xlim(c(0, x_lim))
    }

    if (type == "none")
        return(rst)

    ## ploygon
    if (!is.null(t_dur)) {
        y_dur    <- switch(type,
                           rmst = rbind(c(surv_dur[nrow(surv_dur), 1],
                                          0),
                                        c(0, 0)),
                           rmf  = c(surv_dur[nrow(surv_dur), 1], 1))

        surv_poly <- NULL
        for (i in 1:(nrow(surv_dur) - 1)) {
            surv_poly <- rbind(surv_poly,
                               surv_dur[i, ],
                               c(surv_dur[i + 1, 1], surv_dur[i, 2]))
        }

        surv_poly <- rbind(surv_poly, y_dur)
        rst <- rst + geom_polygon(data = data.frame(x = surv_poly[, 1],
                                                    y = surv_poly[, 2]),
                                  aes(x = x, y = y),
                                  alpha = 0.2)
    }

    rst
}

#' Plot patients by enrollment
#'
#' @export
#'
tb_plt_onstudy <- function(t_enroll, t_time, event, t_dur,
                           add_auc = FALSE, auc_k = 1.2,
                           add_lab = TRUE, size_lab = 8, hjust_lab = -1,
                           h = 0.4) {

    lab_e <- c("Censored", "Event", "Enrolled")
    dat <- data.frame(t_enroll = t_enroll,
                      time     = t_time,
                      event    = event) %>%
        arrange(t_enroll) %>%
        mutate(y     = row_number(),
               event = factor(event, 0:2, lab_e))

    rst <- ggplot(data = dat, aes(x = time, y = y)) +
        geom_point(aes(pch = event, color = event)) +
        geom_vline(xintercept = t_dur, lty = 2) +
        geom_vline(xintercept = 0, lty = 2) +
        geom_text(aes(x = 0, y = 5, label = "Study Started"),
                  angle = 90, vjust = -0.5) +
        geom_text(aes(x = t_dur, y = 5, label = "Study Finished"),
                  angle = 90, vjust = 1) +
        labs(xlim = c(-0.1, t_dur * 1.05), lty = 2) +
        theme_bw() +
        theme(
            axis.line    = element_blank(),
            axis.text    = element_blank(),
            axis.ticks   = element_blank(),
            axis.title   = element_blank(),
            panel.grid   = element_blank(),
            legend.title = element_blank(),
            legend.position = "bottom")

    ## geom_point(data = data.frame(x     = dat$t_enroll,
    ##                              y     = dat$y,
    ##                              event = factor(2, 0:2, lab_e)),
    ##            aes(x     = x,
    ##                y     = y,
    ##                pch   = event,
    ##                color = event)) +

    ## add line
    for (i in 1:nrow(dat)) {
        pt    <- dat[i, "time"]
        t_enr <- dat[i, "t_enroll"]
        pe    <- "Event" == dat[i, "event"]
        pd    <- data.frame(x = c(t_enr, pt),
                            y = c(i, i))
        rst <- rst +
            geom_line(data = pd, aes(x = x, y = y))

        if (add_auc) {
            if (pe) {
                rst <- rst +
                    geom_line(data = data.frame(x = c(pt, pt, t_dur),
                                               y = c(i, i + h, i + h)),
                              aes(x = x, y = y))
                pt_imp <- pt
            } else {
                pt_imp <- runif(1, pt, pt + auc_k * (t_dur - pt))
                pt_imp <- min(pt_imp, t_dur)
                rst    <- rst +
                    geom_line(data = data.frame(
                                  x = c(pt, pt_imp, pt_imp, t_dur),
                                  y = c(i,  i, i + h, i + h)),
                              aes(x = x, y = y), lty = 2)
            }

            rst <- rst +
                geom_polygon(data = data.frame(
                                 x = c(pt_imp, pt_imp, t_dur, t_dur),
                                 y = c(i, i + h, i + h, i)),
                             aes(x = x, y = y),
                             fill = "gray30",
                             alpha = 0.2)
        }
    }

    ## add pt label
    if (add_lab) {
        dat_lab <- data.frame(x    = dat$t_enroll,
                              y    = dat$y,
                              labs = paste("P", dat$y, sep = ""))

        rst <- rst +
            geom_label(data = dat_lab,
                       aes(x = x, y = y, label = labs),
                       hjust = hjust_lab,
                       size  = size_lab)
    }

    ## return
    rst
}
