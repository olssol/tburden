#' Plot a patient
#'
#'
#'@export
#'
tb_plt_ind <- function(pt_his, ylim = NULL, xlim = NULL,
                       type = c("uti", "under"),
                       add_reflines_x = TRUE,
                       add_reflines_y = TRUE,
                       add_poly       = TRUE,
                       add_title      = TRUE,
                       add_line       = TRUE) {

    txt <- c("Progression", "Death", "Analysis")

    f_poly <- function(cur_p, cur_poly) {
        if (!add_poly)
            return(cur_p)

        colnames(cur_poly) <- c("x", "y")
        cur_poly           <- data.frame(cur_poly)

        if (min(cur_poly$y) < 0) {
            col_poly <- "blue"
        } else {
            col_poly <- "brown"
        }

        cur_p + geom_polygon(data = cur_poly, aes(x = x, y = y),
                            fill = col_poly,
                            alpha = 0.2)
    }


    type <- match.arg(type)

    d_ana <- data.frame(pt_his$pri_ana) %>%
        mutate(Type = "Prior") %>%
        rbind(data.frame(pt_his$pos_ana) %>%
              mutate(Type = "Post")) %>%
        mutate(Time = factor(z,
                             levels = 0:3,
                             labels = c("TB", txt)))

    ## limits
    if (is.null(ylim)) {
        ylim <- c(-1, max(d_ana$y) * 1.05)
    }

    if (is.null(xlim)) {
        xlim <- c(0, max(d_ana$x) * 1.05)
    }

    ## d_post
    d_post <- d_ana %>%
        filter(Type == "Post")


    ## plot
    rst <- ggplot(data = d_ana %>% filter(Type == "Prior"),
                  aes(x = x, y = y)) +
         labs(x = "Time", y = "Utility") +
         geom_hline(yintercept = 0, col = "brown")

    ## ploygon
    d_poly <- data.frame(pt_his$pri_ana[, c("x", "y")])
    if ("under" == type) {
        d_xm   <- max(d_poly$x)
        d_poly <- rbind(data.frame(x = 0, y = -1),
                        d_poly,
                        data.frame(x = d_xm, y = -1))

        rst <- rst %>% f_poly(d_poly)
    } else if ("uti" == type) {
        lx       <- d_poly$x[1]
        ly       <- d_poly$y[1]
        cur_poly <- c(lx, ly)
        for (i in 2:nrow(d_poly)) {
            cx <- d_poly$x[i]
            cy <- d_poly$y[i]

            if (cy * ly >= 0) {
                cur_poly <- rbind(cur_poly, c(cx, cy))
            } else {
                cx_0     <- lx + (cx - lx) * abs(ly / (cy - ly))
                cur_poly <- rbind(cur_poly, c(cx_0, 0))
                rst      <- rst %>% f_poly(cur_poly)

                ## reset new poly
                cur_poly <- rbind(c(cx_0, 0),
                                  c(cx,   cy))
            }

            lx <- cx
            ly <- cy
        }

        ## last poly
        cur_poly <- rbind(cur_poly, c(lx, 0))
        rst      <- rst %>% f_poly(cur_poly)
    }

    ## gamma values
    if (add_reflines_y) {
        for (j in pt_his$gamma) {
            rst <- rst +
                geom_hline(yintercept = j, lty = 2, col = "brown")
        }
    }

    ## title
    if (add_title) {
        g_tit <- paste(pt_his$id, " AUC = ", round(pt_his$uti, 2), sep = "")
        rst   <- rst +
            labs(title = g_tit)
    }

    ## lines and points
    rst <- rst +
        geom_point(aes(pch = Time), size = 2) +
        geom_point(data = d_post, aes(pch = Time)) +
        theme_bw() +
        lims(x = xlim, y = ylim)

    if (add_line) {
        rst <- rst +
            geom_line() +
            geom_line(data  = d_post, aes(x = x, y = y), lty = 2)
    }

    ## add vertical lines
    if (add_reflines_x) {
        for (i in txt) {
            cur_d <- d_ana %>% filter(Time == i)

            if (0 == nrow(cur_d))
                next

            cur_x <- cur_d[1, "x"]
            rst   <- rst +
                geom_vline(xintercept = cur_x, lty = 2, col = "brown")
        }
    }

    ## return
    rst
}

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

#' Summarize imputed survival data
#'
#'
#' @export
#'
tb_summary_imp <- function(imp_surv, dat_surv, inx_imp = NULL,
                           by_var = c("ARM")) {
    dat_surv <- imp_surv %>%
        left_join(dat_surv)

    if (!is.null(inx_imp)) {
        dat_surv <- dat_surv %>%
            filter(Imp == inx_imp)
    }

    dat_surv %>%
        group_by(!!as.name(by_var)) %>%
        summarize(Progression_Rate  = mean(is.na(IT_PFS)),
                  Progession_Mean   = mean(IT_PFS,   na.rm = T),
                  Progession_Median = median(IT_PFS, na.rm = T),
                  OS_Mean           = mean(IT_OS,    na.rm = T),
                  OS_Median         = median(IT_OS,  na.rm = T)
                  )
}
