#' Plot a patient
#'
#'
#'@export
#'
tb_plt_ind <- function(pt_his, ylim = NULL, xlim = NULL) {

    d_ana <- data.frame(pt_his$pri_ana) %>%
        mutate(Type = "Prior") %>%
        rbind(data.frame(pt_his$pos_ana) %>%
              mutate(Type = "Post")) %>%
        mutate(Time = factor(z,
                          levels = 0:3,
                          labels = c("TB", "Progression", "Death", "Analysis")))

    ## limits
    if (is.null(ylim)) {
        ylim <- c(-1, max(d_ana$y) * 1.05)
    }

    if (is.null(xlim)) {
        xlim <- c(0, max(d_ana$x) * 1.05)
    }

    ## ploygon
    d_poly <- data.frame(pt_his$pri_ana[, c("x", "y")])
    d_xm   <- max(d_poly$x)
    d_poly <- rbind(data.frame(x = 0, y = -1),
                    d_poly,
                    data.frame(x = d_xm, y = -1))

    ## d_post
    d_post <- d_ana %>%
        filter(Type == "Post")

    ## title
    g_tit <- paste(pt_his$id, " AUC = ", round(pt_his$auc, 2), sep = "")

    ## plot
    rst <- ggplot(data = d_ana %>% filter(Type == "Prior"),
                  aes(x = x, y = y)) +
        geom_line() +
        geom_point(aes(pch = Time)) +
        geom_polygon(data = d_poly, aes(x = x, y = y), alpha = 0.2) +
        geom_line(data = d_post, aes(x = x, y = y), lty = 2) +
        geom_point(data = d_post, aes(pch = Time)) +
        theme_bw() +
        lims(x = xlim, y = ylim) +
        labs(x = "Time", y = "Utility", title = g_tit)

    for (j in 1:length(pt_his$gamma)) {
        rst <- rst +
            geom_hline(yintercept = pt_his$gamma[j], lty = 2, col = "red")
    }

    rst
}

#' Spider plot of tumor burden
#'
#' @export
#'
tb_plt_tb <- function(dat_tb) {
    if (is.null(dat_tb$selected))
        dat_tb$selected <- FALSE

    ggplot(data = dat_tb, aes(x = DAY, y = PCHG)) +
        geom_line(aes(col = selected, group = SUBJID)) +
        facet_wrap(~ARM) +
        theme_bw() +
        theme(legend.position = "none")
}

#' Spider plot of tumor burden
#'
#' @export
#'
tb_plt_tb <- function(dat_tb, sel_ids = NULL) {

    rst <- ggplot(data = dat_tb, aes(x = DAY, y = PCHG)) +
        geom_line(aes(group = SUBJID), col = "brown") +
        facet_wrap(~ARM) +
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

#' Swimmer plot of survival
#'
#' @export
#'
tb_plt_surv <- function(dat_surv) {
}
