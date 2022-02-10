#' Get patient data
#'
#' Get turmor burden, survival and utility data of a specific patient
#'
#' @param date_dbl Date of database lock
#'
#' @export
#'
tb_get_pt <- function(id, imp_surv, dat_tb,
                      date_dbl = NULL,
                      t_ana    = NULL,
                      imp_inx  = 1,
                      reg_tb   = NULL,
                      by_arm   = TRUE,
                      ...) {

    d_surv <- imp_surv %>%
        filter(SUBJID == id &
               Imp    == imp_inx) %>%
        data.frame()

    d_tb <- dat_tb %>%
        filter(SUBJID == id &
               !is.na(PCHG)) %>%
        arrange(DAY)

    time_tb <- d_tb$DAY
    tb      <- d_tb$PCHG
    tb_obs  <- cbind(time_tb, tb)

    ## analysis time
    if (is.null(t_ana)) {
        stopifnot(!is.null(date_dbl))
        t_ana <- as.Date(date_dbl) - as.Date(d_tb$RANDT[1])
        t_ana <- as.numeric(t_ana)
    }


    ## event time
    time_event <- as.numeric(d_surv[1, c("IT_PFS", "IT_OS")])

    ## tumor burden
    if (!is.null(reg_tb)) {
        if (by_arm) {
            cur_reg <- reg_tb[[imp_inx]][[d_tb[1, "ARM"]]]
        } else {
            cur_reg <- reg_tb[[imp_inx]]
        }

        cur_pred <- tb_predict(id, cur_reg, ...)
        time_tb  <- cur_pred$time_tb
        tb       <- cur_pred$tb
    }

    ## get patient history
    rst_hist        <- tb_pt_hist(time_tb, tb, time_event, t_ana = t_ana, ...)
    rst_hist$id     <- id
    rst_hist$tb_obs <- tb_obs

    ## return
    rst_hist
}

#' Generate history and AUC of a patient
#'
#'
#' @export
#'
tb_pt_hist <- function(time_tb, tb, time_event, t_ana,
                       uti_gamma = c(0.2, 0.5),
                       locf = FALSE, ...) {

    ## remove NA pfs from event
    if (is.na(time_event[1])) {
        time_event <- time_event[-1]
        uti_gamma  <- uti_gamma[-1]
        ## os
        ind_event  <- 2
    } else {
        ## pfs os
        ind_event  <- 1:2
    }

    ## remove tb after events
    inx <- which(time_tb > min(time_event, na.rm = TRUE))
    if (0 < length(inx)) {
        time_tb <- time_tb[-inx]
        tb      <- tb[-inx]
    }

    ## tumor burden
    rst_x   <- c(time_tb, time_event)
    rst_y   <- c(tb, uti_gamma)
    rst_ind <- c(rep(0, length(tb)), ind_event)
    ntot    <- length(rst_x)

    ## locf setting
    if (locf) {
        last_t  <- time_tb[length(time_tb)]
        last_tb <- tb[length(tb)]

        rst_y[which(rst_x > last_t)] <- last_tb
    }

    ## adjust duration
    npost <- sum(rst_x >= t_ana)

    if (0 == npost) {
        pri_ana <- cbind(x = c(rst_x, t_ana),
                         y = c(rst_y, rst_y[ntot]),
                         z = c(rst_ind, 3))

        pos_ana <- NULL
    } else {
        inx <- ntot - npost
        y0  <- rst_y[inx]
        y1  <- rst_y[inx + 1]
        x0  <- rst_x[inx]
        x1  <- rst_x[inx + 1]
        y_ana <- y0 + (y1 - y0) * (t_ana - x0) / (x1 - x0)

        pri_ana <- cbind(x = c(rst_x[1:inx],   t_ana),
                         y = c(rst_y[1:inx],   y_ana),
                         z = c(rst_ind[1:inx], 3))

        pos_ana <- cbind(x = c(t_ana, rst_x[- (1:inx)]),
                         y = c(y_ana, rst_y[- (1:inx)]),
                         z = c(3,     rst_ind[- (1:inx)]))
    }

    ## utility
    uti_tb    <- 0
    uti_event <- 0
    for (i in 2:nrow(pri_ana)) {
        y0  <- pri_ana[i - 1, 2] + 1
        y1  <- pri_ana[i,     2] + 1
        x0  <- pri_ana[i - 1, 1]
        x1  <- pri_ana[i,     1]

        cur_u   <- (y1 + y0) * (x1 - x0) / 2

        if (0 == pri_ana[i - 1, 3]) {
            uti_tb  <- uti_tb + cur_u
        } else {
            uti_event <- uti_event + cur_u
        }
    }
    rst_uti <- uti_tb + uti_event

    if (min(time_event) < t_ana) {
        dur_event <- t_ana - min(time_event)
    } else {
        dur_event <- 0
    }
    dur_tb <- t_ana - dur_event

    ## return
    rst <- list(pri_ana     = pri_ana,
                pos_ana     = pos_ana,
                t_ana       = t_ana,
                dur_tb      = dur_tb,
                dur_event   = dur_event,
                utility     = rst_uti - t_ana,
                adj_utility = (rst_uti - t_ana) / t_ana,
                uti_tb      = uti_tb    - dur_tb,
                uti_event   = uti_event - dur_event,
                uti_gamma   = uti_gamma)
}

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
                       add_line       = TRUE,
                       add_tb_obs     = TRUE) {

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
        for (j in pt_his$uti_gamma) {
            rst <- rst +
                geom_hline(yintercept = j, lty = 2, col = "brown")
        }
    }

    ## title
    if (add_title) {
        g_tit <- paste(pt_his$id, " AUC = ", round(pt_his$utility, 2), sep = "")
        rst   <- rst +
            labs(title = g_tit)
    }

    ## lines and points
    rst <- rst +
        geom_point(aes(pch = Time), size = 1) +
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
        x_txt <- NULL
        l_txt <- NULL
        for (i in txt) {
            cur_d <- d_ana %>% filter(Time == i)

            if (0 == nrow(cur_d))
                next

            cur_x <- cur_d[1, "x"]
            x_txt <- c(x_txt, cur_x)
            l_txt <- c(l_txt, i)

            rst   <- rst +
                geom_vline(xintercept = cur_x,
                           lty = 2,
                           col = "brown")
        }

        rst <- rst +
            geom_text(data = data.frame(x = x_txt, l = l_txt),
                      aes(x = x, y = -0.9, label = l),
                      angle = 90, vjust = -0.5)
    }

    ## add observed tb pchg
    if (add_tb_obs) {
        rst <- rst +
            geom_point(data = data.frame(pt_his$tb_obs),
                       aes(x = time_tb, y = tb), pch = 2, col = "red")
    }


    ## return
    rst
}
