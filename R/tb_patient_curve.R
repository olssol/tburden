#' TB linear prediction points
#'
#' TB prediction based on the observed the utility curve
#'
#' @export
#'
tb_pt_curve_emp_single <- function(tb_mat, ts) {
    rst <- matrix(NA, length(ts), 3)
    for (i in seq_len(length(ts))) {
        cur_rst  <- tb_pt_insert(tb_mat, ts[i], 0, method = "extrap")
        rst[i, ] <- cur_rst$ins[1, ]
    }

    rst
}

#' TB linear prediction points
#'
#' TB prediction based on the observed the utility curve
#'
#' @export
#'
tb_pt_curve_emp <- function(tb_mat, ts = 1:300) {
    f_ins <- function(dat, grp) {
        rst <- tb_pt_curve_emp_single(as.matrix(dat[, c("x", "y", "z")]),
                                     ts)

        data.frame(SUBJID = grp$SUBJID[1],
                   imp    = grp$imp[1],
                   x      = rst[, 1],
                   y      = rst[, 2])
    }


    tb_mat %>%
        group_by(SUBJID, imp) %>%
        group_map(.f = ~f_ins(dat = .x, grp = .y),
                  .keep = TRUE) %>%
        rbindlist()
}

#' TB curves summarization
#'
#' Summarize TB curves with Jackknife variances
#'
#' @export
#'
tb_pt_curve_summary <- function(dta_curve, quant = NULL) {

    f_stat <- function(vec) {
        if (is.null(quant)) {
            rst <- mean(vec)
        } else {
            rst <- quantile(vec, quant)
        }

        rst
    }

    f_sum <- function(id, arm) {
        dta_curve %>%
            filter(SUBJID != id &
                   ARM    == arm) %>%
            group_by(ARM, x) %>%
            summarize(jk_y = f_stat(y))
    }

    ## overall mean
    rst <- dta_curve %>%
        group_by(ARM, x) %>%
        summarize(y = f_stat(y)) %>%
        mutate(var = 0)

    ## jackknife
    sid <- dta_curve %>%
        select(ARM, SUBJID) %>%
        distinct()

    for (i in seq_len(nrow(sid))) {
        print(i)
        cur_jk  <- f_sum(sid[i, "SUBJID"],
                         sid[i, "ARM"])

        rst <- rst %>%
            left_join(cur_jk, by = c("ARM", "x"))  %>%
            mutate(var = if_else(is.na(jk_y),
                                 var,
                                 var + (y - jk_y)^2)) %>%
            select(- jk_y)
    }

    rst %>%
        left_join(rst %>%
                  group_by(ARM) %>%
                  summarize(n = n()),
                  by = "ARM") %>%
        mutate(var = var * (n - 1) / n ) %>%
        mutate(UB = y + 1.96 * sqrt(var),
               LB = y - 1.96 * sqrt(var))
}
