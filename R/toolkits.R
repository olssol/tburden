#' Simulate from piecewise constant exponential
#'
#'
#' @examples
#' tkt_sim_pwexp(hazards = c("0.5" = 0.045, "Inf" = 0.025), offset = 0)
#'
#' @export
#'
tkt_simu_pwexp <- function(hazards, offset = 0) {

    if (1 == length(hazards)) {
        tte <- rexp(1,  hazards)
    } else {
        segments <- as.numeric(names(hazards))
        segments <- sort(segments) - offset

        inx    <- min(which(segments > 0))
        cur_t  <- 0
        flag   <- FALSE

        while (!flag) {
            cur_h   <- hazards[inx]
            cur_int <- rexp(1, cur_h)

            if ((cur_t + cur_int) <= segments[inx]) {
                flag <- TRUE
                tte  <- cur_t + cur_int
                break
            }

            cur_int <- segments[inx]
            inx     <- inx + 1
        }
    }

    ## return
    rst <- c(tte    = tte,
             offset = offset)

    return(rst)
}

#' Get P Values
#'
#'
#' @export
#'
tkt_pval <- function(est, sd_est) {
    z     <- est / sd_est
    p_val <- 2 * (1 - pnorm(abs(z)))

    p_val
}

#' Get censoring risk by annual dropout rate
#'
#'@export
#'
#'
tkt_lambda <- function(rate, tp = 1) {
    - log(1 - rate) / tp
}


#' Assign text to numeric vector
#'
#'
#' @export
#'
tkt_assign <- function(txt, prefix = "c(", suffix = ")") {

    rst <- NULL
    txt <- paste("rst <-", prefix, txt, suffix)

    tryCatch({
        eval(parse(text = txt))
    }, error = function(e) {
    })

    rst
}


#' Simulate Enrollment Time
#'
#'
#' @export
#'
tkt_simu_enroll <- function(ntot, enroll_duration, min_fu = 6,
                            date_bos = NULL, mth_to_days = 30.4) {
    rand_enroll <- runif(ntot, 0, enroll_duration)
    day_enroll  <- floor(rand_enroll * mth_to_days)

    rst <- data.frame(sid        = seq_len(ntot),
                      day_enroll = day_enroll)

    ## set up end of study time by minimum follow up
    if (!is.null(min_fu)) {
        day_eos     <- max(rand_enroll) + min_fu
        day_eos     <- day_eos * mth_to_days
        day_eos     <- floor(day_eos)
        rst$day_eos <- day_eos - day_enroll
    }


    ## set up dates in addition to days
    if (!is.null(date_bos)) {
        rst$date_bos    <- date_bos
        rst$date_eos    <- date_bos + day_eos
        rst$date_enroll <- date_bos + day_enroll
    }

    ## return
    rst
}

#' Simulate time to events in days
#'
#' @param ntot  total number of patients
#' @param hazard hazard for event
#' @param median_mth median survival in months
#' @param annual_drop annual drop rate
#'
#' @export
#'
tkt_rexp <- function(ntot,
                     hazard       = NULL,
                     median_mth   = 5,
                     annual_drop  = NULL,
                     mth_to_days  = 30.4,
                     take_floor   = TRUE) {

    if (is.null(hazard)) {
        if (!is.null(median_mth)) {
            hazard <- - log(0.5)  / median_mth
        } else {
            hazard <- - log(1 - annual_drop) / 12
        }
    }

    rand_event <- rexp(ntot, hazard) * mth_to_days

    if (take_floor)
        rand_event <- floor(rand_event)

    rand_event
}


#' Get table of median from survfit results
#'
#'
#' @export
#'
tkt_survfit_median <- function(fit, surv_quants = NULL) {
    rst <- rbind(summary(fit)$table)

    if (is.null(rownames(rst)))
        rownames(rst) <- seq_len(nrow(rst))

    rst <- data.frame(Group  = rownames(rst),
                      N      = rst[, "records"],
                      Quants = 0.5,
                      Days   = rst[, "median"])

    ## quantile
    if (!is.null(surv_quants)) {
        rst_q <- quantile(fit, surv_quants)$quantile
        rst_q <- rbind(rst_q)
        for (i in seq_len(nrow(rst_q))) {
            rst <- rbind(rst,
                         data.frame(Group    = rownames(rst)[i],
                                    N        = NA,
                                    Quants   = surv_quants,
                                    Days     = rst_q[i, ]))
        }
    }

    ## return
    rownames(rst) <- NULL
    rst
}

#' Get subset of data by ID
#'
#'
#' @export
#'
tkt_subset <- function(dta, id = "SUBJID", sub_n = NULL, sub_p = 0.1) {

    if (is.null(sub_p) & is.null(sub_n))
        return(dta)

    d_id <- unique(dta[[id]])
    if (!is.null(sub_p)) {
        sub_n <- ceiling(length(d_id) * sub_p)
    }

    d_id <- sample(d_id, sub_n)
    rst  <- dta %>%
        filter((!!sym(id)) %in% d_id)

    rst
}
