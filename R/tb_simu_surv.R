#' Simulate Survival time for all
#'
#'
#'@export
#'
tb_simu_surv <- function(n, hd_prog = 0.5, hd_death = 1, hd_cens = 1,
                         dur_enroll = 1, dur_study = 3,
                         seed = NULL) {

    ## random seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)


    t_enroll <- runif(n, min = 0, max = dur_enroll)
    t_prog   <- sapply(1:n, function(x) {
        rst <- tkt_simu_pwexp(hd_prog)
        rst[1]
    })

    t_death  <- sapply(1:n, function(x) {
        rst <- tkt_simu_pwexp(hd_death)
        rst[1]
    })

    t_cens   <- sapply(1:n, function(x) {
        rst <- tkt_simu_pwexp(hd_cens)
        rst[1]
    })

    ## set seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## summarize
    rst <- cbind(t_enroll, t_cens, t_prog, t_death)
    rst <- apply(rst, 1, function(x) {
        o_dur <- dur_study - x[1]

        if (x[3] < min(x[2], x[4], o_dur)) {
            o_prog <- x[3]
        } else {
            o_prog <- NA
        }

        if (x[4] < min(x[2], o_dur)) {
            o_death <- x[4]
        } else {
            o_death <- NA
        }

        if (is.na(o_death)) {
            o_cens <- min(x[2], o_dur)
        } else {
            o_cens <- NA
        }

        c(x, o_prog, o_death, o_cens, o_dur)
    })

    ## return
    rst           <- t(rst)
    colnames(rst) <- c("T_Enroll", "T_Cens",  "T_Prog", "T_Death",
                       "O_Prog",   "O_Death", "O_Cens", "O_Dur")
    rownames(rst) <- NULL

    data.frame(rst)
}
