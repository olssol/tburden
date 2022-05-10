## ----------------------------------------------------------
##
##              FUNCTIONS FOR STUDY DESIGN PRESENTATION
##
## ----------------------------------------------------------

#' Response frequency table
#'
#' @export
#'
tb_des_freq_or <- function(dat_sub) {
    dat_sub <- dat_sub %>%
        select(ARM, SUBJID, OR, PosSub) %>%
        distinct()

    tbl_or <- dat_sub %>%
        count(ARM, OR) %>%
        group_by(ARM) %>%
        mutate(Rate = prop.table(n))

    tbl_pos <- dat_sub %>%
        count(ARM, PosSub, OR) %>%
        group_by(ARM, PosSub) %>%
        mutate(Rate = prop.table(n))

    list(tbl_or  = tbl_or,
         tbl_pos = tbl_pos)
}

#' TB summary
#'
#' @export
#'
tb_des_summary <- function(dat_sub) {
    dat_sub %>%
        group_by(ARM, DAY) %>%
        summarize(Mean = mean(PCHG),
                  SD   = sd(PCHG)) %>%
        filter(DAY > 0)
}


#' Spider plot of tumor burden
#'
#' @export
#'
tb_des_plt_tb_mean <- function(par_tb_0, par_tb_1, add_pred = FALSE) {

    fr <- function(par_tb) {
        if (add_pred) {
            rst_1 <- data.frame(DAY       = par_tb$vec_t,
                                PCHG      = par_tb$pred_y,
                                ARM       = par_tb$label,
                                Group     = "Fitted Mean")
        } else {
            rst_1 <- NULL
        }

        rst_2 <- data.frame(DAY       = par_tb$vec_t,
                            PCHG      = par_tb$mean_pchg,
                            ARM       = par_tb$label,
                            Group     = "Elicited Mean")

        rbind(rst_1, rst_2)
    }

    rst    <- rbind(fr(par_tb_0), fr(par_tb_1))
    by_var <- "ARM"
    s_fml  <- paste("~", paste(by_var, collapse = "+"))
    rst    <- ggplot(data = rst, aes(x = DAY)) +
        facet_wrap(as.formula(s_fml)) +
        theme_bw() +
        theme(legend.position = "bottom") +
        geom_point(aes(y = PCHG, col = Group)) +
        geom_line(aes(y = PCHG, col = Group)) +
        geom_hline(yintercept = 0, lty = 2) +
        geom_hline(yintercept = -0.3, lty = 2)

    rst
}


#' Spider plot of tumor burden
#'
#' @export
#'
tb_des_plt_tb_sub <- function(dat_sub, by_var = "ARM", tb_cut = -0.3,
                              sample_n = NULL) {


    s_fml <- paste("~", paste(by_var, collapse = "+"))

    if (!is.null(sample_n)) {
        smp_tb <- dat_sub %>%
            select(ARM, SUBJID) %>%
            distinct()

        inx <- sample(seq_len(nrow(smp_tb)), sample_n)
        dat_sub <- smp_tb[inx, ] %>%
            left_join(dat_sub,
                      by = c("ARM", "SUBJID"))
    }

    rst <- ggplot(data = dat_sub, aes(x = DAY, y = PCHG)) +
        facet_wrap(as.formula(s_fml)) +
        theme_bw() +
        theme(legend.position = "bottom") +
        geom_hline(yintercept = 0,    lty = 2) +
        geom_hline(yintercept = tb_cut, lty = 2) +
        geom_line(aes(group = SUBJID))


    rst
}


#' Waterfall plot of tumor burden
#'
#' @export
#'
tb_des_plt_tb_resp_waterf <- function(dat_sub, by_var = "ARM") {
    s_fml <- paste("~", paste(by_var, collapse = "+"))

    dat <- dat_sub %>%
        select(c("SUBJID", "PCHG_OR", by_var)) %>%
        distinct() %>%
        arrange(desc(PCHG_OR)) %>%
        group_by(!!as.name(by_var)) %>%
        mutate(Patients = row_number())

    rst   <- ggplot(data = dat, aes(x = Patients, y = PCHG_OR)) +
        facet_wrap(as.formula(s_fml)) +
        theme_bw() +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0,    lty = 2) +
        geom_hline(yintercept = -0.3, lty = 2) +
        labs( y = "Best PCHG")

    rst
}

#' Plot survival by response
#'
#'
#'
#' @export
#'
tb_des_plt_pfs <- function(dat_surv, dat_tb, ..., mth_to_days = 30.4) {
    dat_surv <- dat_surv %>%
        left_join(dat_tb %>%
                  select(ARM, SUBJID, OR, PosSub) %>%
                  distinct(),
                  by = c("SUBJID", "ARM")) %>%
        mutate(status  = 1,
               T_Event = T_Event / mth_to_days)

    print(
        dat_surv %>%
        group_by(ARM, PosSub) %>%
        summarize(n  = n(),
                  m  = median(T_Event),
                  m2 = mean(T_Event),
                  m3 = mean(surv_lambda)))

    by_trt      <- plot_km(dat_surv, ..., by_var = c("ARM"))
    by_trt_resp <- plot_km(dat_surv, ..., by_var = c("ARM", "OR"))
    by_trt_pos  <- plot_km(dat_surv, ..., by_var = c("ARM", "PosSub"))

    list(by_trt      = by_trt,
         by_trt_resp = by_trt_resp,
         by_trt_pos  = by_trt_pos)
}

#' Plot survival by response
#'
#'
#'
#' @export
#'
tb_des_plt_pfs_tbcut <- function(dta_tb, dta_surv, tb_cut = -0.3,
                                 var_tb      = "PCHG",
                                 var_id      = "SUBJID",
                                 var_time    = "PFS_DAYS",
                                 var_status  = "PFS_CNSR",
                                 ...) {

    dta <- dta_surv %>%
        left_join(dta_tb %>%
                  group_by(!!as.name(var_id)) %>%
                  summarize(PCHG_Min = min(!!as.name(var_tb),
                                           na.rm = TRUE)),
                  by = var_id) %>%
        mutate(BestTb = if_else(PCHG_Min < tb_cut,
                                  paste(">", tb_cut, sep = ""),
                                  paste("<", tb_cut, sep = "")))


    plot_km(dta, var_time, var_status,
            by_var = c("ARM", "BestTb"),
            ...)
}
