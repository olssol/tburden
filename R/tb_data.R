#' Get Analysis Data
#'
#'
#' @export
#'
tb_get_data <- function(raw_dat_rs, raw_dat_te,
                        first_n = 99999,
                        days_fu = 99999,
                        ...) {

    ## censor events
    f_censor <- function(dat, cut_date = cut_date_fu) {
        cut_days <- cut_date - dat[["RANDT"]] + 1
        inx      <- which(dat[["AVAL"]] > cut_days)

        dat[inx, "AVAL"]     <- cut_days[inx]
        dat[inx, "CNSR"]     <- 5
        dat[inx, "EVNTDESC"] <- "TRUNCATED AT ANALYSIS"

        dat
    }

    ## cut point
    dat_rs <- raw_dat_rs %>%
        filter(MITT1FL == "Y" &
               DTYPE   == ""  &
               EVAL    == "INDEPENDENT ASSESSOR" &
               PARAMCD == "SUMTLDLO")

    dat_id <- dat_rs %>%
        select(SUBJID, RANDT, ARM,
               BASE, AGE, SEX,
               ECOGGR1, STRATA1, P1TERTL) %>%
        distinct()  %>%
        mutate(SCL_AGE  = scale(AGE),
               SCL_BASE = scale(BASE)) %>%
        arrange(RANDT)

    cut_date_enroll <- dat_id[min(nrow(dat_id), first_n), "RANDT"]
    cut_date_fu     <- cut_date_enroll + days_fu

    ## tumor burden
    dat_tb <- dat_rs %>%
        select(SUBJID, VISIT, AVISITN, PCHG) %>%
        arrange(SUBJID, AVISITN, PCHG) %>%
        mutate(PCHG = if_else(1 == AVISITN, 0, PCHG / 100),
               DAY  = if_else(1 == AVISITN, 0, AVISITN)) %>%
        left_join(dat_id, by = "SUBJID") %>%
        filter(DAY <= cut_date_fu - RANDT + 1) %>%
        filter(RANDT <= cut_date_enroll) %>%
        filter(!is.na(BASE))

    ## survival
    dat_os <- raw_dat_te %>%
        filter(MITT1FL == "Y" &
               PARAMCD == "OS") %>%
        select(SUBJID, CNSR, EVNTDESC, AVAL) %>%
        f_censor() %>%
        rename(OS_CNSR  = CNSR,
               OS_EVENT = EVNTDESC,
               OS_DAYS  = AVAL)

    dat_pfs <- raw_dat_te %>%
        filter(MITT1FL == "Y"   &
               PARAMCD == "PFS" &
               EVAL    == "INDEPENDENT ASSESSOR") %>%
        select(SUBJID, CNSR, EVNTDESC, AVAL) %>%
        f_censor() %>%
        rename(PFS_CNSR  = CNSR,
               PFS_EVENT = EVNTDESC,
               PFS_DAYS  = AVAL)

    dat_surv <- dat_id %>%
        na.omit %>%
        left_join(dat_pfs, by = "SUBJID") %>%
        left_join(dat_os,  by = "SUBJID") %>%
        filter(RANDT <= cut_date_enroll)

    ## missing data
    warning(paste("Patients with missing covariates are",
                  "excluded from the result dataset.\n",
                  " Consider imputing the missing covariates",
                  "first (e.g., by mice) \n",
                  " before calling this function."))

    ## permuate data if necessary
    rst <- tb_data_permute(dat_tb, dat_surv, ...)

    ## add date_fu
    rst$cut_date_fu <- cut_date_fu

    rst
}


#' Permute data
#'
#'
#'
#' @export
#'
tb_data_permute <- function(dat_tb, dat_surv, permute = FALSE, seed = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    if (permute) {
        d_sub <- dat_tb %>%
            select(SUBJID, ARM, RANDT) %>%
            distinct()

        d_sub$SUBJID <- sample(d_sub$SUBJID)

        dat_tb <- dat_tb %>%
            select(-ARM, -RANDT) %>%
            left_join(d_sub, by = "SUBJID")

        dat_surv <- dat_surv %>%
            select(-ARM, -RANDT) %>%
            left_join(d_sub, by = "SUBJID")
    }

    if (!is.null(seed))
        set.seed(old_seed)

    list(dat_tb   = dat_tb,
         dat_surv = dat_surv)
}
