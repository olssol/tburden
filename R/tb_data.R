#' Get Analysis Data
#'
#'
#' @export
#'
tb_get_data <- function(raw_dat_rs, raw_dat_te, first_n = 99999, days_fu = 99999) {

    ## censor events
    f_censor <- function(dat, prefix = "PFS") {
        fp <- function(vn) {
            paste(prefix, "_", vn, sep = "")
        }

        inx <- which(dat[[fp("DAYS")]] > days_fu)

        dat[inx, fp("DAYS")]  <- days_fu
        dat[inx, fp("CNSR")]  <- 5
        dat[inx, fp("EVENT")] <- "TRUNCATED AT ANALYSIS"

        dat
    }

    ## cut point
    dat_rs <- raw_dat_rs %>%
        filter(MITT1FL == "Y" &
               DTYPE   == ""  &
               EVAL    == "INDEPENDENT ASSESSOR" &
               PARAMCD == "SUMTLDLO")

    dat_id <- dat_rs %>%
        select(SUBJID, RANDT) %>%
        distinct()  %>%
        arrange(RANDT)

    cut_date <- dat_id[min(nrow(dat_id), first_n), "RANDT"]
    dat_rs   <- dat_rs %>%
        filter(RANDT <= cut_date)

    ## tumor burden
    dat_tb <- dat_rs %>%
        select(ARM, SUBJID, RANDT, VISIT, AVISITN, PCHG) %>%
        arrange(ARM, SUBJID, AVISITN, PCHG) %>%
        mutate(PCHG = if_else(1 == AVISITN, 0, PCHG / 100),
               DAY  = if_else(1 == AVISITN, 0, AVISITN)) %>%
        filter(DAY <= days_fu)

    ## survival
    dat_os <- raw_dat_te %>%
        filter(MITT1FL == "Y" &
               PARAMCD == "OS") %>%
        select(SUBJID, CNSR, EVNTDESC, AVAL) %>%
        rename(OS_CNSR  = CNSR,
               OS_EVENT = EVNTDESC,
               OS_DAYS  = AVAL) %>%
        f_censor(prefix = "OS")

    dat_pfs <- raw_dat_te %>%
        filter(MITT1FL == "Y"   &
               PARAMCD == "PFS" &
               EVAL    == "INDEPENDENT ASSESSOR") %>%
        select(SUBJID, CNSR, EVNTDESC, AVAL) %>%
        rename(PFS_CNSR  = CNSR,
               PFS_EVENT = EVNTDESC,
               PFS_DAYS  = AVAL) %>%
        f_censor(prefix = "PFS")

    dat_surv <- dat_rs %>%
        select(SUBJID, RANDT, ARM, BASE, AGE, SEX, ECOGGR1, STRATA1, P1TERTL) %>%
        distinct() %>%
        na.omit %>%
        left_join(dat_pfs) %>%
        left_join(dat_os) %>%
        select(SUBJID, RANDT, ARM,
               BASE, AGE, SEX, ECOGGR1, STRATA1, P1TERTL,
               PFS_DAYS,  OS_DAYS,
               PFS_CNSR,  OS_CNSR,
               PFS_EVENT, OS_EVENT)

    list(dat_tb   = dat_tb,
         dat_surv = dat_surv)
}
