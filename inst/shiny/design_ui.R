##-------------------------------------------------------------
##           UI FUNCTIONS
##-------------------------------------------------------------
tab_present <- function() {
    tabPanel("AUC",
             fluidRow(column(7,
                             wellPanel(h4("Survival Data"),
                                       div(DT::dataTableOutput("dt_surv"),
                                           style = "font-size:90%")),
                             wellPanel(
                                 h4("Survival Outcome"),
                                 DT::dataTableOutput("dt_impsurv")
                             ),
                             wellPanel(
                                 h4("Baseline"),
                                 DT::dataTableOutput("dt_cov"),
                                 h4("Tumor Burden"),
                                 DT::dataTableOutput("dt_tb")
                             )),
                      column(5,
                             wellPanel(
                                 h4("History and AUC"),
                                 plotOutput("pltPt", height = "500px"),
                                 sliderInput("inXlim",
                                             label = "",
                                             value = 0, min = 0, max = 5000, step = 100)
                             ),
                             wellPanel(h4("Options for Utility Plot"),
                                       fluidRow(
                                           column(5,
                                                  radioButtons("inAnaTime",
                                                               "Time for the Final Analysis",
                                                               choices = c("Calendar Time" = 1,
                                                                           "Fixed Time"    = 2)),
                                                  textInput("inDBL",
                                                            label = "Date of Analysis for (Calendar Time",
                                                            value = "2020-03-01"),
                                                  numericInput("inTana",
                                                               label = "Time for Analysis in Months (Fixed Time)",
                                                               value = 36)
                                                  ),
                                           column(3,
                                                  numericInput("inGammaPFS",
                                                               label = "Utility post PFS",
                                                               value = 0.2),
                                                  numericInput("inGammaOS",
                                                               label = "Utility post OS",
                                                               value = 0.5),
                                                  checkboxInput("inLocf",
                                                                label = "LOCF",
                                                                value = FALSE)
                                                  )
                                       )),
                             wellPanel(h4("Utility Details"),
                                       verbatimTextOutput("txtHist"))
                             )))
}

tab_upload <- function() {
    tabPanel("Upload Data",
             wellPanel(h4("Select R Data to Upload"),
                       fluidRow(
                           column(4,
                                  fileInput(inputId = 'inRdata',
                                            label   = 'Choose the analysis result R data file',
                                            accept  = '.Rdata')
                                  ))
                       ),
             wellPanel(h4("Create Pseudo Study"),
                       numericInput("inFirstn",
                                    label = "Keep the first N patients",
                                    value = 999999,
                                    width = "400px"),
                       numericInput("inFudays",
                                    label = "Follow up days since the last enrollment",
                                    value = 999999,
                                    width = "400px"),
                       )
             )
}

tab_results <- function() {
    tabPanel("Results",
             wellPanel(h4("Estimate and Confidence Interval"),
                       DT::dataTableOutput("dt_rst")
                       )
             )
}

tab_survival <- function() {
    tabPanel("TB and Survival",
             wellPanel(h4("Options for Tumor Burden and Survival Plot"),
                       selectInput(inputId = "inByvar",
                                   label   = "Group by",
                                   choices = c("ARM", "SEX",
                                               "STRATA1", "P1TERTL"),
                                   multiple = TRUE,
                                   selected = "ARM",
                                   width    = "400px")
                       ),
             wellPanel(h4("Tumor Burden by Time"),
                       plotOutput("pltTb")),
             wellPanel(fluidRow(
                           column(6,
                                  h4("Progression Free Survival"),
                                  plotOutput("pltPFS")),
                           column(6,
                                  h4("Overall Survival"),
                                  plotOutput("pltOS"))
                       )),
             wellPanel(h4("Imputed Survival"),
                       DT::dataTableOutput("dt_impsurv_summary")
                       ),
             wellPanel(h4("MSM Model Fitting Results"),
                       verbatimTextOutput("txtMsm"))
             )
}

##define the main tabset for beans
tab_main <- function() {
    tabsetPanel(type = "pills",
                id   = "mainpanel",
                tab_upload(),
                tab_present(),
                tab_survival(),
                tab_results()
                )
}



##-------------------------------------------------------------
##           DATA FUNCTIONS
##-------------------------------------------------------------

get_imp_data <- function(dat_surv, fml_surv) {
    ## multistate survival data
    msm_surv <- tb_msm_set_surv(dat_surv) %>%
        mutate(time = max(time, 10))

    ## fit imputation model
    msm_fit <- tb_msm_fit(msm_surv, fml_surv)

    ## imputation
    imp_surv <- tb_msm_imp(msm_fit, imp_m = imp_m)

    imp_surv
}

get_data <- reactive({
    ## ss <- load("./www/imp_data.Rdata")
    rst <- userLog$data
    if (is.null(rst)) {
        return(NULL)
    }

    first_n <- input$inFirstn
    days_fu <- input$inFudays
    if (!is.null(first_n) &
        !is.null(days_fu)) {

        ## only impute if smaller study created
        if (first_n < 1000 |
            days_fu < 10000) {
            ana_data <- tb_get_data(rst$raw_dat_rs,
                                    rst$raw_dat_te,
                                    first_n,
                                    days_fu)

            rst$dat_tb   <- ana_data$dat_tb
            rst$dat_surv <- ana_data$dat_surv
            rst$imp_surv <- get_imp_data(rst$dat_surv,
                                         rst$formula_surv)
        }
    }

    rst
})

## upload simulated results
observe({
    in_file <- input$inRdata

    if (!is.null(in_file)) {
        ss  <- load(in_file$datapath)
        isolate({
            userLog$data <- list(imp_surv      = rst_all$rst_orig$imp_surv,
                                 fit_msm       = rst_all$rst_orig$msm_fit$mdl_fit,
                                 dat_tb        = rst_all$rst_orig$params$dat_tb,
                                 dat_surv      = rst_all$rst_orig$params$dat_surv,
                                 formula_surv  = rst_all$rst_orig$params$fml_surv,
                                 results       = rst_all$summary,
                                 raw_dat_rs    = raw_dat_rs,
                                 raw_dat_te    = raw_dat_te)
        })
    }
})


get_cur_imp_surv <- reactive({
    dat <- get_data()
    if (is.null(dat))
        return(NULL)

    id <- get_cur_id()
    if (is.null(id))
        return(NULL)

    dat$imp_surv %>%
        filter(SUBJID == id) %>%
        select(Imp, SUBJID, IT_PFS, IT_OS)
})

get_cur_id <- reactive({
    dat <- get_data()
    if (is.null(dat))
        return(NULL)

    s <- input$dt_surv_rows_selected
    if (is.null(s))
        return(NULL)

    id <- dat$dat_surv[s, "SUBJID"]

    id
})

get_cur_imp_inx <- reactive({
    dat <- get_cur_imp_surv()
    if (is.null(dat))
        return(NULL)

    s <- input$dt_impsurv_rows_selected
    if (is.null(s))
        return(NULL)

    id <- dat[s, "Imp"]

    id
})


get_cur_tb <- reactive({
    dat <- get_data()
    if (is.null(dat))
        return(NULL)

    id <- get_cur_id()
    if (is.null(id))
        return(NULL)

    dat$dat_tb %>%
        filter(SUBJID == id)
})

get_cur_hist <- reactive({
    id <- get_cur_id()

    if (is.null(id))
        return(NULL)

    imp_inx <- get_cur_imp_inx()
    if (is.null(imp_inx))
        imp_inx <- 1

    dat <- get_data()
    if (is.null(dat))
        return(NULL)

    time_dbl  <- input$inDBL
    gamma_pfs <- input$inGammaPFS
    gamma_os  <- input$inGammaOS

    if (1 == input$inAnaTime) {
        t_ana <- NULL
    } else {
        t_ana <- input$inTana * 365.25 / 12
    }

    d_pt      <- tb_get_pt(id, dat$imp_surv, dat$dat_tb,
                           imp_inx  = imp_inx,
                           t_ana    = t_ana,
                           date_dbl = time_dbl,
                           gamma    = c(gamma_pfs, gamma_os),
                           locf     = input$inLocf)

    d_pt
})

## history of a patient
get_cur_plt <- reactive({
    cur_his <- get_cur_hist()
    if (is.null(cur_his))
        return(NULL)
    rst   <- tb_plt_ind(cur_his, type = "uti")
    x_max <- input$inXlim
    if (!is.na(x_max)) {
        if (x_max > 0)
            rst <- rst + coord_cartesian(xlim = c(0, x_max))
    }

    rst
})

get_impsurv_summary <- reactive({
    dat <- get_data()
    if (is.null(dat))
        return(NULL)

    by_var <- input$inByvar
    if (is.null(by_var))
        return(NULL)

    dat$imp_surv %>%
        left_join(dat$dat_tb) %>%
        group_by_(.dots = by_var) %>%
        summarize(Progression_Rate  = mean(is.na(IT_PFS)),
                  Progession_Mean   = mean(IT_PFS,   na.rm = T),
                  Progession_Median = median(IT_PFS, na.rm = T),
                  OS_Mean           = mean(IT_OS,   na.rm = T),
                  OS_Median         = median(IT_OS, na.rm = T)
                  )
})
