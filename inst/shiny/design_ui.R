##-------------------------------------------------------------
##           UI FUNCTIONS
##-------------------------------------------------------------
tab_upload <- function() {
    tabPanel("Upload Data",
             wellPanel(h4("Select R Data to Upload"),
                       fluidRow(
                           column(4,
                                  fileInput(inputId = 'inRdata',
                                            label   = 'Choose result R data file',
                                            accept  = '.Rdata')
                                  ))
                       ),
             conditionalPanel(condition = "output.loadcomplete",
                              wellPanel(h4("Results: Estimate and Confidence Interval"),
                                        DT::dataTableOutput("dt_rst")
                                        ),
                              wellPanel(h4("Settings"),
                                        verbatimTextOutput("txtSetting")))
             ## wellPanel(h4("Create Pseudo Study"),
             ##           numericInput("inFirstn",
             ##                        label = "Keep the first N patients",
             ##                        value = 999999,
             ##                        width = "400px"),
             ##           numericInput("inFudays",
             ##                        label = "Follow up days since the last enrollment",
             ##                        value = 999999,
             ##                        width = "400px"),
             ##           )
             )
}

panel_auc_options <- function() {
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
              ))
}

tab_present <- function() {
    tabPanel("AUC",
             fluidRow(
                 wellPanel(h4("Patient Survival Data"),
                           div(DT::dataTableOutput("dt_surv"),
                               style = "font-size:90%"))),
             fluidRow(column(6,
                             wellPanel(
                                 h4("Baseline"),
                                 DT::dataTableOutput("dt_cov")),
                             wellPanel(
                                 h4("Survival Outcome"),
                                 DT::dataTableOutput("dt_impsurv")),
                             wellPanel(
                                 h4("Tumor Burden"),
                                 DT::dataTableOutput("dt_tb"))
                             ),
                      column(6,
                             wellPanel(
                                 h4("History and AUC"),
                                 plotOutput("pltPt", height = "500px"),
                                 sliderInput("inXlim",
                                             label = "",
                                             value = 0,
                                             min   = 0,
                                             max   = 5000,
                                             step  = 100),
                                 checkboxInput("inRegTb",
                                               "By observed TB",
                                               value = FALSE)
                             ))),
             fluidRow(
                 wellPanel(h4("Utility Details"),
                           verbatimTextOutput("txtHist"))
             ))
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
             wellPanel(h4("Observed Survival"),
                       fluidRow(
                           column(6,
                                  h4("Progression Free Survival"),
                                  plotOutput("pltPFS")),
                           column(6,
                                  h4("Overall Survival"),
                                  plotOutput("pltOS"))
                       )),
             wellPanel(h4("Imputed Survival"),
                       fluidRow(
                           column(4,
                                  selectInput(inputId = "inImpInx",
                                              label   = "Index of Imputation",
                                              choices = 1:5,
                                              width   = "400px")),
                           column(4,
                                  sliderInput("inSurvXlim",
                                              label = "",
                                              value = 0, min = 0,
                                              max = 8000, step = 100)
                                  )),
                       fluidRow(
                           column(6,
                                  h4("Imputed Progression Free Survival"),
                                  plotOutput("pltImpPFS")),
                           column(6,
                                  h4("Imputed Overall Survival"),
                                  plotOutput("pltImpOS"))
                       ),
                       DT::dataTableOutput("dt_impsurv_summary")
                       ),
             wellPanel(h4("MSM Model Fitting Results"),
                       verbatimTextOutput("txtMsm"))
             )
}

tab_corr <- function() {
    xx <- c("utility", "adj_utility",
            "uti_tb", "uti_event", "t_ana")

    tabPanel("Correlation",
             wellPanel(h4("Select Subject-Level Measurements"),
                       fluidRow(
                           column(3,
                                  selectInput(inputId = "inCorX",
                                              label   = "X",
                                              choices = xx,
                                              selected = "uti_tb")
                                  ),
                           column(3,
                                  selectInput(inputId = "inCorY",
                                              label   = "Y",
                                              choices = xx,
                                              selected = "uti_event")))
                       ),
             wellPanel(h4("Correlation"),
                       plotOutput("pltCorr", height = "800px"))
             )
}


##define the main tabset for beans
tab_main <- function() {
    tabsetPanel(type = "pills",
                id   = "mainpanel",
                tab_upload(),
                tab_survival(),
                tab_present(),
                tab_corr()
                ## ,tab_results()
                )
}



##-------------------------------------------------------------
##           DATA FUNCTIONS
##-------------------------------------------------------------
## upload simulated results
observe({
    in_file <- input$inRdata

    if (!is.null(in_file)) {
        ss  <- load(in_file$datapath)
        print("load data...")
        isolate({
            userLog$data <- tb_extract_rst(rst_all)
        })
    }
})

get_data <- reactive({
    rst <- userLog$data
    if (is.null(rst)) {
        return(NULL)
    }

    if (0) {
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
    }

    rst
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

## get patient history
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

    if (0) {
        time_dbl  <- input$inDBL
        gamma_pfs <- input$inGammaPFS
        gamma_os  <- input$inGammaOS

        if (1 == input$inAnaTime) {
            t_ana <- NULL
        } else {
            t_ana <- input$inTana * 365.25 / 12
        }
    }

    if (input$inRegTb) {
        reg_tb <- NULL
    } else {
        reg_tb <- dat$reg_tb
    }

    d_pt <- tb_get_pt(id,
                      imp_surv  = dat$imp_surv,
                      dat_tb    = dat$dat_tb,
                      imp_inx   = imp_inx,
                      t_ana     = NULL,
                      date_dbl  = dat$date_dbl,
                      uti_gamma = dat$uti_gamma,
                      reg_tb    = reg_tb)

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

    inx_imp <- input$inImpInx
    tb_summary_imp(dat$imp_surv, dat$dat_surv, inx_imp, by_var)
})



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
