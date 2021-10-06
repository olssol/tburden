##-------------------------------------------------------------
##           UI FUNCTIONS
##-------------------------------------------------------------
tab_present <- function() {
    tabPanel("Data Visualization",
             fluidRow(column(7,
                             wellPanel(h4("Survival Data"),
                                       div(DT::dataTableOutput("dt_surv"),
                                           style = "font-size:90%")),
                             wellPanel(h4("Tumor Burden by Time"),
                                       plotOutput("pltTb")),
                             wellPanel(
                                 h4("Baseline"),
                                 DT::dataTableOutput("dt_cov"),
                                 h4("Tumor Burden"),
                                 DT::dataTableOutput("dt_tb")
                             )),
                      column(5,
                             wellPanel(
                                 h4("Survival Outcome"),
                                 DT::dataTableOutput("dt_impsurv")
                             ),
                             wellPanel(
                                 h4("History and AUC"),
                                 plotOutput("pltPt")),
                             wellPanel(
                                 fluidRow(
                                     column(4,
                                            textInput("inDBL",
                                                      label = "DBL Date",
                                                      value = "2020-03-01"),
                                            numericInput("inXlim",
                                                         label = "x-axis Range",
                                                         value = 0)
                                            ),
                                     column(4,
                                            numericInput("inGammaPFS",
                                                         label = "Utility post PFS",
                                                         value = 0.2),
                                            numericInput("inGammaOS",
                                                         label = "Utility post OS",
                                                         value = 0.5)
                                            )
                                 )
                             )
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
                       ))
}

tab_results <- function() {
    tabPanel("Results",
             wellPanel(h4("Estimate and Confidence Interval"),
                       DT::dataTableOutput("dt_rst")
                       ))
}

##define the main tabset for beans
tab_main <- function() {
    tabsetPanel(type = "pills",
                id   = "mainpanel",
                tab_upload(),
                tab_present(),
                tab_results()
                )
}


##-------------------------------------------------------------
##           DATA FUNCTIONS
##-------------------------------------------------------------
get_data <- reactive({
    ## ss <- load("./www/imp_data.Rdata")
    userLog$data
})

## upload simulated results
observe({
    in_file <- input$inRdata

    if (!is.null(in_file)) {
        ss  <- load(in_file$datapath)
        print(ss)
        isolate({
            userLog$data <- list(imp_surv = rst_all$rst_orig$imp_surv,
                                 dat_tb   = rst_all$rst_orig$params$dat_tb,
                                 dat_surv = rst_all$rst_orig$params$dat_surv,
                                 results  = rst_all$summary)
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

    d_pt      <- tb_get_pt(id, dat$imp_surv, dat$dat_tb, date_dbl = time_dbl,
                           imp_inx = imp_inx,
                           gamma = c(gamma_pfs, gamma_os))

    d_pt
})

get_cur_plt <- reactive({
    cur_his <- get_cur_hist()
    if (is.null(cur_his))
        return(NULL)
    rst   <- tb_plt_ind(cur_his)
    x_max <- input$inXlim
    if (!is.na(x_max)) {
        if (x_max > 0)
            rst <- rst + xlim(0, x_max)
    }

    rst
})
