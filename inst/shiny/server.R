options(shiny.maxRequestSize = 200 * 1024 ^ 2)

require(ggplot2)

shinyServer(function(input, output, session) {

    source("design_ui.R", local = TRUE);

    userLog          <- reactiveValues();
    userLog$data     <- NULL;


    ##--------------------------------------
    ##---------main-------------------------
    ##--------------------------------------
    output$mainpage <- renderUI({
       tab_main()
    })

    ##--------------------------------------
    ##---------exit-------------------------
    ##--------------------------------------
    observeEvent(input$close, {
        stopApp()})

    ##--------------------------------------
    ##---------data-------------------------
    ##--------------------------------------
    output$dt_surv <- DT::renderDataTable({
        dta <- get_data()

        if (is.null(dta))
            return(NULL)

        dta$dat_surv %>%
            select(SUBJID, ARM, RANDT, PFS_EVENT, PFS_DAYS, OS_EVENT, OS_DAYS)
    },
    selection = 'single',
    server    = TRUE,
    options   = list())

    output$dt_tb <- DT::renderDataTable({
        dat <- get_cur_tb()

        if (is.null(dat))
            return(NULL)
        dat %>%
            select(SUBJID, VISIT, DAY, PCHG)
    }, options = list(dom = 't'))

    output$dt_cov <- DT::renderDataTable({
        dat <- get_data()
        if (is.null(dat))
            return(NULL)

        id <- get_cur_id()
        if (is.null(id))
            return(NULL)

        dat$dat_surv %>%
            filter(SUBJID == id) %>%
            select(SUBJID, ARM, BASE, AGE, SEX, STRATA1, P1TERTL)
    }, options = list(dom = 't'))

    output$dt_impsurv <- DT::renderDataTable({
        get_cur_imp_surv()
    },
    selection = 'single',
    server    = TRUE,
    options = list(dom = 't'))

    ##--------------------------------------
    ##---------PLOTS------------------------
    ##--------------------------------------

    output$pltPt <- renderPlot({
        get_cur_plt()
    })

    output$pltTb <- renderPlot({
        dta <- get_data()

        if (is.null(dta))
            return(NULL)

        dat_tb <- dta$dat_tb
        id     <- get_cur_id()

        tb_plt_tb(dat_tb, id, by_var = input$inByvar)
    })

    ## observed survival
    output$pltPFS <- renderPlot({
        dta <- get_data()

        if (is.null(dta))
            return(NULL)

        tb_plt_km(dta$dat_surv, "PFS", by_var = input$inByvar)
    })

    output$pltOS <- renderPlot({
        dta <- get_data()

        if (is.null(dta))
            return(NULL)

        tb_plt_km(dta$dat_surv, "OS", by_var = input$inByvar)
    })

    ## imputed survival
    output$pltImpPFS <- renderPlot({
        dta <- get_data()

        if (is.null(dta))
            return(NULL)

        tb_plt_km_imp(dta$imp_surv, dta$dat_surv,
                      inx_imp = input$inImpInx, type = "PFS",
                      by_var = input$inByvar,
                      lim_x  = input$inSurvXlim)
    })

    output$pltImpOS <- renderPlot({
        dta <- get_data()

        if (is.null(dta))
            return(NULL)

        tb_plt_km_imp(dta$imp_surv, dta$dat_surv,
                      inx_imp = input$inImpInx, type = "OS",
                      by_var = input$inByvar,
                      lim_x  = input$inSurvXlim)
    })

    ##--------------------------------------
    ##---------TEXT-------------------------
    ##--------------------------------------

    output$txtHist <- renderPrint({
        print(get_cur_hist())
    })

    output$txtMsm <- renderPrint({
        dta <- get_data()

        if (is.null(dta))
            return(NULL)

        print(dta$fit_msm)
    })

    ##--------------------------------------
    ##---------SURVIVAL---------------------
    ##--------------------------------------
    output$dt_impsurv_summary <- DT::renderDataTable({
        get_impsurv_summary()
    }, options = list(dom = 't'))

    ##--------------------------------------
    ##---------Results----------------------
    ##--------------------------------------
    output$dt_rst <- DT::renderDataTable({
        dat <- get_data()
        if (is.null(dat))
            return(NULL)

        dat$results
    }, options = list(dom = 't'))


})
