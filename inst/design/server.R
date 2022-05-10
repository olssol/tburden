require(ggplot2)

## Server
shinyServer(function(input, output, session) {

    ##---------------------------------------------------------------
    ##---------LOAD FUNCTIONS----------------------------------------
    ##---------------------------------------------------------------
    ## source("shiny_tools.R",  local = TRUE)

    ## ----------------------------------------------------
    ##           VALUES
    ## ----------------------------------------------------
    userLog          <- reactiveValues()
    userLog$labels   <- NULL
    userLog$sub_rand <- cbind(rnorm(1000), rnorm(1000))
    userLog$par_tb   <- NULL
    userLog$par_surv <- NULL

    ## read tumor burden parameters
    observeEvent(input$btnUpdateTb, {
        if (0 == input$btnUpdateTb)
            return(NULL)

        par_tb               <- list()
        par_tb$vec_t         <- tkt_assign(input$inVecT)
        par_tb$mean_pchg_ctl <- tkt_assign(input$inMeanCtl)
        par_tb$mean_pchg_trt <- tkt_assign(input$inMeanTrt)
        par_tb$sd_pchg       <- tkt_assign(input$inTbSig)
        par_tb$tb_rand_ctl   <- tkt_assign(input$inTbRndCtl)
        par_tb$tb_rand_trt   <- tkt_assign(input$inTbRndTrt)
        par_tb$labels        <- c(input$inLblCtl,
                                  input$inLblTrt)

        ## check if there is null value
        if (length(par_tb) < 7)
            return(NULL)

        userLog$par_tb <- par_tb
    })

    ## tumor burden curves
    par_tb_dta <- reactive({
        pb    <- userLog$par_tb
        srand <- userLog$sub_rand

        if (is.null(pb) |
            is.null(srand))
            return(NULL)

        par_0 <- tb_des_tpar(pb$mean_pchg_ctl,
                             vec_t = pb$vec_t,
                             label = pb$labels[1])
        par_1 <- tb_des_tpar(pb$mean_pchg_trt,
                             vec_t = pb$vec_t,
                             label = pb$labels[2])

        rst_0 <- tb_des_simu_tb(vec_t       = pb$vec_t,
                                par_tb      = par_0$par_tb,
                                sub_rand    = srand[, 1],
                                par_tb_rand = pb$tb_rand_ctl,
                                par_tb_sig  = pb$sd_pchg,
                                label       = par_0$label)

        rst_1 <- tb_des_simu_tb(vec_t       = pb$vec_t,
                                par_tb      = par_1$par_tb,
                                sub_rand    = srand[, 2],
                                par_tb_rand = pb$tb_rand_trt,
                                par_tb_sig  = pb$sd_pchg,
                                label       = par_1$label)

        rst_tb <- rbind(rst_0, rst_1)

        ## return
        rst <- list(par_tb_ctl = par_0,
                    par_tb_trt = par_1,
                    dta_tb     = rst_tb)
        rst
    })

    ## sample tumor burden curves
    par_tb_smps <- reactive({
        ptb <- par_tb_dta()
        if (is.null(ptb))
            return(NULL)

        n      <- input$inN
        smp_tb <- ptb$dta_tb %>%
            select(ARM, SUBJID) %>%
            distinct()

        inx <- sample(seq_len(nrow(smp_tb)), n)
        rst <- smp_tb[inx, ] %>%
            left_join(ptb$dta_tb,
                      by = c("ARM", "SUBJID"))

        rst
    })

    ## tb summary tables
    tbl_tb <- reactive({
        ptb <- par_tb_dta()
        if (is.null(ptb))
            return(NULL)

        tb_des_freq_or(ptb$dta_tb)
    })

    ## read survival parameters
    observeEvent(input$btnUpdateSurv, {
        if (0 == input$btnUpdateSurv)
            return(NULL)

        par_surv               <- list()
        par_surv$median_ctl    <- tkt_assign(input$inMedianCtl)
        par_surv$median_trt    <- tkt_assign(input$inMedianTrt)
        par_surv$surv_rand_ctl <- tkt_assign(input$inSurvRndCtl)
        par_surv$surv_rand_trt <- tkt_assign(input$inSurvRndTrt)

        ## check if there is null value
        if (length(par_surv) < 4)
            return(NULL)

        userLog$par_surv <- par_surv
    })

    ## survival results
    par_surv_dta <- reactive({
        pb    <- userLog$par_tb
        pv    <- userLog$par_surv
        srand <- userLog$sub_rand

        if (is.null(pb) |
            is.null(pv) |
            is.null(srand))
            return(NULL)

        par_surv_0   <- tb_des_spar(median_surv   = pv$median_ctl,
                                    par_surv_rand = pv$surv_rand_ctl,
                                    label         = pb$labels[1])
        par_surv_1   <- tb_des_spar(median_surv   = pv$median_trt,
                                    par_surv_rand = pv$surv_rand_trt,
                                    label         = pb$labels[2])

        rst_surv_0   <- tb_des_simu_pfs(par_surv_0,
                                        sub_rand = srand[, 1])
        rst_surv_1   <- tb_des_simu_pfs(par_surv_1,
                                        sub_rand = srand[, 2])

        rst_surv     <- rbind(rst_surv_0, rst_surv_1)

        ## return
        list(par_surv_ctl = par_surv_0,
             par_surv_trt = par_surv_1,
             dta_surv     = rst_surv)
    })


    ## survival curves
    par_surv_plt <- reactive({
        psurv <- par_surv_dta()
        ptb   <- par_tb_dta()
        xlim  <- tkt_assign(input$inSurvXlim)

        if (is.null(xlim)) {
            xlim <- 54
        }

        if (is.null(psurv) | is.null(ptb))
            return(NULL)

        rst <- tb_des_plt_pfs(psurv$dta_surv,
                              ptb$dta_tb,
                              var_time   = "T_Event",
                              var_status = "status",
                              ind_event  = 1,
                              lim_x      = xlim)

        rst
    })

    ## --------------------------------------------------------
    ##                    UI Update
    ## --------------------------------------------------------

    ## tb mean curves
    output$pltMeanTb <- renderPlot({
        ptb <- par_tb_dta()
        if (is.null(ptb))
            return(NULL)

        tb_des_plt_tb_mean(ptb$par_tb_ctl,
                           ptb$par_tb_trt)
    })

    ## tb individual curves
    output$pltSubTb <- renderPlot({
        smp_tb <- par_tb_smps()
        if (is.null(smp_tb))
            return(NULL)

        tb_des_plt_tb_sub(smp_tb)
    })

    ## waterfall
    output$pltSubTbWf <- renderPlot({
        smp_tb <- par_tb_smps()
        if (is.null(smp_tb))
            return(NULL)

        tb_des_plt_tb_resp_waterf(smp_tb %>%
                                  filter(OR == "Yes"))
    })


    ## response rates
    output$tblResp <- renderTable({
        tbl <- tbl_tb()
        if (is.null(tbl))
            return(NULL)

        tbl$tbl_or
    })

    ## response rates
    output$tblRespPos <- renderTable({
        tbl <- tbl_tb()
        if (is.null(tbl))
            return(NULL)

        tbl$tbl_pos
    })

    ## tb summary
    output$tblTb <- renderTable({
        ptb <- par_tb_dta()
        if (is.null(ptb))
            return(NULL)

        tb_des_summary(ptb$dta_tb)
    })

    output$pltSurvArm <- renderPlot({
        plt <- par_surv_plt()
        if (is.null(plt))
            return(NULL)

        plt$by_trt$plot
    })

    output$pltSurvArmResp <- renderPlot({
        plt <- par_surv_plt()
        if (is.null(plt))
            return(NULL)

        plt$by_trt_resp$plot
    })

    output$pltSurvArmPos <- renderPlot({
        plt <- par_surv_plt()
        if (is.null(plt))
            return(NULL)

        plt$by_trt_pos$plot
    })

    output$tblSurvArm <- renderTable({
        plt <- par_surv_plt()
        if (is.null(plt))
            return(NULL)

        tkt_survfit_median(plt$by_trt$fit)
    })

    output$tblSurvArmResp <- renderTable({
        plt <- par_surv_plt()
        if (is.null(plt))
            return(NULL)

        tkt_survfit_median(plt$by_trt_resp$fit)
    })

    output$tblSurvArmPos <- renderTable({
        plt <- par_surv_plt()
        if (is.null(plt))
            return(NULL)

        tkt_survfit_median(plt$by_trt_pos$fit)
    })

    output$txtInfo <- renderPrint({
        print("-----------TB---------------")
        print(userLog$par_tb)
        print("-----------------------------")
        print(par_tb_dta()$par_tb_ctl)
        print(par_tb_dta()$par_tb_trt)
        print("----------PFS----------------")
        print(userLog$par_surv)
        print("-----------------------------")
        print(par_surv_dta()$par_surv_ctl)
        print(par_surv_dta()$par_surv_trt)
    })
})
