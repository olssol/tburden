##    ----------------------------------------------------------------------
##    Copyright (C) 2015  Daniel O. Scharfstein and Chenguang Wang
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.

##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##    -----------------------------------------------------------------------

options(shiny.maxRequestSize = 200*1024^2)

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
            select(SUBJID, ARM, PFS_EVENT, PFS_DAYS, OS_EVENT, OS_DAYS)
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
            select(SUBJID, ARM, AGE, SEX, STRATA1, P1TERTL)
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

        tb_plt_tb(dat_tb, id)
    })

    output$txtHist <- renderPrint({
        print(get_cur_hist())
    })

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
