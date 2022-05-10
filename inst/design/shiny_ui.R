
## ------------------------------------------------------------------
##
##                     TUMOR BURDEN
##
## ------------------------------------------------------------------
box_tb_input <-
    shinydashboard::box(title  = "Tumor Burden Input",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width       = 12,
                        height      = 1150,
                        textInput("inVecT",
                                  label = "Time Points (Days)",
                                  value = "9,18,27,36,45,54,63"),
                        textInput("inMeanCtl",
                                  label = "Mean PCHG: Control",
                                  value = "-0.2, -0.2,  -0.25, -0.2, -0.15, -0.1, -0.05"),
                        textInput("inMeanTrt",
                                  label = "Mean PCHG: Treatment",
                                  value = "-0.2, -0.25,  -0.3, -0.3, -0.2, -0.2, -0.1"),
                        textInput("inTbSig",
                                  label = "SD PCHG",
                                  value = "0.1"),
                        textInput("inTbRndCtl",
                                  label = "Random: Control",
                                  value = "0.2"),
                        textInput("inTbRndTrt",
                                  label = "Random: Treatment",
                                  value = "0.2"),
                        textInput("inLblCtl",
                                  label = "Label: Control",
                                  value = "Control"),
                        textInput("inLblTrt",
                                  label = "Label: Treatment",
                                  value = "Treatment"),
                        textInput("inN",
                                  label = "Number of subjects",
                                  value = 100),
                        actionButton("btnUpdateTb",
                                     label = "Update Tumor Burden Results",
                                     class = "btn-success")
                        )

box_tb_output <-
    shinydashboard::box(title  = "Tumor Burden Output",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width       = 12,
                        height      = 1150,
                        fluidRow(
                            column(8,
                                   wellPanel(
                                       h5("Mean Tumor Burden"),
                                       plotOutput("pltMeanTb", height = "275px")
                                   ),
                                   wellPanel(
                                       h5("Individual Tumor Burden"),
                                       plotOutput("pltSubTb", height = "275px")),
                                   wellPanel(
                                       h5("Best PCHG of Responders"),
                                       plotOutput("pltSubTbWf", height = "275px"))
                                   ),
                            column(4,
                                   wellPanel(
                                       fluidRow(
                                           column(12,
                                                  h5("Pseudo-Response Rate"),
                                                  tableOutput("tblResp"),
                                                  h5("Tumor Burden PCHG"),
                                                  tableOutput("tblTb")
                                                  )
                                           ## ,column(0,
                                           ##        h5("Pseudo-Response Rate by Random Effect"),
                                           ##        tableOutput("tblRespPos")
                                           ##        )
                                       )
                                   ))
                        )
)


## ------------------------------------------------------------------
##
##                     SURVIVAL
##
## ------------------------------------------------------------------

box_surv_input <-
    shinydashboard::box(title  = "PFS Input",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width       = 12,
                        height      = 850,
                        textInput("inMedianCtl",
                                  label = "Median PFS (Months): Control",
                                  value = "12"),
                        textInput("inMedianTrt",
                                  label = "Median PFS (Months): Treatment",
                                  value = "16"),
                        textInput("inSurvRndCtl",
                                  label = "Random (Association): Control",
                                  value = "0.1"),
                        textInput("inSurvRndTrt",
                                  label = "Random (Association): Treatment",
                                  value = "0.1"),
                        textInput("inSurvXlim",
                                  label = "Limit of Time (Months)",
                                  value = "54"),
                        actionButton("btnUpdateSurv",
                                     label = "Update Survival Results",
                                     class = "btn-success")
                        )

box_surv_output <-
    shinydashboard::box(title  = "PFS Output",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width       = 12,
                        height      = 850,
                        fluidRow(
                            column(8,
                                   wellPanel(
                                       h5("PFS by ARM"),
                                       plotOutput("pltSurvArm",
                                                  height = "300px")
                                   ),
                                   wellPanel(
                                       h5("PFS by ARM and Response"),
                                       plotOutput("pltSurvArmResp",
                                                  height = "300px"))
                                   ),
                            column(4,
                                   ## wellPanel(
                                       ## h5("PFS by ARM and Random Effect"),
                                       ## plotOutput("pltSurvArmPos",
                                       ##           height = "300px")),
                                   wellPanel(
                                       fluidRow(
                                           column(12,
                                                  h5("PFS by ARM"),
                                                  tableOutput("tblSurvArm"),
                                                  h5("PFS by ARM and Response"),
                                                  tableOutput("tblSurvArmResp"))
                                           ## ,column(6,
                                           ##        h5("PFS by ARM and Random Effect"),
                                           ##        tableOutput("tblSurvArmPos"))
                                       )
                                   ))
                        ))

box_info <-
    shinydashboard::box(title  = "Info",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width       = 12,
                        verbatimTextOutput("txtInfo"))

## ------------------------------------------------------------------
##
##                     MAJOR FRAME
##
## ------------------------------------------------------------------

## Header
header <- shinydashboard::dashboardHeader(
                              title = "Study Design Scenarios for Applying TB-Integrated Survival Analysis",
                              titleWidth = 800
                          )

## Sidebar
sidebar <- shinydashboard::dashboardSidebar(disable = TRUE)

## Body
body <- shinydashboard::dashboardBody(
    tags$style(".shiny-file-input-progress {display: none}"),
    tags$head(tags$link(
                       rel = "stylesheet",
                       type = "text/css", href = "styles.css"
                   )
              ),
    height = 2500,
    fluidRow(
        column(3, box_tb_input),
        column(9, box_tb_output)),
    fluidRow(
        column(3, box_surv_input),
        column(9, box_surv_output))
    ## ,box_info
    )
