#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source("setup.R")

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    shinyjs::useShinyjs(),
    
    hidden(colourInput("col1", "Col1")),#necessary for any other colourPicker to work
    # Application title
    titlePanel("seqtsne", windowTitle = "seqtsne"),
    div(
        id = "libs-content",
        h4("Loading R libraries...")
    ),
    hidden(div(
        id = "waiting-content",
        h4("Dataset has not been selected.")
    )),
    hidden(div(
        id = "loading-content",
        h4("Loading", id = "loading-status", style = "word-wrap:break-word")
        #,
        
    )),
    hidden(
        div(
            id = "app-content",
            tabsetPanel(
                tabPanel("Basic Plot",
                         sidebarLayout(
                             sidebarPanel(
                                 fluidRow(
                                     column(
                                         width = 4,
                                         radioButtons(
                                             "globalViewType",
                                             label = "Global View",
                                             choices = c(
                                                 GLOBAL_VIEW_POINTS,
                                                 GLOBAL_VIEW_POINTS_COMPARISON,
                                                 GLOBAL_VIEW_POINTS_RGB,
                                                 GLOBAL_VIEW_PROFILES_FAST,
                                                 GLOBAL_VIEW_PROFILES_SLOW#,
                                                 # GLOBAL_VIEW_DENSITY
                                             )
                                         )
                                     ),
                                     radioButtons("selFacetType", label = "Facet Type", choices = c("Free", "Cell/Mark", "Mark/Cell"), inline = TRUE),
                                     shinyjs::hidden(
                                         tags$div(id = "app-tall-wide-sel",
                                                  # shinyjs::hidden(
                                                  column(id = "app-tall-sel",
                                                         width = 4,
                                                         uiOutput("ui_global_facet_top")
                                                         # )
                                                  ),
                                                  column(id = "app-wide-sel",
                                                         width = 4,
                                                         uiOutput("ui_global_facet_bot")
                                                  )
                                         )
                                     ),
                                     shinyjs::hidden(
                                         tags$div(id = "app-compare-sel",
                                                  # shinyjs::hidden(
                                                  column(#id = "app-tall-sel",
                                                      width = 6,
                                                      uiOutput("ui_compare_marks")
                                                      # )
                                                  )
                                         )
                                     ),
                                     shinyjs::hidden(
                                         tags$div(id = "app-rgb-sel",
                                                  # shinyjs::hidden(
                                                  column(#id = "app-tall-sel",
                                                      width = 8,
                                                      uiOutput("ui_rgb_marks")
                                                      # )
                                                  )
                                         )
                                     )
                                 ),
                                 radioButtons("selGlobalColoring", 
                                              label = "Color By", 
                                              choices = c("mark", "cell", "both")),
                                 actionButton("btnCustomColors", label = "Customize Colors"),
                                 radioButtons("selNumPlotted", 
                                              label = "Number of points plotted:", 
                                              inline = TRUE,
                                              choices = c("500", "5000", "50000", "all"), 
                                              selected = "5000"),
                                 sliderInput("numBins",
                                             "Number of bins:",
                                             min = 2,
                                             max = 16,
                                             value = 8),
                                 verbatimTextOutput("globalDebug"),
                                 actionButton("btnDebug", "browser()")
                             ),
                             
                             # Show a plot of the generated distribution
                             mainPanel(
                                 fluidRow(
                                     withSpinner(plotOutput("globalPlot", 
                                                            width = "600px", 
                                                            height = "600px",
                                                            brush = brushOpts("global_brush", 
                                                                              delay = 500, 
                                                                              delayType = "debounce"), 
                                                            click = "global_click")),
                                     actionButton("btnZoom", "Zoom"),
                                     actionButton("btnReset", "Reset"),
                                     actionButton("btnDlGlobal", "Download Image")
                                 )
                             )
                         )
                ), 
                tabPanel("Gene Table",
                         sidebarLayout(
                             sidebarPanel = sidebarPanel(
                                 actionButton("dlGeneTable", label = "Download"),
                                 checkboxInput("dlUnique", label = "as unique gene list", value = TRUE),
                                 textOutput("textCountUnique")
                             ),
                             mainPanel = mainPanel(
                                 tabsetPanel(id = "tabset_gene_tables",
                                             tabPanel(title = "Genes", withSpinner(DT::dataTableOutput("tableGenes", width = "600px"))), 
                                             tabPanel(title = "GO", withSpinner(DT::dataTableOutput("tableGO", width = "600px"))))
                                 
                             )
                         )
                ),
                tabPanel("Gene Query",
                         textInput("textGeneQ", label = "Gene Query", value = "Runx1 Runx2"),
                         withSpinner(plotOutput("geneQPlot", 
                                                width = "600px", 
                                                height = "600px",
                                                brush = brushOpts("geneq_brush", 
                                                                  delay = 500, 
                                                                  delayType = "debounce"), 
                                                click = "geneq_click"))
                ),
                tabPanel("Define Sets",
                         sidebarLayout(
                             sidebarPanel = sidebarPanel(
                                 uiOutput("uiSelector"),
                                 # selectInput("selSelector", label = "Select Data Source", choices = names(DATA()))
                                 textInput("newSetName", label = "New Set Name", value = "new set"),
                                 actionButton("btnAddSet", label = "Add Set"),
                                 withSpinner(plotOutput("setPreviewGlobal", 
                                                        width = "600px", 
                                                        height = "600px")),
                                 withSpinner(plotOutput("setPreviewZoom", 
                                                        width = "600px", 
                                                        height = "600px"))
                             ),  
                             mainPanel = mainPanel(
                                 withSpinner(DT::dataTableOutput("tableSelect", width = "600px"))
                                 
                             )
                             
                         )
                )),
            tags$br(),
            sidebarLayout(
                sidebarPanel(
                    fluidRow(
                        column(
                            width = 4,
                            uiOutput("ui_zoom_cells")
                        ),
                        column(
                            width = 4,
                            uiOutput("ui_zoom_marks")
                        ),
                        column(
                            width = 4,
                            shinyjs::hidden(radioButtons("detail_file_type", "File Type", choices = c("bigwig", "bam"))),
                            shinyjs::hidden(radioButtons("detail_bam_plot_ype", "Plot Type", choices = c("stranded", "pileup", "scc", "shiftmap"))),
                            numericInput("n_detail", "Number of regions to plot", 
                                         value = 5, min = 1, max = 20, step = 1),
                            numericInput("detail_view_size", "View size", 
                                         value = 1000, min = 100, max = 2e4, step = 100),
                            radioButtons("sel_detail_type", label = "Select Type", choices = c("sample", "aggregate"))
                        )
                    )
                ),
                mainPanel(
                    fluidRow(
                        column(
                            width = 4,
                            withSpinner(plotOutput("detailPlot", 
                                                   width = "600px", 
                                                   height = "600px")),
                            actionButton("btnDlDetail", "Download Image")
                        )
                    )
                    
                )
            )
        ))
))