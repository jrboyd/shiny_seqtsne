#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
source("setup.R")

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    selectInput("selDataset", "Select Dataset", choices = UI_DATASETS),
    # Application title
    titlePanel("Runx1 bookmarking with seqtsne", windowTitle = "seqtsne"),
    
    # Sidebar with a slider input for number of bins
    
    tabsetPanel(tabPanel("Basic Plot",
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
                                                 GLOBAL_VIEW_PROFILES_FAST,
                                                 GLOBAL_VIEW_PROFILES_SLOW,
                                                 GLOBAL_VIEW_DENSITY
                                             )
                                         )
                                     ),
                                     column(
                                         width = 4,
                                         uiOutput("ui_global_cells")
                                         # checkboxGroupInput(
                                         #     "selCells",
                                         #     "Select Cells",
                                         #     choices = UI_CELLS,
                                         #     selected = UI_CELLS
                                         # )
                                     ),
                                     column(
                                         width = 4,
                                         uiOutput("ui_global_marks")
                                         # checkboxGroupInput(
                                         #     "selMarks",
                                         #     "Select Marks",
                                         #     choices = UI_MARKS,
                                         #     selected = UI_MARKS[1]
                                         # )
                                     )
                                 ),
                                 
                                 sliderInput("numBins",
                                             "Number of bins:",
                                             min = 2,
                                             max = 16,
                                             value = 8),
                                 # selectInput("selGenes", "Select Genes", choices = UI_GENES, selected = "RUNX1"),
                                 verbatimTextOutput("globalDebug")
                             ),
                             
                             # Show a plot of the generated distribution
                             mainPanel(
                                 fluidRow(
                                     plotOutput("globalPlot", width = "600px", height = "600px",
                                                brush = brushOpts("global_brush", delay = 500, delayType = "debounce"), click = "global_click"),
                                     actionButton("btnZoom", "Zoom"),
                                     actionButton("btnReset", "Reset")
                                 )#,
                                 # fluidRow(
                                 #     plotOutput("genePlot", width = "400px", height = "400px"),
                                 #     plotOutput("profilePlot", width = "400px", height = "400px")
                                 # ),
                                 # fluidRow(
                                 #     plotOutput("pairArrows", width = "400px", height = "400px"),
                                 #     plotOutput("pairKey", width = "400px", height = "400px")
                                 # 
                                 # )
                             )
                         )
    ), 
    tabPanel("Gene Table",
             sidebarLayout(
                 sidebarPanel = sidebarPanel(),
                 mainPanel = mainPanel(
                     DT::dataTableOutput("tableGenes", width = "600px")
                 )
             )
    )
    ,
    tabPanel("Gene Query",
             textInput("textGeneQ", label = "Gene Query", value = "Runx1 Runx2"),
             plotOutput("geneQPlot", width = "600px", height = "600px")
    )
    ),
    sidebarLayout(
        sidebarPanel(
            fluidRow(
                column(
                    width = 4,
                    uiOutput("ui_zoom_cells")
                    # checkboxGroupInput(
                    #     "selCellsDetail",
                    #     "Select Cells",
                    #     choices = UI_CELLS,
                    #     selected = UI_CELLS
                    # )
                ),
                column(
                    width = 4,
                    uiOutput("ui_zoom_marks")
                    # checkboxGroupInput(
                    #     "selMarksDetail",
                    #     "Select Marks",
                    #     choices = UI_MARKS,
                    #     selected = UI_MARKS
                    # )
                ),
                column(
                    width = 4,
                    numericInput("n_detail", "Number of regions to plot", value = 5, min = 1, max = 20, step = 1)  
                )
            )
        ),
        mainPanel(
            plotOutput("detailPlot", width = "600px", height = "600px")
        )
    )
)
)
