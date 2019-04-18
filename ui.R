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

    # Application title
    titlePanel("Runx1 bookmarking with seqtsne", windowTitle = "seqtsne"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            radioButtons("globalViewType", label = "Global View",
                         choices = c(GLOBAL_VIEW_POINTS,
                                     GLOBAL_VIEW_PROFILES_FAST,
                                     GLOBAL_VIEW_PROFILES_SLOW,
                                     GLOBAL_VIEW_DENSITY)),

            sliderInput("bins",
                        "Number of bins:",
                        min = 2,
                        max = 16,
                        value = 8),
            # sliderInput("xrng", "x-range", -.5, .5, value = c(-.5, .5), dragRange = TRUE),
            # sliderInput("yrng", "y-range", -.5, .5, value = c(-.5, .5), dragRange = TRUE),
            # actionButton("doZoom", label = "Zoom"),
            selectInput("selCells", "Select Cells", choices = UI_CELLS, selected = UI_CELLS[1:4], multiple = TRUE),
            selectInput("selGenes", "Select Genes", choices = UI_GENES, selected = "RUNX1")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(
                plotOutput("globalPlot", width = "400px", height = "400px") ,
                plotOutput("zoomPlot", width = "400px", height = "400px")
            ),
            fluidRow(
                plotOutput("genePlot", width = "400px", height = "400px"),
                plotOutput("profilePlot", width = "400px", height = "400px")
            ),
            fluidRow(
                plotOutput("pairArrows", width = "400px", height = "400px"),
                plotOutput("pairKey", width = "400px", height = "400px")

            )
        )
    )
))
