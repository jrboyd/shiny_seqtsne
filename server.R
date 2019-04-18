#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    output$globalPlot <- renderPlot({
        typ = input$globalViewType
        if(typ == GLOBAL_VIEW_POINTS){
            p = ggplot(tsne_tp, aes(x = tx, y = ty)) +
                geom_point()
        }else if(typ == GLOBAL_VIEW_DENSITY){
            p = ggplot(tsne_dt, aes(x = tx, y = ty)) +
                geom_density2d()
        }else if(typ == GLOBAL_VIEW_PROFILES_FAST){
            p = stsPlotSummaryProfiles(profile_dt, tsne_dt, x_points = 8)
        }else if(typ == GLOBAL_VIEW_PROFILES_SLOW){
            # n_points = 10
            # piles_img_res = make_tsne_img(profiles_dt = profile_dt,
            #                               position_dt = tsne_dt, #force_rewrite = TRUE,
            #                               apply_norm = FALSE,
            #                               ylim = c(0,10),
            #                               # xrng = zoom_x,
            #                               # yrng = zoom_y,
            #                               n_points = n_points,
            #                               line_colors = c(
            #                                   "H3K4me3" = "forestgreen",
            #                                   "H3K27me3" = "firebrick")
            # )
            # 
            # p_basic = make_img_plots(img_results = list(piles_img_res),
            #                          qcell = NULL,
            #                          min_size = 0,
            #                          N_ceiling = NULL,
            #                          as_facet = FALSE)[[1]]
            # p = p_basic
            p = stsPlotSummaryProfiles(profile_dt, tsne_dt, x_points = 8, plot_type = "raster")
        }

        if(any(input$xrng != c(-.5, .5)) | any(input$yrng != c(-.5, .5))){
            p = p + annotate("rect",
                             xmin = min(input$xrng), xmax = max(input$xrng),
                             ymin = min(input$yrng), ymax = max(input$yrng),
                             fill = "#00FF0055", color = "black")
        }
        p +
            coord_fixed() + theme_classic() + labs(x = "", y = "")

    })

    output$zoomPlot <- renderPlot({
        zimg_res = make_tsne_img(profile_dt[cell %in%  input$selCells],
                                 tsne_dt, n_points = input$bins,
                                 xrng = zoom_xrng(), yrng = zoom_yrng())
        p = make_img_plots(zimg_res, min_size = .3, qcell = input$selCells,
                           xrng = zoom_xrng(),
                           yrng = zoom_yrng(),
                           as_facet = TRUE ) +
            theme_classic()
        p
    })

    output$imgPlot <- renderPlot({
        p_basic +
            coord_fixed() + theme_classic() + labs(x = "", y = "")
    })

    output$genePlot <- renderPlot({
        plot_velocity_arrows_selected(tsne_dt,
                                      query_gr,
                                      input$selCells,
                                      tss_ids = input$selGenes) +
            coord_fixed() + theme_classic() + labs(x = "", y = "")
    })

    output$profilePlot <- renderPlot({
        plot_profiles_selected(profile_dt,
                               query_gr,
                               input$selCells,
                               tss_ids = input$selGenes) +
            theme_classic() + labs(x = "", y = "")
    })

    cellPair = reactiveVal(NULL)
    observeEvent({
        input$selCells
    }, {
        if(length(input$selCells) > 1){
            sc = input$selCells[1:2]
            if(!all(sc == cellPair())){
                cellPair(sc)
            }
        }

    })


    output$pairArrows = renderPlot({
        # req(cellPair())
        # cp = cellPair()
        cp = c("H7", "CD34")
        if(length(cp) == 2){
            vel_plots = plot_velocity_arrows(tsne_dt, cp[1], cp[2])
        }
        vel_plots[[1]] + theme_classic()
    })

    output$pairKey = renderPlot({
        # req(cellPair())
        # cp = cellPair()
        cp = c("H7", "CD34")
        if(length(cp) == 2){
            vel_plots = plot_velocity_arrows(tsne_dt, cp[1], cp[2])
        }
        vel_plots[[2]]
    })

    # output$pairArrows = renderPlot({
    #
    #         vel_pl
    #     plot_velocity_arrows(tsne_dt, cp[1], cp[2])[[1]] + theme_classic()
    # })

    zoom_xrng = reactiveVal(c(-.5, .5))
    zoom_yrng = reactiveVal(c(-.5, .5))

    observeEvent(input$doZoom, {
        req(input$xrng)
        req(input$yrng)
        if(any(input$xrng != zoom_xrng())){
            zoom_xrng(input$xrng)
        }
        if(any(input$yrng != zoom_yrng())){
            zoom_yrng(input$yrng)
        }
    })

})
