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
        point_size = .1
        xrng = plot_zoom_xrng()
        yrng = plot_zoom_yrng()
        if(typ == GLOBAL_VIEW_POINTS){
            if(length(input$selMarks) >= 1){
                if(length(input$selMarks) == 2){
                    corr_dt = dcast(
                        agg_dt[cell %in% input$selCells & mark %in% input$selMarks], 
                        id+cell+tx+ty~mark, value.var = "value")
                    corr_dt$diff = corr_dt[, 6] - corr_dt[,5]
                    p = ggplot(corr_dt, aes(x = tx, y = ty, color = diff)) +
                        geom_point(size = point_size) + 
                        scale_color_gradientn(colors = c("blue", "white", "red"), 
                                              breaks = c(-1, -.5, 0, .5, 1),
                                              labels = c(colnames(corr_dt)[5], "", "-", "", colnames(corr_dt)[6]),
                                              limits = c(-1, 1)) +
                        coord_cartesian(xlim = xrng, ylim = yrng) +
                        labs(color = "difference")
                }else{
                    p = ggplot(
                        agg_dt[cell %in% input$selCells & mark %in% input$selMarks], 
                        aes(x = tx, y = ty, color = value)) +
                        geom_point(size = point_size) +  
                        coord_cartesian(xlim = xrng, ylim = yrng) +
                        facet_wrap("mark")    
                }
            }else{
                p = ggplot(tsne_dt[cell %in% input$selCells], aes(x = tx, y = ty)) +
                    coord_cartesian(xlim = xrng, ylim = yrng) +
                    geom_point(size = point_size)
            }
            
        }else if(typ == GLOBAL_VIEW_DENSITY){
            p = ggplot(tsne_dt[cell %in% input$selCells], aes(x = tx, y = ty)) +
                coord_cartesian(xlim = xrng, ylim = yrng) +
                geom_density2d() + 
                facet_wrap("cell") 
        }else if(typ == GLOBAL_VIEW_PROFILES_FAST){
            p = stsPlotSummaryProfiles(profile_dt,
                                       tsne_dt, 
                                       q_cells = input$selCells,
                                       q_marks = input$selMarks,
                                       x_points = input$numBins,
                                       xrng = xrng,
                                       yrng = yrng,
                                       line_color_mapping = color_mapping
                                        )
            p
        }else if(typ == GLOBAL_VIEW_PROFILES_SLOW){
            p = stsPlotSummaryProfiles(profile_dt, tsne_dt, 
                                       q_cells = input$selCells,
                                       q_marks = input$selMarks,
                                       x_points = input$numBins, 
                                       xrng = xrng,
                                       yrng = yrng,
                                       line_color_mapping = color_mapping,
                                       plot_type = "raster")
        }
        
        # if(any(input$xrng != c(-.5, .5)) | any(input$yrng != c(-.5, .5))){
        #     p = p + annotate("rect",
        #                      xmin = min(input$xrng), xmax = max(input$xrng),
        #                      ymin = min(input$yrng), ymax = max(input$yrng),
        #                      fill = "#00FF0055", color = "black")
        # }
        p +
            # coord_fixed() + 
            theme_classic() + 
            labs(x = "", y = "")
        
    })
    
    output$zoomPlot <- renderPlot({
        ggplot()
        # zimg_res = make_tsne_img(profile_dt[cell %in%  input$selCells],
        #                          tsne_dt, n_points = input$bins,
        #                          xrng = sel_zoom_xrng(), yrng = sel_zoom_yrng())
        # p = make_img_plots(zimg_res, min_size = .3, qcell = input$selCells,
        #                    xrng = sel_zoom_xrng(),
        #                    yrng = sel_zoom_yrng(),
        #                    as_facet = TRUE ) +
        #     theme_classic()
        # p
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
    
    output$globalDebug = renderText({
        if(!is.null(input$global_click) | !is.null(input$global_brush)){
            # browser()
        }
        msg_head = "globalPlot interaction:"
        if(!is.null(input$global_click)){
            click_val = lapply(input$global_click[c("x", "y")], round, digits = 3)
            msg_click = paste0("  Click:", click_val$x, ", ", click_val$y)
        }else{
            msg_click = "  Click: -"
        }
        if(!is.null(input$global_brush)){
            brush_val = lapply(input$global_brush[c("xmin", "xmax", "ymin", "ymax")], round, digits = 3)
            msg_brush = paste0("  Brush:", brush_val$xmin, " to ", brush_val$xmax, 
                               ", ", brush_val$ymin, " to ", brush_val$ymax)
        }else{
            msg_brush = "  Brush: -"
        }
        if(!is.null(sel_zoom_xrng())){
            msg_xrng = paste0("  select xrng:", paste(round(sel_zoom_xrng(), digits = 3), collapse = " to "))
        }else{
            msg_xrng = "  select xrng: -"
        }
        if(!is.null(sel_zoom_yrng())){
            msg_yrng = paste0("  select yrng:", paste(round(sel_zoom_yrng(), digits = 3), collapse = " to "))
        }else{
            msg_yrng = "  select yrng: -"
        }
        
        if(!is.null(plot_zoom_xrng())){
            msg_pxrng = paste0("  plot xrng:", paste(round(plot_zoom_xrng(), digits = 3), collapse = " to "))
        }else{
            msg_pxrng = "  plot xrng: -"
        }
        if(!is.null(plot_zoom_yrng())){
            msg_pyrng = paste0("  plot yrng:", paste(round(plot_zoom_yrng(), digits = 3), collapse = " to "))
        }else{
            msg_pyrng = "  plot yrng: -"
        }
        
        # tags$span(tags$br(msg_head), tags$br(msg_click), tags$br(msg_brush))
        paste(sep = "\n",
              msg_head, msg_click, msg_brush, 
              msg_xrng, msg_yrng, 
              msg_pxrng, msg_pyrng
              )
        
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
    
    sel_zoom_xrng = reactiveVal(c(-.5, .5))
    sel_zoom_yrng = reactiveVal(c(-.5, .5))
    
    plot_zoom_xrng = reactiveVal(c(-.5, .5))
    plot_zoom_yrng = reactiveVal(c(-.5, .5))
    
    observeEvent(input$global_brush, {
        if(is.null(input$global_brush)){
            browser()
            xrng = c(-.5, .5)
            yrng = c(-.5, .5)
        }else{
            xrng = c(input$global_brush$xmin, input$global_brush$xmax)
            yrng = c(input$global_brush$ymin, input$global_brush$ymax)    
        }
        if(any(xrng != sel_zoom_xrng())){
            sel_zoom_xrng(xrng)
        }
        if(any(yrng != sel_zoom_yrng())){
            sel_zoom_yrng(yrng)
        }
    })
    
    observeEvent(input$btnZoom, {
        plot_zoom_xrng(sel_zoom_xrng())
        plot_zoom_yrng(sel_zoom_yrng())
    })
    
    observeEvent(input$btnReset, {
        plot_zoom_xrng(c(-.5, .5))
        plot_zoom_yrng(c(-.5, .5))
    })
    
    
    output$tableGenes = DT::renderDT(expr = {
        data.frame(a = 1:3, b = 4:6)
    })
})
