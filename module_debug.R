server_debug = function(input, output, session, 
                        plot_zoom_xrng, plot_zoom_yrng, 
                        sel_zoom_xrng, sel_zoom_yrng){
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
    

}