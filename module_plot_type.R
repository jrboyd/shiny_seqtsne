server_plot_type = function(input,
                            output,
                            session,
                            DATA,
                            UI_TALLV,
                            UI_WIDEV,
                            plot_zoom_xrng,
                            plot_zoom_yrng,
                            get_curr_col,
                            GLOBAL_PLOT) {
    
    
    output$globalPlot <- renderPlot({
        req(input$selTallVars)
        req(input$selWideVars)
        typ = input$globalViewType
        
        xrng = plot_zoom_xrng()
        yrng = plot_zoom_yrng()
        frac_shown = max(c(diff(xrng), diff(yrng)))
        point_size = .1 / frac_shown
        in_view_id =  DATA()$tsne_dt[
            tx >= min(xrng) & 
                tx <= max(xrng) & 
                ty >= min(yrng) & 
                ty <= max(yrng)]$id
        if(!input$selNumPlotted == "all"){
            set.seed(0)
            in_view_id = sampleCap(in_view_id, as.numeric(input$selNumPlotted))
        }
        if(typ == GLOBAL_VIEW_POINTS){
            if(length(input$selWideVars) >= 1){
                nc = ceiling(sqrt(length(input$selWideVars)))
                p = ggplot(
                    DATA()$agg_dt[id %in% in_view_id][tall_var %in% input$selTallVars & wide_var %in% input$selWideVars], 
                    aes(x = tx, y = ty, color = value)) +
                    geom_point(size = point_size) +  
                    scale_color_viridis_c() +
                    coord_fixed(xlim = xrng, ylim = yrng, ratio = diff(xrng)/diff(yrng)) +
                    facet_wrap("wide_var", ncol = nc)    
                # }
            }else{
                p = ggplot(DATA()$tsne_dt[id %in% in_view_id][tall_var %in% input$selTallVars], aes(x = tx, y = ty)) +
                    coord_fixed(xlim = xrng, ylim = yrng, ratio = diff(xrng)/diff(yrng)) +
                    geom_point(size = point_size)
            }
        }else if(typ == GLOBAL_VIEW_POINTS_COMPARISON){
            req(input$selCompare1)
            req(input$selCompare2)
            p = plot_agree_diff(
                p_dt = DATA()$agg_dt,
                sel_id = in_view_id,
                sel_tall = input$selTallVars,
                sel_wide1 = input$selCompare1,
                sel_wide2 = input$selCompare2,
                point_size = point_size,
                xrng = xrng,
                yrng = yrng
            )
        }else if(typ == GLOBAL_VIEW_POINTS_RGB){
            req(input$selCompareR)
            req(input$selCompareG)
            req(input$selCompareB)
            p = plot_rgb(
                p_dt = DATA()$agg_dt,
                sel_id = in_view_id,
                sel_tall = input$selTallVars,
                sel_wideR = input$selCompareR,
                sel_wideG = input$selCompareG,
                sel_wideB = input$selCompareB,
                point_size = point_size,
                xrng = xrng,
                yrng = yrng
            )
        }else if(typ == GLOBAL_VIEW_DENSITY){
            p = ggplot(DATA()$tsne_dt[tall_var %in% input$selTallVars], aes(x = tx, y = ty)) +
                coord_cartesian(xlim = xrng, ylim = yrng) +
                geom_density2d() + 
                facet_wrap("tall_var") 
        }else if(typ == GLOBAL_VIEW_PROFILES_FAST){
            cm = get_curr_col()
            p = stsPlotSummaryProfiles(DATA()$profile_dt[id %in% in_view_id, ],
                                       DATA()$tsne_dt[id %in% in_view_id, ], 
                                       q_tall_vars = input$selTallVars,
                                       q_wide_vars = input$selWideVars,
                                       x_points = input$numBins, 
                                       min_size = 1,
                                       N_ceiling = 1,
                                       xrng = xrng,
                                       yrng = yrng,
                                       apply_norm = FALSE,
                                       ylim = DATA()$tylim,
                                       line_color_mapping = cm
            )
            p
        }else if(typ == GLOBAL_VIEW_PROFILES_SLOW){
            cm = get_curr_col()
            p = stsPlotSummaryProfiles(DATA()$profile_dt[id %in% in_view_id, ], 
                                       DATA()$tsne_dt[id %in% in_view_id, ], 
                                       q_tall_vars = input$selTallVars,
                                       q_wide_vars = input$selWideVars,
                                       x_points = input$numBins, 
                                       min_size = 1,
                                       N_ceiling = 1,
                                       xrng = xrng,
                                       yrng = yrng,
                                       apply_norm = FALSE,
                                       ylim = DATA()$tylim,
                                       line_color_mapping = cm,
                                       plot_type = "raster")
        }else{
            stop("unrecognized input$globalViewType ", typ)
        }
        p = p +
            theme_classic() + 
            labs(x = "", y = "")
        GLOBAL_PLOT(p)
        p
        
    })
    
    observe({
        req(UI_TALLV())
        req(input$selTallVars)
        if(length(UI_TALLV()) == 1){
            # showNotification("hide it")
            shinyjs::hide("app-tall-sel")
        }
    })
    
    
    output$ui_global_facet_top = renderUI({
        req(UI_TALLV())
        checkboxGroupInput(
            "selTallVars",
            "Select Aspect",
            choices = UI_TALLV(),
            selected = UI_TALLV()
        )
    })
    
    output$ui_global_facet_bot = renderUI({
        req(UI_WIDEV())
        checkboxGroupInput(
            "selWideVars",
            "Select Variables",
            choices = UI_WIDEV(),
            selected = UI_WIDEV()[1]
        )
    })
    
    output$ui_compare_marks = renderUI({ #TODO 2 and exactly 2 marks
        req(UI_WIDEV())
        fluidRow(
            column(width = 6, 
                   radioButtons(
                       "selCompare1",
                       "Select 1",
                       choices = UI_WIDEV(),
                       selected = UI_WIDEV()[1]
                   )),
            column(width = 6,
                   radioButtons(
                       "selCompare2",
                       "Select 2",
                       choices = UI_WIDEV(),
                       selected = UI_WIDEV()[2]
                   )
            )
        )
    })
    
    output$ui_rgb_marks = renderUI({ #TODO at least 1 up to 3
        req(UI_WIDEV())
        wide_v = c(UI_WIDEV(), "none")
        qr = min(length(wide_v), 1)
        qg = min(length(wide_v), 2)
        qb = min(length(wide_v), 3)
        fluidRow(
            column(width = 4, 
                   radioButtons(
                       "selCompareR",
                       "Select R",
                       choices = wide_v,
                       selected = wide_v[qr]
                   )),
            column(width = 4,
                   radioButtons(
                       "selCompareG",
                       "Select G",
                       choices = wide_v,
                       selected = wide_v[qg]
                   )),
            column(width = 4,
                   radioButtons(
                       "selCompareB",
                       "Select B",
                       choices = wide_v,
                       selected = wide_v[qb]
                   ))
            
        )
    })
    
    observeEvent({#manage conditional UI elements based on globalViewType
        input$globalViewType
    }, {
        if(input$globalViewType %in% c(GLOBAL_VIEW_POINTS, 
                                       GLOBAL_VIEW_PROFILES_FAST, 
                                       GLOBAL_VIEW_PROFILES_SLOW)){
            shinyjs::show("app-tall-wide-sel")
            shinyjs::hide("app-compare-sel")
            shinyjs::hide("app-rgb-sel")
        }else if(input$globalViewType == GLOBAL_VIEW_POINTS_COMPARISON){
            shinyjs::hide("app-tall-wide-sel")
            shinyjs::show("app-compare-sel")
            shinyjs::hide("app-rgb-sel")
        }else  if(input$globalViewType == GLOBAL_VIEW_POINTS_RGB){
            shinyjs::hide("app-tall-wide-sel")
            shinyjs::hide("app-compare-sel")
            shinyjs::show("app-rgb-sel")
        }else{
            stop("unrecognized input$globalViewType")
        }
    })
    
    
    
}

plot_rgb = function(p_dt,
                    sel_id,
                    sel_tall,
                    sel_wideR,
                    sel_wideG,
                    sel_wideB,
                    point_size,
                    xrng,
                    yrng
) {
    p_dt = p_dt[id %in% sel_id]
    sel_wide = unique(c(sel_wideR, sel_wideG, sel_wideB))
    sel_wide = sel_wide[sel_wide != "none"]
    p_dt = p_dt[tall_var %in% sel_tall &
                    wide_var %in% sel_wide]
    compare_dt = dcast(p_dt, id + tall_var + tx + ty ~ wide_var, value.var = "value")
    compare_dt$none = 0
    compare_dt[, color := rgb(get(sel_wideR), get(sel_wideG), get(sel_wideB))]
    p = ggplot() +
        geom_point(
            data = compare_dt,
            aes(x = tx, y = ty, color = color),
            size = point_size
        ) +
        scale_color_identity() +
        coord_fixed(xlim = xrng,
                    ylim = yrng,
                    ratio = diff(xrng) / diff(yrng)) +
        facet_wrap("tall_var") 
}

plot_agree_diff = function(p_dt,
                           sel_id,
                           sel_tall,
                           sel_wide1,
                           sel_wide2,
                           point_size,
                           xrng,
                           yrng) {
    compare_dt = dcast(p_dt[id %in% sel_id][tall_var %in% sel_tall &
                                             wide_var %in% c(sel_wide1, sel_wide2)],
                    id + tall_var + tx + ty ~ wide_var, value.var = "value")
    compare_dt$difference = compare_dt[[sel_wide2]] - compare_dt[[sel_wide1]]
    compare_dt$agreement = pmin(compare_dt[[sel_wide2]], compare_dt[[sel_wide1]])
    val_dt = melt(
        compare_dt[, .(tx, ty, difference, agreement)],
        id.vars = c("tx", "ty"),
        measure.vars = c("difference", "agreement")
    )
    p = ggplot() +
        geom_point(
            data = val_dt[variable == "difference"],
            aes(x = tx, y = ty, color = value),
            size = point_size) +
        geom_point(
            data = val_dt[variable == "agreement"],
            aes(x = tx, y = ty, fill = value),
            size = point_size * 1.5,
            shape = 21,
            color = "#00000000"
        ) +
        scale_color_gradientn(
            colors = c("blue", "white", "red"),
            breaks = c(-1, -.5, 0, .5, 1),
            labels = c(sel_wide1, "", "-", "", sel_wide2),
            limits = c(-1, 1)
        ) +
        scale_fill_gradientn(
            colors = c("gray80", "gray0"),
            breaks = c(0, .5, 1),
            labels = c(0, .5, 1),
            limits = c(0, 1)
        ) +
        coord_fixed(xlim = xrng,
                    ylim = yrng,
                    ratio = diff(xrng) / diff(yrng)) +
        facet_wrap("variable") +
        labs(color = "difference",
             fill = paste("min of", sel_wide1,
                          "\nand", sel_wide2))
}