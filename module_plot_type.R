server_plot_type = function(input, output, session, DATA, UI_TALLV, UI_WIDEV, plot_zoom_xrng, plot_zoom_yrng){
    
    output$globalPlot <- renderPlot({
        req(input$selCells)
        req(input$selMarks)
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
        # tp = 
        if(typ == GLOBAL_VIEW_POINTS){
            if(length(input$selMarks) >= 1){
                #     if(length(input$selMarks) == 2){
                #         corr_dt = dcast(
                #             DATA()$agg_dt[id %in% in_view_id][tall_var %in% input$selCells & wide_var %in% input$selMarks],
                #             id+tall_var+tx+ty~wide_var, value.var = "value")
                #         corr_dt$difference = corr_dt[, 6] - corr_dt[,5]
                #         corr_dt$agreement = pmin(corr_dt[, 6], corr_dt[,5])
                #         val_dt = melt(corr_dt[, .(tx, ty, difference, agreement)],
                #                       id.vars = c("tx", "ty"),
                #                       measure.vars = c("difference", "agreement"))
                #         p = ggplot() +
                #             geom_point(data = val_dt[variable == "difference"],
                #                        aes(x = tx, y = ty, color = value),
                #                        size = point_size) +
                #             geom_point(
                #                 data = val_dt[variable == "agreement"],
                #                 aes(x = tx, y = ty, fill = value),
                #                 size = point_size * 1.5,
                #                 shape = 21,
                #                 color = "#00000000"
                #             ) +
                #             scale_color_gradientn(
                #                 colors = c("blue", "white", "red"),
                #                 breaks = c(-1,-.5, 0, .5, 1),
                #                 labels = c(colnames(corr_dt)[5], "", "-", "", colnames(corr_dt)[6]),
                #                 limits = c(-1, 1)
                #             ) +
                #             scale_fill_gradientn(
                #                 colors = c("gray80", "gray0"),
                #                 breaks = c(0, .5, 1),
                #                 labels = c(0, .5, 1),
                #                 limits = c(0, 1)
                #             ) +
                #             coord_fixed(xlim = xrng, ylim = yrng, ratio = diff(xrng)/diff(yrng)) +
                #             facet_wrap("variable") +
                #             labs(color = "difference",
                #                  fill = paste("min of", colnames(corr_dt)[5],
                #                               "\nand", colnames(corr_dt)[6]))
                #     }else{
                nc = ceiling(sqrt(length(input$selMarks)))
                p = ggplot(
                    DATA()$agg_dt[id %in% in_view_id][tall_var %in% input$selCells & wide_var %in% input$selMarks], 
                    aes(x = tx, y = ty, color = value)) +
                    geom_point(size = point_size) +  
                    # coord_cartesian(xlim = xrng, ylim = yrng) +
                    coord_fixed(xlim = xrng, ylim = yrng, ratio = diff(xrng)/diff(yrng)) +
                    facet_wrap("wide_var", ncol = nc)    
                # }
            }else{
                p = ggplot(DATA()$tsne_dt[id %in% in_view_id][tall_var %in% input$selCells], aes(x = tx, y = ty)) +
                    # coord_cartesian(xlim = xrng, ylim = yrng) +
                    coord_fixed(xlim = xrng, ylim = yrng, ratio = diff(xrng)/diff(yrng)) +
                    geom_point(size = point_size)
            }
        }else if(typ == GLOBAL_VIEW_POINTS_COMPARISON){
            plot_compare(DATA()$agg_dt, in_view_id, input$selCells, input$selCompare)
            # corr_dt = dcast(
            #     DATA()$agg_dt[id %in% in_view_id][tall_var %in% input$selCells & wide_var %in% input$selCompare],
            #     id+tall_var+tx+ty~wide_var, value.var = "value")
            # corr_dt$difference = corr_dt[, 6] - corr_dt[,5]
            # corr_dt$agreement = pmin(corr_dt[, 6], corr_dt[,5])
            # val_dt = melt(corr_dt[, .(tx, ty, difference, agreement)],
            #               id.vars = c("tx", "ty"),
            #               measure.vars = c("difference", "agreement"))
            # p = ggplot() +
            #     geom_point(data = val_dt[variable == "difference"],
            #                aes(x = tx, y = ty, color = value),
            #                size = point_size) +
            #     geom_point(
            #         data = val_dt[variable == "agreement"],
            #         aes(x = tx, y = ty, fill = value),
            #         size = point_size * 1.5,
            #         shape = 21,
            #         color = "#00000000"
            #     ) +
            #     scale_color_gradientn(
            #         colors = c("blue", "white", "red"),
            #         breaks = c(-1,-.5, 0, .5, 1),
            #         labels = c(colnames(corr_dt)[5], "", "-", "", colnames(corr_dt)[6]),
            #         limits = c(-1, 1)
            #     ) +
            #     scale_fill_gradientn(
            #         colors = c("gray80", "gray0"),
            #         breaks = c(0, .5, 1),
            #         labels = c(0, .5, 1),
            #         limits = c(0, 1)
            #     ) +
            #     coord_fixed(xlim = xrng, ylim = yrng, ratio = diff(xrng)/diff(yrng)) +
            #     facet_wrap("variable") +
            #     labs(color = "difference",
            #          fill = paste("min of", colnames(corr_dt)[5],
            #                       "\nand", colnames(corr_dt)[6]))
        }else if(typ == GLOBAL_VIEW_DENSITY){
            p = ggplot(DATA()$tsne_dt[tall_var %in% input$selCells], aes(x = tx, y = ty)) +
                coord_cartesian(xlim = xrng, ylim = yrng) +
                geom_density2d() + 
                facet_wrap("tall_var") 
        }else if(typ == GLOBAL_VIEW_PROFILES_FAST){
            cm = get_curr_col()
            p = stsPlotSummaryProfiles(DATA()$profile_dt[id %in% in_view_id, ],
                                       DATA()$tsne_dt[id %in% in_view_id, ], 
                                       q_tall_vars = input$selCells,
                                       q_wide_vars = input$selMarks,
                                       x_points = input$numBins,
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
                                       q_tall_vars = input$selCells,
                                       q_wide_vars = input$selMarks,
                                       x_points = input$numBins, 
                                       xrng = xrng,
                                       yrng = yrng,
                                       apply_norm = FALSE,
                                       ylim = DATA()$tylim,
                                       line_color_mapping = cm,
                                       plot_type = "raster")
        }
        p +
            theme_classic() + 
            labs(x = "", y = "")
        
    })
    
    observe({
        req(UI_TALLV())
        req(input$selCells)
        if(length(UI_TALLV()) == 1){
            # showNotification("hide it")
            shinyjs::hide("app-tall-sel")
        }
    })
    
    
    output$ui_global_cells = renderUI({
        req(UI_TALLV())
        checkboxGroupInput(
            "selCells",
            "Select Aspect",
            choices = UI_TALLV(),
            selected = UI_TALLV()
        )
    })
    
    output$ui_global_marks = renderUI({
        req(UI_WIDEV())
        checkboxGroupInput(
            "selMarks",
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
        fluidRow(
            column(width = 4, 
                   radioButtons(
                       "selCompareR",
                       "Select R",
                       choices = UI_WIDEV(),
                       selected = UI_WIDEV()[1]
                   )),
            column(width = 4,
                   radioButtons(
                       "selCompareG",
                       "Select G",
                       choices = UI_WIDEV(),
                       selected = UI_WIDEV()[2]
                   )),
            column(width = 4,
                   radioButtons(
                       "selCompareB",
                       "Select B",
                       choices = UI_WIDEV(),
                       selected = UI_WIDEV()[3]
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

plot_compare = function(p_dt, sel_id, sel_tall, sel_wide){
    corr_dt = dcast(
        p_dt[id %in% sel_id][tall_var %in% sel_tall & wide_var %in% sel_wide],
        id+tall_var+tx+ty~wide_var, value.var = "value")
    corr_dt$difference = corr_dt[, 6] - corr_dt[,5]
    corr_dt$agreement = pmin(corr_dt[, 6], corr_dt[,5])
    val_dt = melt(corr_dt[, .(tx, ty, difference, agreement)],
                  id.vars = c("tx", "ty"),
                  measure.vars = c("difference", "agreement"))
    p = ggplot() +
        geom_point(data = val_dt[variable == "difference"],
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
            breaks = c(-1,-.5, 0, .5, 1),
            labels = c(colnames(corr_dt)[5], "", "-", "", colnames(corr_dt)[6]),
            limits = c(-1, 1)
        ) +
        scale_fill_gradientn(
            colors = c("gray80", "gray0"),
            breaks = c(0, .5, 1),
            labels = c(0, .5, 1),
            limits = c(0, 1)
        ) +
        coord_fixed(xlim = xrng, ylim = yrng, ratio = diff(xrng)/diff(yrng)) +
        facet_wrap("variable") +
        labs(color = "difference",
             fill = paste("min of", colnames(corr_dt)[5],
                          "\nand", colnames(corr_dt)[6]))
}