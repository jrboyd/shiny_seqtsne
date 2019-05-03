#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    #stop app when tab closed
    session$onSessionEnded(stopApp)
    #load libs with a progress bar.
    load_libs_withProgress(
        libs = c(
            "GenomicRanges", 
            "data.table",
            "seqsetvis",
            "seqtsne"
        ),
        session = session
    )
    
    UI_TALLV = reactiveVal(NULL)
    UI_WIDEV = reactiveVal(NULL)
    UI_CELLS = reactiveVal(NULL)
    UI_MARKS = reactiveVal(NULL)
    UI_GENES = reactiveVal(NULL)
    DATA = reactiveVal(NULL)
    
    server_loading(input, output, session, 
                   DATA, 
                   UI_TALLV,
                   UI_WIDEV,
                   UI_CELLS,
                   UI_MARKS,
                   UI_GENES)
    custom_colors = server_custom_colors(input, ouptut, session, DATA)
    server_plot_type(input, output, session, 
                     DATA,
                     UI_TALLV, UI_WIDEV,
                     plot_zoom_xrng, plot_zoom_yrng, 
                     get_curr_col = custom_colors$get_curr_col)
    server_debug(input, output, session, 
                 plot_zoom_xrng, plot_zoom_yrng,
                 sel_zoom_xrng, sel_zoom_yrng)
    
    output$ui_zoom_cells = renderUI({
        req(UI_CELLS())
        checkboxGroupInput(
            "selCellsDetail",
            "Select Cells",
            choices = UI_CELLS(),
            selected = UI_CELLS()
        )
    })
    
    output$ui_zoom_marks = renderUI({
        req(UI_MARKS())
        checkboxGroupInput(
            "selMarksDetail",
            "Select Marks",
            choices = UI_MARKS(),
            selected = UI_MARKS()
        )
    })
    
    output$genePlot <- renderPlot({
        plot_velocity_arrows_selected(DATA()$tsne_dt,
                                      DATA()$query_gr,
                                      input$selCells,
                                      tss_ids = input$selGenes) +
            coord_fixed() + theme_classic() + labs(x = "", y = "")
    })
    
    sel_zoom_xrng = reactiveVal(c(-.5, .5))
    sel_zoom_yrng = reactiveVal(c(-.5, .5))
    
    plot_zoom_xrng = reactiveVal(c(-.5, .5))
    plot_zoom_yrng = reactiveVal(c(-.5, .5))
    
    observeEvent(input$global_brush, {
        if(is.null(input$global_brush)){
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
    
    zoom_id = reactiveVal(NULL)
    
    observe({
        req(DATA())
        xrng = sel_zoom_xrng()
        yrng = sel_zoom_yrng()
        samp_id = sampleCap(DATA()$tsne_dt[tx >= min(xrng) & tx <= max(xrng) & ty >= min(yrng) & ty <= max(yrng), ]$id, 100)
        zoom_id(samp_id)
    })
    
    output$detailPlot = renderPlot({
        req(input$n_detail)
        req(input$selMarksDetail)
        req(input$selCellsDetail)
        n_detail = input$n_detail
        view_size = input$detail_view_size
        qdt = config_dt[mark %in% input$selMarksDetail & 
                            cell %in% input$selCellsDetail, ]
        qdt$norm_factor = 1
        if(input$sel_detail_type == "sample"){
            samp_id = zoom_id()[seq(n_detail)]
            qgr = subset(DATA()$query_gr, id %in% samp_id)
            if(is.null(qdt$file)){
                prof_dt = DATA()$profile_dt[id %in% samp_id]
            }else{
                # browser()
                prof_dt = stsFetchTsneInput(qdt[, .(file, tall_var, wide_var, cell, mark, norm_factor)], 
                                            cap_value = Inf,
                                            qgr = resize(qgr, view_size, fix = "center"), 
                                            qmet = "sample",
                                            # qwin = 50, 
                                            qwin = round(view_size / 100),
                                            skip_checks = TRUE)$bw_dt    
            }
            
            prof_dt$id = factor(prof_dt$id, levels = samp_id)
            ggplot(prof_dt, aes(x = x, y = y, color = wide_var, group = paste(id, cell, mark))) +
                geom_path() +
                scale_color_manual(values = custom_colors$get_curr_col()) + 
                facet_grid("cell~id") +
                scale_x_continuous(breaks = 0) +
                labs(x = "relative position", y = "fold-enrichment") +
                theme(axis.text.x = element_blank())
        }else if(input$sel_detail_type == "aggregate"){
            samp_id = zoom_id()
            qgr = subset(DATA()$query_gr, id %in% samp_id)
            if(is.null(qdt$file)){
                prof_dt = DATA()$profile_dt
            }else{
                # browser()
                prof_dt = stsFetchTsneInput(qdt[, .(file, tall_var, wide_var, cell, mark, norm_factor)], 
                                            cap_value = Inf,
                                            qgr = resize(qgr, view_size, fix = "center"), 
                                            qmet = "sample",
                                            # qwin = 50, 
                                            qwin = round(view_size / 100),
                                            skip_checks = TRUE)$bw_dt    
            }
            # browser()
            agg_dt = prof_dt[, .(y = mean(y)), .(x, cell, mark, wide_var)]
            # ggplot(agg_dt, aes)
            ggplot(agg_dt, aes(x = x, y = y, color = wide_var, group = paste(cell, mark))) +
                geom_path() +
                scale_color_manual(values = custom_colors$get_curr_col()) + 
                facet_grid("cell~.") +
                scale_x_continuous(breaks = 0) +
                labs(x = "relative position", y = "fold-enrichment") +
                theme(axis.text.x = element_blank())
        }else{
            stop("unrecognized input$sel_detail_type")
        }
    })
    
    # })
    gene_dt = reactiveVal(NULL)
    
    observeEvent({
        sel_zoom_xrng()
        sel_zoom_yrng()
        DATA()
    }, {
        xrng = sel_zoom_xrng()
        yrng = sel_zoom_yrng()
        reg_id = DATA()$tsne_dt[tx >= min(xrng) & tx <= max(xrng) & ty >= min(yrng) & ty <= max(yrng), ]$id
        dt = unique(DATA()$annotation_dt[id %in% reg_id][distance < 1e6][order(distance)])
        gene_dt(dt)
    })
    
    output$tableGenes = DT::renderDT(expr = {
        DT::datatable(gene_dt(), 
                      extensions = c('Scroller'),
                      options = list(pageLength = 15,
                                     scrollY = 400,
                                     scroller = TRUE), 
                      filter = "top")
    }, server = TRUE)
    
    
    output$geneQPlot = renderPlot({
        max_dist = 3e4
        gene_text = input$textGeneQ
        gene_list = strsplit(toupper(gene_text), "[ ,]+")[[1]]
        hit_id = unique(DATA()$annotation_dt[gene_name %in% gene_list & distance <= max_dist]$id)
        ggplot() + 
            annotate("point", x =  DATA()$tsne_dt$tx, y = DATA()$tsne_dt$ty, size = .2, color = "gray") +
            annotate("point", x =  DATA()$tsne_dt[id %in% hit_id]$tx, y = DATA()$tsne_dt[id %in% hit_id]$ty, size = .2, color = "red")
    })
    
    observeEvent({
        input$dlGeneTable
    }, {
        showModal(
            modalDialog(title = "Dowload Gene Table", 
                        footer = fluidRow(
                            downloadButton("dlCsv", label = "Download csv"),
                            downloadButton("dlXlsx", label = "Download xlsx")
                        ),
                        textInput("txtPrefix", 
                                  label = "File Prefix", 
                                  value = "tsne_selected_genes")
            )
        )
    })
    
    output$dlCsv = downloadHandler(
        filename = function(){
            pre = input$txtPrefix
            if(pre == ""){
                pre = "tsne_selected_genes"
            }
            pre = paste0(pre, ".csv")
            pre
        }, 
        content = function(file){
            chk = input$dlUnique
            if(chk){
                odt = unique(gene_dt()[input$tableGenes_rows_all, .(gene_name)])
                
            }else{
                odt = gene_dt()[input$tableGenes_rows_all,]
            }
            fwrite(odt, file)
            browser()
            
            removeModal()
        }
    )
    
    output$dlXlsx = downloadHandler(
        filename = function(){
            pre = input$txtPrefix
            if(pre == ""){
                pre = "tsne_selected_genes"
            }
            pre = paste0(pre, ".xlsx")
            pre
        }, 
        content = function(file){
            chk = input$dlUnique
            if(chk){
                odt = unique(gene_dt()[input$tableGenes_rows_all, .(gene_name)])
                
            }else{
                odt = gene_dt()[input$tableGenes_rows_all,]
            }
            openxlsx::write.xlsx(odt, file)
            removeModal()
        }
    )
    
    output$textCountUnique = renderText({
        paste(length(gene_dt()[input$tableGenes_rows_all,]$gene_name), "unique genes")
    })
    
    
})

