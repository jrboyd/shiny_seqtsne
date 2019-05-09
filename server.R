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
    
    shinyjs::hide("libs-content")#, anim = TRUE)#, animType = "slide")
    shinyjs::show("waiting-content")#, animType = "slide")
    
    UI_TALLV = reactiveVal(NULL)
    UI_WIDEV = reactiveVal(NULL)
    UI_CELLS = reactiveVal(NULL)
    UI_MARKS = reactiveVal(NULL)
    UI_GENES = reactiveVal(NULL)
    DATA = reactiveVal(NULL)
    GLOBAL_PLOT = reactiveVal(NULL)
    DETAIL_PLOT = reactiveVal(NULL)
    
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
                     get_curr_col = custom_colors$get_curr_col,
                     GLOBAL_PLOT)
    server_debug(input, output, session, 
                 plot_zoom_xrng, plot_zoom_yrng,
                 sel_zoom_xrng, sel_zoom_yrng)
    observeEvent({
        input$btnDebug
    }, {
        browser()
    })
    
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
    samp_id = reactiveVal(NULL)
    detail_id = reactiveVal(NULL)
    
    observe({
        req(DATA())
        xrng = sel_zoom_xrng()
        yrng = sel_zoom_yrng()
        # samp_id = sampleCap(DATA()$tsne_dt[tx >= min(xrng) & tx <= max(xrng) & ty >= min(yrng) & ty <= max(yrng), ]$id, 100)
        # zoom_id(samp_id)
        z_id = DATA()$tsne_dt[tx >= min(xrng) & tx <= max(xrng) & ty >= min(yrng) & ty <= max(yrng), ]$id
        zoom_id(z_id)
        samp_id(sampleCap(z_id, 100))
    })
    
    observe({
        req(samp_id())
        req(input$n_detail)
        detail_id(samp_id()[seq(input$n_detail)])
    })
    
    output$detailPlot = renderPlot({
        req(input$n_detail)
        req(input$selMarksDetail)
        req(input$selCellsDetail)
        req(input$detail_file_type)
        n_detail = input$n_detail
        view_size = input$detail_view_size
        qdt = config_dt[mark %in% input$selMarksDetail & 
                            cell %in% input$selCellsDetail, ]
        qdt$norm_factor = 1
        if(input$sel_detail_type == "sample"){
            samp_id = detail_id()
        }else{
            samp_id = zoom_id()
        }
        qgr = subset(DATA()$query_gr, id %in% samp_id)
        qgr = resize(qgr, view_size, fix = "center")
        
        if(input$detail_file_type == "bam"){
            if(input$sel_detail_type == "sample"){
                
                qdt[, file := bam_file]
                if(is.null(qdt$file)){
                    prof_dt = DATA()$profile_dt[id %in% samp_id]
                    ggplot(prof_dt, aes(
                        x = x,
                        y = y,
                        color = wide_var,
                        group = paste(id, cell, mark, strand)
                    )) +
                        geom_path() +
                        scale_color_manual(values = custom_colors$get_curr_col()) + 
                        facet_grid(cell + mark ~ id, scales = "free_y") +
                        scale_x_continuous(breaks = 0) +
                        labs(x = "relative position", y = "fold-enrichment") +
                        theme(axis.text.x = element_blank())
                }else{
                    p = switch(input$detail_bam_plot_ype,
                               stranded = {
                                   prof_dt = stsFetchTsneInput(
                                       qdt[, .(file, tall_var, wide_var, cell, mark, norm_factor)], 
                                       qgr = qgr,
                                       cap_value = Inf,
                                       qmet = "summary",
                                       fetch_FUN = seqtsne::fetch_bam_stranded_dt,
                                       # qwin = 50,
                                       qwin = 100, 
                                       high_on_right = FALSE,
                                       # qwin = round(view_size / 100),
                                       skip_checks = TRUE
                                   )$bw_dt    
                                   
                                   
                                   prof_dt$id = factor(prof_dt$id, levels = samp_id)
                                   ggplot(prof_dt, aes(
                                       x = x,
                                       y = y,
                                       color = strand,
                                       group = paste(id, cell, mark, strand)
                                   )) +
                                       geom_path() +
                                       scale_color_manual(values = c("+" = "red", "-" = "black")) +
                                       facet_grid(cell + mark ~ id, scales = "free_y") +
                                       scale_x_continuous(breaks = 0) +
                                       labs(x = "relative position", y = "fold-enrichment") +
                                       theme(axis.text.x = element_blank())
                               },
                               pileup = {
                                   prof_dt = stsFetchTsneInput(
                                       qdt[, .(file, tall_var, wide_var, cell, mark, norm_factor)],
                                       qgr = qgr,
                                       cap_value = Inf,
                                       qmet = "summary",
                                       # fetch_FUN = seqtsne::fetch_bam_stranded_dt,
                                       # qwin = 50,
                                       qwin = 100, 
                                       high_on_right = FALSE,
                                       # qwin = round(view_size / 100),
                                       skip_checks = TRUE
                                   )$bw_dt    
                                   
                                   
                                   prof_dt$id = factor(prof_dt$id, levels = samp_id)
                                   ggplot(prof_dt, aes(
                                       x = x,
                                       y = y,
                                       color = wide_var,
                                       group = paste(id, cell, mark, strand)
                                   )) +
                                       geom_path() +
                                       scale_color_manual(values = custom_colors$get_curr_col()) + 
                                       facet_grid(cell + mark ~ id, scales = "free_y") +
                                       scale_x_continuous(breaks = 0) +
                                       labs(x = "relative position", y = "fold-enrichment") +
                                       theme(axis.text.x = element_blank())
                               }, 
                               scc = {
                                   names(qgr) = qgr$id
                                   shift_dt = mclapply(qdt$file, function(f){
                                       rl = seqsetvis:::getReadLength(f, query_gr = qgr)
                                       shift_dt = seqsetvis::crossCorrByRle(
                                           f,
                                           query_gr = qgr,
                                           fragment_sizes = seq(100, 500, by = 5),
                                           read_length = rl
                                       )
                                       shift_dt$file = f
                                       shift_dt = merge(shift_dt, qdt, by = "file")
                                       shift_dt$file = NULL
                                       shift_dt
                                   }) %>% rbindlist
                                   
                                   shift_dt$id = factor(shift_dt$id, levels = samp_id)
                                   shift_max = shift_dt[, .(fragLen = shift[which.max(correlation)], correlation = max(correlation)), by = .(id, cell, mark, wide_var)]
                                   
                                   ggplot() + 
                                       geom_path(data = shift_dt, 
                                                 aes(
                                                     x = shift,
                                                     y = correlation,
                                                     group = cell,
                                                     color = wide_var
                                                 )) + 
                                       geom_point(data = shift_max, 
                                                  aes(x = fragLen, y = correlation),
                                                  color = "black") + 
                                       scale_color_manual(values = custom_colors$get_curr_col()) + 
                                       coord_cartesian(ylim = c(0, 1)) +
                                       facet_grid(mark+cell~id)
                               }, 
                               shiftmap = {
                                   
                                   names(qgr) = qgr$id
                                   shift_dt = mclapply(qdt$file, function(f){
                                       rl = tryCatch({
                                           seqsetvis:::getReadLength(f, query_gr = qgr)
                                       }, error = function(e)NULL)
                                       if(is.null(rl)){
                                           shift_dt = data.table(shift = 1, id = factor(qgr$id),  correlation = 0, file = f)
                                           shift_dt = merge(shift_dt, qdt, by = "file")
                                           shift_dt$file = NULL
                                           shift_dt[numeric()]
                                       }else{
                                           shift_dt = seqsetvis::crossCorrByRle(
                                               f,
                                               query_gr = qgr,
                                               fragment_sizes = seq(100, 500, by = 5),
                                               read_length = rl
                                           )
                                           shift_dt$file = f
                                           shift_dt = merge(shift_dt, qdt, by = "file")
                                           shift_dt$file = NULL
                                           shift_dt
                                       }
                                   }) %>% rbindlist
                                   
                                   
                                   shift_max = shift_dt[, .(fragLen = shift[which.max(correlation)], correlation = max(correlation)), by = .(id, cell, mark, wide_var)]
                                   
                                   step = 10
                                   prof_dt = stsFetchTsneInput(qdt, 
                                                               qwin = step/2, 
                                                               qmet = "sample",
                                                               qgr = qgr, 
                                                               cap_value = Inf, 
                                                               high_on_right = FALSE,
                                                               fetch_FUN = fetch_bam_stranded_dt)$bw_dt
                                   qshift = seq(0, 500/2, step/2)
                                   all_shift = mclapply(qshift, function(s){
                                       off_dt = dcast(prof_dt[, 
                                                              .(cell, mark, wide_var, id, 
                                                                x = ifelse(strand == "+", x + s, x - s), y, strand)], 
                                                      cell+mark+id+x+wide_var~strand, value.var = "y")
                                       off_dt[, delta := as.numeric(`+` - `-`)]
                                       # off_dt[!is.na(delta)]
                                       off_dt$shift = s
                                       off_dt
                                   }) %>% rbindlist()
                                   
                                   all_shift$id = factor(all_shift$id, levels = samp_id)
                                   
                                   # all_shift[, delta := delta / max(abs(delta), na.rm = TRUE), .(cell, mark)]
                                   ggplot() +
                                       geom_raster(data = all_shift, 
                                                   aes(x = x, y = shift * 2, fill = delta)) +
                                       geom_segment(data = shift_max,
                                                    aes(
                                                        x = min(all_shift$x),
                                                        xend = max(all_shift$x),
                                                        y = fragLen,
                                                        yend = fragLen,
                                                        color = wide_var
                                                    ),
                                                    color = "red") +
                                       scale_color_manual(values = custom_colors$get_curr_col()) + 
                                       facet_grid(cell + mark ~ id) +
                                       scale_fill_viridis_c() + 
                                       labs(x = "bp", y = "shift")
                               }, {
                                   stop("unrecognized plot type: ", input$detail_bam_plot_ype)
                               })
                }
            }else if(input$sel_detail_type == "aggregate"){
                samp_id = zoom_id()
                if(is.null(qdt$file)){
                    prof_dt = DATA()$profile_dt
                }else{
                    prof_dt = stsFetchTsneInput(qdt[, .(file, tall_var, wide_var, cell, mark, norm_factor)], 
                                                cap_value = Inf,
                                                qgr = qgr, 
                                                qmet = "summary",
                                                qwin = 100,
                                                # qwin = round(view_size / 100),
                                                skip_checks = TRUE)$bw_dt    
                }
                agg_dt = prof_dt[, .(y = mean(y)), .(x, cell, mark, wide_var)]
                p = ggplot(agg_dt, aes(x = x, y = y, color = wide_var, group = paste(cell, mark))) +
                    geom_path() +
                    scale_color_manual(values = custom_colors$get_curr_col()) + 
                    facet_grid("cell~.") +
                    scale_x_continuous(breaks = 0) +
                    labs(x = "relative position", y = "fold-enrichment") +
                    theme(axis.text.x = element_blank())
            }else{
                stop("unrecognized input$sel_detail_type")
            }
        }else if(input$detail_file_type == "bigwig"){
            if(input$sel_detail_type == "sample"){
                if(is.null(qdt$file)){
                    prof_dt = DATA()$profile_dt[id %in% samp_id]
                }else{
                    prof_dt = stsFetchTsneInput(qdt[, .(file, tall_var, wide_var, cell, mark, norm_factor)], 
                                                cap_value = Inf,
                                                qgr = qgr, 
                                                qmet = "summary",
                                                # qwin = 50, 
                                                qwin = 100,
                                                # qwin = round(view_size / 100),
                                                skip_checks = TRUE)$bw_dt    
                }
                
                prof_dt$id = factor(prof_dt$id, levels = samp_id)
                p = ggplot(prof_dt, aes(x = x, y = y, color = wide_var, group = paste(id, cell, mark, strand))) +
                    geom_path() +
                    scale_color_manual(values = custom_colors$get_curr_col()) + 
                    facet_grid("cell~id") +
                    scale_x_continuous(breaks = 0) +
                    labs(x = "relative position", y = "fold-enrichment") +
                    theme(axis.text.x = element_blank())
            }else if(input$sel_detail_type == "aggregate"){
                samp_id = zoom_id()
                if(is.null(qdt$file)){
                    prof_dt = DATA()$profile_dt
                }else{
                    prof_dt = stsFetchTsneInput(qdt[, .(file, tall_var, wide_var, cell, mark, norm_factor)], 
                                                cap_value = Inf,
                                                qgr = qgr, 
                                                qmet = "summary",
                                                qwin = 100,
                                                # qwin = round(view_size / 100),
                                                skip_checks = TRUE)$bw_dt    
                }
                agg_dt = prof_dt[, .(y = mean(y)), .(x, cell, mark, wide_var)]
                p = ggplot(agg_dt, aes(x = x, y = y, color = wide_var, group = paste(cell, mark))) +
                    geom_path() +
                    scale_color_manual(values = custom_colors$get_curr_col()) + 
                    facet_grid("cell~.") +
                    scale_x_continuous(breaks = 0) +
                    labs(x = "relative position", y = "fold-enrichment") +
                    theme(axis.text.x = element_blank())
            }else{
                stop("unrecognized input$sel_detail_type")
            }
        }else{
            stop("unrecognized input$detail_file_type : ", input$detail_file_type)
        }
        
        DETAIL_PLOT(p)
        p
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
        dt = unique(DATA()$annotation_dt[id %in% reg_id][distance < 1e5][order(distance)])
        gene_dt(dt)
    })
    
    output$tableGenes = DT::renderDataTable(expr = {
        req(gene_dt())
        DT::datatable(gene_dt(), 
                      extensions = c('Scroller'),
                      options = list(pageLength = 15,
                                     scrollY = 400,
                                     scroller = TRUE), 
                      filter = "top")
    }, server = TRUE)
    
    go_dt = reactiveVal(NULL)
    observeEvent({
        genes_unique()
        input$tabset_gene_tables
    }, {
        req(genes_unique())
        if(input$tabset_gene_tables != "GO"){
            return(NULL)
        }
        df = data.frame(genes_unique(), 1)
        
        if(nrow(df) > 5000){
            df = data.frame("more than 5000 genes selected.")
        }else{
            df = data.frame("disabled")
            # # browser()
            # # library(topGO)
            # # sampleGOdata <- new("topGOdata",
            # #                     description = "Simple session", ontology = "BP",
            # #                     allGenes = unique(annotation_dt$gene_name), 
            # #                     geneSel = as.character(df[[1]]),
            # #                     nodeSize = 10,
            # #                     annot = annFUN.db, affyLib = affyLib)
            # # 
            # # 
            # suppressWarnings({
            #     go_res = GOfuncR::go_enrich(df, silent = TRUE, n_randsets = 1)    
            # })
            # go_dt = as.data.table(go_res$results)
            # df = go_dt[FWER_overrep < .05]
            # # browser()
        }
        
        go_dt(df)
    })
    
    output$tableGO = DT::renderDataTable(expr = {
        req(go_dt())
        dt = go_dt()
        DT::datatable(dt, 
                      extensions = c('Scroller'),
                      options = list(pageLength = 15,
                                     scrollY = 400,
                                     scrollX = ncol(dt),
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
                            actionButton("btnCancelModal", "Cancel"),
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
    
    genes_unique = reactiveVal(NULL)
    
    observe({
        uniq = unique(gene_dt()[input$tableGenes_rows_all,]$gene_name)
        genes_unique(uniq)
        
        output$textCountUnique = renderText({
            paste(length(gene_dt()[input$tableGenes_rows_all,]$gene_name), "unique genes")
        })
    })
    
    observeEvent({
        input$btnDlGlobal
    },{
        def_col = ifelse(is.null(input$colorSelection), 
                         "99ccff55", 
                         input$colorSelection)
        showModal(modalDialog(
            fluidRow(
                radioButtons("selImgFormat", "Image Format", choices = c("png", "pdf")),
                textInput("txtImgPrefix", "Prefix", value = "global_image"),
                colourInput("colorSelection", "Selection Color", value = def_col, allowTransparent = TRUE),
                numericInput("numImgWidth", "Width (in)", value = 6, min = 1, max = 20),
                numericInput("numImgHeight", "Height (in)", value = 6, min = 1, max = 20),
                numericInput("numImgDpi", "DPI", value = 150, min = 50, max = 600)
            ),
            title = "Save Current Global Image",
            footer = fluidRow(
                actionButton("btnCancelModal", label = "Cancel" ),
                downloadButton("btnSaveGlobalImg")
            )
        ))
    })
    
    observeEvent({
        input$selImgFormat
    }, {
        if(input$selImgFormat == "pdf"){
            shinyjs::disable("numImgDpi")
        }else{
            shinyjs::enable("numImgDpi")
        }
    })
    
    observe({
        req(input$btnCancelModal)
        input$btnCancelModal
        removeModal()
    })
    
    output$btnSaveGlobalImg = downloadHandler(
        filename = function(){
            paste0(input$txtImgPrefix, ".", input$selImgFormat)
        }, 
        content = function(file){
            p = GLOBAL_PLOT()
            plot_rect = c(plot_zoom_xrng(), plot_zoom_yrng())
            sel_rect = c(sel_zoom_xrng(), sel_zoom_yrng())
            if(!all(plot_rect == sel_rect)){
                box_col = input$colorSelection
                if(nchar(box_col) == 7){
                    box_col = paste0(box_col, "FF")
                }
                line_col = paste0(substr(box_col, 0, 7), "FF")
                p = p + annotate("rect", 
                                 xmin = sel_rect[1],
                                 xmax = sel_rect[2],
                                 ymin = sel_rect[3],
                                 ymax = sel_rect[4],
                                 fill = box_col,
                                 color = line_col)
            }
            ggsave(
                file,
                plot = p,
                units = "in",
                width = input$numImgWidth,
                height = input$numImgHeight,
                dpi = input$numImgDpi
            )
            removeModal()
        }
    )
    
    observeEvent({
        input$btnDlDetail
    }, {
        showModal(modalDialog(
            fluidRow(
                radioButtons("selImgFormat", "Image Format", choices = c("png", "pdf")),
                textInput("txtImgPrefix", "Prefix", value = "detail_image"),
                numericInput("numImgWidth", "Width (in)", value = 6, min = 1, max = 20),
                numericInput("numImgHeight", "Height (in)", value = 6, min = 1, max = 20),
                numericInput("numImgDpi", "DPI", value = 150, min = 50, max = 600)
            ),
            title = "Save Current Detail Image",
            footer = fluidRow(
                actionButton("btnCancelModal", label = "Cancel" ),
                downloadButton("btnSaveDetailImg")
            )
        ))
    })
    
    output$btnSaveDetailImg = downloadHandler(
        filename = function(){
            paste0(input$txtImgPrefix, ".", input$selImgFormat)
        }, 
        content = function(file){
            p = DETAIL_PLOT()
            ggsave(
                file,
                plot = p,
                units = "in",
                width = input$numImgWidth,
                height = input$numImgHeight,
                dpi = input$numImgDpi
            )
            removeModal()
        }
    )
    
    observeEvent({
        DATA()
    }, {
        DATA()$config_dt
        if(!is.null(DATA()$config_dt$bam_file)){
            shinyjs::show("detail_file_type")
        }
        
    })
    
    
    observeEvent({
        input$detail_file_type
    }, {
        if(input$detail_file_type == "bam"){
            shinyjs::show("detail_bam_plot_ype")
        }else{
            shinyjs::hide("detail_bam_plot_ype")
        }
        
    })
    
    output$uiSelector = renderUI({
        selectInput("selSelector", "Select Data", choices = names(DATA()$selector))
    })
    
    output$tableSelect = DT::renderDataTable({
        req(input$selSelector)
        dt = DATA()$selector[[input$selSelector]]
        dt = dt[id %in% zoom_id()]
        
        rnd_cn = colnames(dt)[sapply(dt, is.numeric)]
        
        dt = DT::datatable(dt, 
                      extensions = c('Scroller'),
                      options = list(pageLength = 15,
                                     scrollY = 400,
                                     scrollX = ncol(dt),
                                     scroller = TRUE), 
                      filter = "top")
        dt = DT::formatRound(dt, rnd_cn, 3)
        dt
    }, server = TRUE)

# observeEvent({
#     input$detail_file_type
# }, {
#     showNotification(input$detail_file_type)
# })
})

