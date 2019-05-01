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
    
    showModal(modalDialog(    
        selectInput("selDataset", "Select Dataset", choices = UI_DATASETS),
        footer = tagList(
            actionButton("btnLoadDataset", label = "Load")    
        )
    ))
    
    
    # library(GenomicRanges)
    # library(data.table)
    # library(seqsetvis)
    # library(seqtsne)
    
    UI_TALLV = reactiveVal(NULL)
    UI_WIDEV = reactiveVal(NULL)
    UI_CELLS = reactiveVal(NULL)
    UI_MARKS = reactiveVal(NULL)
    UI_GENES = reactiveVal(NULL)
    DATA = reactiveVal(NULL)
    
    observeEvent({
        input$btnLoadDataset
    }, {
        # runjs("var today = new Date(); alert(today);")
        shinyjs::runjs({
            'function sleep(ms) {
               return new Promise(resolve => setTimeout(resolve, ms));
             }
             async function myFunction() {
             var text = "";
             var i = 0;
             document.getElementById("loading-status").innerHTML = "Loading";
             tps = 5
             while (1 < 10) {
               await sleep(1000/tps);
               txt = document.getElementById("loading-status").innerHTML
               nxt = "."
               if(i % tps == tps - 1) nxt = (i + 1) / tps;
               document.getElementById("loading-status").innerHTML = txt + nxt;
               if(getComputedStyle(document.getElementById("loading-content")).display == "none")
                 break;
               i++;
             }
             }
             myFunction()
             '
        })
        removeModal()
        shinyjs::hide("waiting-content")
        shinyjs::show("loading-content")
        sel = input$selDataset
        message(paste("start", sel))
        src = UI_DATASOURCES[sel]
        res = load_dataset(src)
        DATA(res)
        shinyjs::hide("loading-content")
        shinyjs::show("app-content")
        # hidden = sapply(req_vars, function(x)remove(x))
        message(paste("end", sel))
        # sapply(req_vars, function(x)exists(x))
    })
    
    observeEvent({
        DATA()
    }, {
        if(is.null(DATA())){
            shinyjs::show("loadPrompt")
        }else{
            shinyjs::hide("loadPrompt")
        }
    })
    
    observeEvent({
        DATA()
    }, {
        UI_TALLV(unique(DATA()$config_dt$tall_var))
        UI_WIDEV(unique(DATA()$config_dt$wide_var))
        UI_GENES(unique(DATA()$query_gr$gene_name))
        UI_CELLS(unique(DATA()$config_dt$cell))
        UI_MARKS(unique(DATA()$config_dt$mark))
    })
    
    output$ui_global_cells = renderUI({
        req(UI_TALLV())
        checkboxGroupInput(
            "selCells",
            "Select Cells",
            choices = UI_TALLV(),
            selected = UI_TALLV()
        )
    })
    
    output$ui_global_marks = renderUI({
        req(UI_WIDEV())
        checkboxGroupInput(
            "selMarks",
            "Select Marks",
            choices = UI_WIDEV(),
            selected = UI_WIDEV()[1]
        )
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
    
    observeEvent({#manage conditional UI elements based on globalViewType
        input$globalViewType
    }, {
        if(input$globalViewType %in% c(GLOBAL_VIEW_POINTS, 
                                       GLOBAL_VIEW_PROFILES_FAST, 
                                       GLOBAL_VIEW_PROFILES_SLOW)){
            shinyjs::show("app-tall-wide-sel")
        }else{
            shinyjs::hide("app-tall-wide-sel")
        }
    })
    
    output$globalPlot <- renderPlot({
        req(input$selCells)
        req(input$selMarks)
        typ = input$globalViewType
        
        xrng = plot_zoom_xrng()
        yrng = plot_zoom_yrng()
        frac_shown = max(c(diff(xrng), diff(yrng)))
        point_size = .1 / frac_shown
        in_view_id =  DATA()$tsne_dt[tx >= min(xrng) & tx <= max(xrng) & ty >= min(yrng) & ty <= max(yrng)]$id
        if(!input$selNumPlotted == "all"){
            set.seed(0)
            in_view_id = sampleCap(in_view_id, as.numeric(input$selNumPlotted))
        }
        # tp = 
        if(typ == GLOBAL_VIEW_POINTS){
            if(length(input$selMarks) >= 1){
                if(length(input$selMarks) == 2){
                    corr_dt = dcast(
                        DATA()$agg_dt[id %in% in_view_id][tall_var %in% input$selCells & wide_var %in% input$selMarks], 
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
                        # coord_cartesian(xlim = xrng, ylim = yrng) +
                        facet_wrap("variable") +
                        labs(color = "difference",
                             fill = paste("min of", colnames(corr_dt)[5],
                                          "\nand", colnames(corr_dt)[6]))
                    
                    # p = ggplot(corr_dt, aes(x = tx, y = ty, color = diff)) +
                    #     geom_point(size = point_size) + 
                    #     scale_color_gradientn(colors = c("blue", "white", "red"), 
                    #                           breaks = c(-1, -.5, 0, .5, 1),
                    #                           labels = c(colnames(corr_dt)[5], "", "-", "", colnames(corr_dt)[6]),
                    #                           limits = c(-1, 1)) +
                    #     coord_cartesian(xlim = xrng, ylim = yrng) +
                    #     labs(color = "difference")
                    # 
                    # ggplot(corr_dt, aes(x = tx, y = ty, color = agree)) +
                    #     geom_point(size = point_size) + 
                    #     scale_color_gradientn(colors = c("blue", "white", "red"), 
                    #                           breaks = c(-1, -.5, 0, .5, 1),
                    #                           labels = c(colnames(corr_dt)[5], "", "-", "", colnames(corr_dt)[6]),
                    #                           limits = c(-1, 1)) +
                    #     coord_cartesian(xlim = xrng, ylim = yrng) +
                    #     labs(color = "difference")
                }else{
                    nc = ceiling(sqrt(length(input$selMarks)))
                    p = ggplot(
                        DATA()$agg_dt[id %in% in_view_id][tall_var %in% input$selCells & wide_var %in% input$selMarks], 
                        aes(x = tx, y = ty, color = value)) +
                        geom_point(size = point_size) +  
                        coord_cartesian(xlim = xrng, ylim = yrng) +
                        facet_wrap("wide_var", ncol = nc)    
                }
            }else{
                p = ggplot(DATA()$tsne_dt[id %in% in_view_id][tall_var %in% input$selCells], aes(x = tx, y = ty)) +
                    coord_cartesian(xlim = xrng, ylim = yrng) +
                    geom_point(size = point_size)
            }
            
        }else if(typ == GLOBAL_VIEW_DENSITY){
            p = ggplot(DATA()$tsne_dt[tall_var %in% input$selCells], aes(x = tx, y = ty)) +
                coord_cartesian(xlim = xrng, ylim = yrng) +
                geom_density2d() + 
                facet_wrap("tall_var") 
        }else if(typ == GLOBAL_VIEW_PROFILES_FAST){
            # if(input$selGlobalColoring == "mark"){
            #     tmp = unique(DATA()$profile_dt[, .(wide_var, mark)])
            #     cm = DATA()$color_mapping_byMark[tmp$mark]
            #     names(cm) = tmp$wide_var    
            # }else if(input$selGlobalColoring == "both"){
            #     # brows er()
            #     # tmp = unique(DATA()$profile_dt[, .(wide_var, cell, mark)])
            #     cm = DATA()$color_mapping_byMark_byWide[unique(DATA()$profile_dt$wide_var)]
            #     names(cm) = unique(DATA()$profile_dt$wide_var)
            # }else{
            #     tmp = unique(DATA()$profile_dt[, .(wide_var, cell)])
            #     cm = DATA()$color_mapping_byMark_byCell[tmp$cell]
            #     names(cm) = tmp$wide_var    
            # }
            cm = get_curr_col()
            p = stsPlotSummaryProfiles(DATA()$profile_dt[id %in% in_view_id, ],
                                       DATA()$tsne_dt[id %in% in_view_id, ], 
                                       q_tall_vars = input$selCells,
                                       q_wide_vars = input$selMarks,
                                       x_points = input$numBins,
                                       xrng = xrng,
                                       yrng = yrng,
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
                                       ylim = DATA()$tylim,
                                       line_color_mapping = cm,
                                       plot_type = "raster")
        }
        p +
            theme_classic() + 
            labs(x = "", y = "")
        
    })
    
    output$zoomPlot <- renderPlot({
        ggplot()
        # zimg_res = make_tsne_img(DATA()$profile_dt[cell %in%  input$selCells],
        #                          DATA()$tsne_dt, n_points = input$bins,
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
        plot_velocity_arrows_selected(DATA()$tsne_dt,
                                      DATA()$query_gr,
                                      input$selCells,
                                      tss_ids = input$selGenes) +
            coord_fixed() + theme_classic() + labs(x = "", y = "")
    })
    
    output$profilePlot <- renderPlot({
        plot_profiles_selected(DATA()$profile_dt,
                               DATA()$query_gr,
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
            vel_plots = plot_velocity_arrows(DATA()$tsne_dt, cp[1], cp[2])
        }
        vel_plots[[1]] + theme_classic()
    })
    
    output$pairKey = renderPlot({
        # req(cellPair())
        # cp = cellPair()
        cp = c("H7", "CD34")
        if(length(cp) == 2){
            vel_plots = plot_velocity_arrows(DATA()$tsne_dt, cp[1], cp[2])
        }
        vel_plots[[2]]
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
                scale_color_manual(values = get_curr_col()) + 
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
                scale_color_manual(values = get_curr_col()) + 
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
    
    observeEvent({
        input$btnCustomColors
    }, {
        curr_col = get_curr_col()
        showModal(modalDialog(
            gen_color_picker_ui(
                picker_names = names(curr_col), 
                initial_colors = curr_col
            ),
            title = "Customize Color", 
            footer = fluidRow(
                actionButton("btnCancelCustomColors", label = "Cancel"),
                actionButton("btnApplyCustomColors", label = "OK")
            )
        ))
    })
    
    get_curr_col = function(){
        if(input$selGlobalColoring == "mark"){
            cm = DATA()$color_mapping_byMark
        }else if(input$selGlobalColoring == "both"){
            cm = DATA()$color_mapping_byBoth
        }else if(input$selGlobalColoring == "cell"){
            cm = DATA()$color_mapping_byCell
        }else{
            stop("unable to get_curr_col, unrecognized input$selGlobalColoring")
        }
        cm
    }
    
    set_curr_col = function(cm){
        res = DATA()
        if(input$selGlobalColoring == "mark"){
            stopifnot(names(res$color_mapping_byMark) == names(cm))
            res$color_mapping_byMark = cm
        }else if(input$selGlobalColoring == "both"){
            stopifnot(names(res$color_mapping_byMark_byWide) == names(cm))
            res$color_mapping_byMark_byWide = cm
        }else if(input$selGlobalColoring == "cell"){
            stopifnot(names(res$color_mapping_byMark_byCell) == names(cm))
            res$color_mapping_byMark_byCell = cm
        }else{
            stop("unable to set_curr_col, unrecognized input$selGlobalColoring")
        }
        DATA(res)
    }
    
    observeEvent({
        input$btnCancelCustomColors
    }, {
        removeModal()
    })
    
    observeEvent({
        input$btnApplyCustomColors
    }, {
        cm = get_curr_col()
        cm_new = fetch_color_picker_ui(names(cm), input)
        set_curr_col(cm_new)
        removeModal()
    })
})

