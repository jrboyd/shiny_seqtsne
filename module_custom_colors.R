
server_custom_colors = function(input, output, session, DATA){
    get_curr_col = function(){
        if(input$selGlobalColoring == "mark"){
            cm = DATA()$color_mapping_byMark
        }else if(input$selGlobalColoring == "cell"){
            cm = DATA()$color_mapping_byCell
        }else if(input$selGlobalColoring == "both"){
            cm = DATA()$color_mapping_byBoth
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
        }else if(input$selGlobalColoring == "cell"){
            stopifnot(names(res$color_mapping_byCell) == names(cm))
            res$color_mapping_byCell = cm
        }else if(input$selGlobalColoring == "both"){
            stopifnot(names(res$color_mapping_byBoth) == names(cm))
            res$color_mapping_byBoth = cm
        }else{
            stop("unable to set_curr_col, unrecognized input$selGlobalColoring")
        }
        DATA(res)
    }
    
    
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
    
    return(list(get_curr_col = get_curr_col, set_curr_col = set_curr_col))
}

#' generate a series of color picker inputs with colourpicker::colourInput
#'
#' @param picker_names character vector of unique
#' @param initial_colors hex value colors to initialize color pickers with.  must be same length as picker_names or NULL.
#' @param possible_colors 
#' @param brew_name character, valid RColorBrewer palettte name. only used if possible_colors is NULL and is_free_color is FALSE.
#' @param is_free_color if TRUE, overrides possible_colors and brew_name to allow full color palette selection.
#'
#' @return
#' @export
#'
#' @examples
gen_color_picker_ui = function(picker_names, 
                               initial_colors = NULL, 
                               possible_colors = NULL,
                               brew_name = "Dark2", 
                               is_free_color = FALSE){
    if(any(duplicated(picker_names))){
        stop("all picker_names must be unique")
    }
    if(is.null(possible_colors)){
        possible_colors = unique(safeBrew(50, brew_name))    
    }
    stopifnot(length(picker_names) == length(initial_colors) || is.null(initial_colors))
    if(is.null(initial_colors)){
        initial_colors = possible_colors[seq_along(picker_names)]
        #initial_colors = safeBrew(length(picker_names), brew_name)    
    }
    
    possible_colors = union(possible_colors, initial_colors)
    
    #special case for gradients where all color choices aren't in possible colors
    if(!all(initial_colors %in% possible_colors)){
        clen = length(possible_colors)
        glen = length(picker_names) + 1
        gi = round((1:glen-1) * ((clen-1) / (glen - 1)) + 1)
        initial_colors = possible_colors[gi]
    }
    
    # initial_colors = apply(col2rgb(initial_colors),2, function(x)paste(x, collapse = ","))
    names(initial_colors) = picker_names
    # rv$key2color = initial_colors
    inputs = character(length(picker_names))
    
    for(i in seq_along(picker_names)){
        g = as.character(picker_names[i])
        gkey = gsub("[/ -]", "_", g) #remove illegal chars
        
        if(is_free_color){
            ginput = colourpicker::colourInput(inputId = paste0("color_", gkey),
                                               label = g,
                                               value = initial_colors[i])
        }else{
            ginput = colourpicker::colourInput(inputId = paste0("color_", gkey),
                                               label = g,
                                               value = initial_colors[i],
                                               palette = "limited",
                                               allowedCols = possible_colors)
        }
        
        inputs[i] = as.character(ginput)
    }
    HTML(paste(inputs, collapse = "</br>"))
}

fetch_color_picker_ui = function(picker_names, input){
    gkey = gsub("[/ -]", "_", picker_names) #remove illegal chars
    ids = paste0("color_", gkey)
    
    cols = sapply(ids, function(id){
        input[[id]]
    })
    names(cols) = picker_names
    cols
}
