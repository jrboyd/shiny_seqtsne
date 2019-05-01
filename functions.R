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
