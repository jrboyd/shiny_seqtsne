gen_color_picker_ui = function(grps, brew_name = "Dark2", is_free_color = FALSE){
    browser()
    possibleColors = unique(safeBrew(50, brew_name))
    colorChoices = safeBrew(length(grps), brew_name)
    #special case for gradients where all color choices aren't in possible colors
    if(!all(colorChoices %in% possibleColors)){
        clen = length(possibleColors)
        glen = length(grps) + 1
        gi = round((1:glen-1) * ((clen-1) / (glen - 1)) + 1)
        colorChoices = possibleColors[gi]
    }
    
    # colorChoices = apply(col2rgb(colorChoices),2, function(x)paste(x, collapse = ","))
    names(colorChoices) = grps
    # rv$key2color = colorChoices
    inputs = character(length(grps))
    
    for(i in seq_along(grps)){
        g = as.character(grps[i])
        g = gsub("/", "_", g)
        
        if(is_free_color){
            ginput = colourInput(inputId = paste0("color_", g),
                                 label = g,
                                 value = colorChoices[i])
        }else{
            ginput = colourInput(inputId = paste0("color_", g),
                                 label = g,
                                 value = colorChoices[i],
                                 palette = "limited",
                                 allowedCols = possibleColors)
        }
        
        inputs[i] = as.character(ginput)
    }
    HTML(paste(inputs, collapse = "</br>"))
}

fetch_color_picker_ui = function(grps){
    ids = paste0("color_", grps)
    cols = sapply(ids, function(id){
        input[[id]]
    })
    names(cols) = grps
    cols
}