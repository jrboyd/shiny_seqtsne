server_loading = function(input, output, session, 
                          DATA,
                          UI_TALLV,
                          UI_WIDEV,
                          UI_CELLS,
                          UI_MARKS,
                          UI_GENES){
    
    observeEvent({#watch data and update reactives used by UI elements
        DATA()
    }, {
        UI_TALLV(unique(DATA()$config_dt$tall_var))
        UI_WIDEV(unique(DATA()$config_dt$wide_var))
        UI_GENES(unique(DATA()$query_gr$gene_name))
        UI_CELLS(unique(DATA()$config_dt$cell))
        UI_MARKS(unique(DATA()$config_dt$mark))
    })
    
    showModal(modalDialog(    
        selectInput("selDataset", "Select Dataset", choices = UI_DATASETS),
        footer = tagList(
            actionButton("btnLoadDataset", label = "Load")    
        )
    ))
    
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
        message(paste("end", sel))
    })
    
    observeEvent({#show/hide loadPrompt based on if DATA is set
        DATA()
    }, {
        if(is.null(DATA())){
            shinyjs::show("loadPrompt")
        }else{
            shinyjs::hide("loadPrompt")
        }
    })
}