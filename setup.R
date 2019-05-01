# if(!exists("LOADED")){
# LOADED = TRUE

suppressPackageStartupMessages({
    library(shiny)
    library(shinyjs)
    library(magrittr)
    library(shinycssloaders)
    library(colourpicker)
    # library(GenomicRanges)
    # library(data.table)
    # library(seqsetvis)
    # library(seqtsne)
})

source("functions.R")

load_libs_withProgress = function(libs, session){
    tmp.sink = suppressPackageStartupMessages({
        # libs = c(
        #     "data.table",
        #     "rtracklayer",
        #     "GenomicRanges",
        #     "ggplot2",
        #     "cowplot",
        #     "seqsetvis"
        # )
        withProgress(session = session, 
                     message = 'Loading libraries', 
                     value = 0, 
                     max = length(libs),  
                     expr = {
                         for(i in seq_along(libs)){
                             incProgress(session = session,
                                         amount = 1, 
                                         message = libs[i],
                                         detail = paste0("(", i, "/", length(libs), ")"))
                             library(libs[i], character.only = TRUE)
                         }
                         
                         
                     })
    })
}


# source("generate_data.R")

options(mc.cores = 16)
message(getwd())

D_EXAMPLE = "example"
D_HESC_CD34 = "hESC and CD34 bivalency"
D_BCELL = "B cell bivalency"
D_MSC_TIME = "MSC timecourse bivalency"
D_WALDRON_CONSENSUS = "Waldron hESC+CD34 consensus"
UI_DATASETS = c(D_EXAMPLE, D_WALDRON_CONSENSUS, D_HESC_CD34, D_BCELL, D_MSC_TIME)
data_dir = "~/ShinyApps/shiny_seqtsne"
UI_DATASOURCES = c(file.path(data_dir, "dataset_scripts/generate_example.R"),
                   file.path(data_dir, "dataset_scripts/generate_data_hESC_and_CD34_consensus.R"),
                   file.path(data_dir, "dataset_scripts/generate_data_hESC_and_CD34_biv_wide_full.R"),
                   file.path(data_dir, "dataset_scripts/generate_data_Bcell_biv_wide_full.R"),
                   file.path(data_dir, "dataset_scripts/generate_data_MSC_biv_wide_full.R"))
stopifnot(file.exists(UI_DATASOURCES))
names(UI_DATASOURCES)= UI_DATASETS

UI_DATASETS = UI_DATASETS

GLOBAL_VIEW_POINTS = "points"
GLOBAL_VIEW_PROFILES_FAST = "profiles (fast)"
GLOBAL_VIEW_PROFILES_SLOW = "profiles (slow)"
GLOBAL_VIEW_DENSITY = "density"

# biv19 = fread("PMC5354816_bivalent_hg19.tsv", col.names = c("seqnames", "start", "end")) %>% GRanges
# k4me319 = fread("PMC5354816_k4m3_hg19.tsv", col.names = c("seqnames", "start", "end")) %>% GRanges
# k27me319 = fread("PMC5354816_k27me3_hg19.tsv", col.names = c("seqnames", "start", "end")) %>% GRanges
#
# ch = rtracklayer::import.chain("/slipstream/galaxy/data/hg19/liftOver/hg19ToHg38.over.chain")
# biv38 = rtracklayer::liftOver(biv19, ch) %>% unlist
# k27me338 = rtracklayer::liftOver(k27me319, ch) %>% unlist
# k4me338 = rtracklayer::liftOver(k4me319, ch) %>% unlist
#
# biv_ids = subsetByOverlaps(qgr, biv38, minoverlap = 500)$id
# tsne_dt[, is_biv := id %in% biv_ids]
#
# k27me3_ids = subsetByOverlaps(qgr, k27me338, minoverlap = 500)$id
# tsne_dt[, is_k27me3 := id %in% k27me3_ids]
#
# k4me3_ids = subsetByOverlaps(qgr, k4me338, minoverlap = 500)$id
# tsne_dt[, is_k4me3 := id %in% k4me3_ids]
#
# fmemb = ssvFactorizeMembTable(tsne_dt[, .(is_biv, is_k27me3, is_k4me3)])
# tsne_dt$group = fmemb$group
#
# ggplot(tsne_dt[cell %in% c("H7", "CD34")], aes(x = tx, y = ty, color = cell)) +
#     facet_grid(".~group") + geom_point(alpha = .1, shape = 16, size = .5) +
#     theme_classic()

# }

load_dataset = function(src){
    source(src)
    if(!exists("color_mapping") || is.null(color_mapping)){
        cm = unique(config_dt$mark)
        color_mapping = safeBrew(length(cm), pal = "Set1")
        names(color_mapping) = cm
    }
    if(setequal(names(color_mapping), unique(config_dt$mark)) &&
       !setequal(names(color_mapping), unique(config_dt$wide_var))){
            tmp = unique(config_dt[, .(wide_var, mark)])
            color_mapping = color_mapping[tmp$mark]
            names(color_mapping) = tmp$wide_var
    }
    stopifnot(setequal(names(color_mapping), unique(config_dt$wide_var)))
    
    if(!exists("color_mapping.alt") || is.null(color_mapping.alt)){
        cm_alt = unique(config_dt$cell)
        color_mapping.alt = safeBrew(length(cm_alt), pal = "Dark2")
        names(color_mapping.alt) = cm_alt
    }
    if(setequal(names(color_mapping.alt), unique(config_dt$cell)) &&
       !setequal(names(color_mapping.alt), unique(config_dt$wide_var))){
        tmp = unique(config_dt[, .(wide_var, cell)])
        color_mapping.alt = color_mapping.alt[tmp$cell]
        names(color_mapping.alt) = tmp$wide_var
    }
    stopifnot(setequal(names(color_mapping.alt), unique(config_dt$wide_var)))
    
    if(!exists("color_mapping.full") || is.null(color_mapping.full)){
        cm_full = unique(config_dt$wide_var)
        color_mapping.full = safeBrew(n = length(cm_full), pal = "Set3")
        names(color_mapping.full) = cm_full
    }
    stopifnot(setequal(names(color_mapping.full), unique(config_dt$wide_var)))
    
    
    if(!exists("tylim") || is.null(color_mapping.full)){
        tylim = range(profile_dt$y)
    }
    
    setkey(profile_dt, "id")
    setkey(tsne_dt, "id")
    res = list(
        profile_dt = profile_dt,
        tsne_dt = tsne_dt,
        query_gr = query_gr,
        agg_dt = agg_dt,
        annotation_dt = annotation_dt,
        config_dt = config_dt,
        color_mapping_byMark = color_mapping,
        color_mapping_byCell = color_mapping.alt,
        color_mapping_byBoth = color_mapping.full,
        tylim = tylim
    )
    res
}
