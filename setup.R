if(!exists("LOADED")){
    LOADED = TRUE
    library(shiny)
    library(data.table)
    library(seqsetvis)
    library(seqtsne)
    library(magrittr)
    library(GenomicRanges)
    source("generate_data.R")
    options(mc.cores = 16)
    message(getwd())
    
    UI_CELLS = unique(profile_dt$cell)
    
    UI_MARKS = unique(profile_dt$mark)
    
    UI_GENES = unique(query_gr$gene_name)
    
    n_tp = 5000
    set.seed(1)
    UI_TP = sample(unique(tsne_dt$id), n_tp / length(UI_CELLS))
    tsne_tp = tsne_dt[id %in% UI_TP]
    
    GLOBAL_VIEW_POINTS = "points"
    GLOBAL_VIEW_PROFILES_FAST = "profiles (fast)"
    GLOBAL_VIEW_PROFILES_SLOW = "profiles (slow)"
    GLOBAL_VIEW_DENSITY = "density"
    
    color_mapping = safeBrew(length(UI_MARKS))
    names(color_mapping) = UI_MARKS
    
    mdt = prep_summary(profile_dt, tsne_dt, 12)
    glyph_df = GGally::glyphs(mdt, x_major = "bx", x_minor = "x", y_major = "by", y_minor = "y")
    ggplot(glyph_df, aes(gx, gy, group = paste(gid, mark), color = mark)) +
        geom_path() +
        scale_color_manual(values =  color_mapping)
    
    names(query_gr) = query_gr$id
    
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
    
}
