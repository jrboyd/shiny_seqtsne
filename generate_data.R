library(magrittr)
library(seqtsne)
library(seqsetvis)
library(GenomicRanges)
library(cowplot)
bfc = BiocFileCache::BiocFileCache()

res = bfcif(bfc, "KZ_seqtsne_bookmarking", function(){
    #bigwigs matching Josh's published peak call
    wd = "/slipstream/galaxy/uploads/working/qc_framework/output_JR_bookmarking_published"
    bw_book = dir(wd, full.names = TRUE, pattern = "ed$") %>% dir(pattern = "FE.bw$", full.names = TRUE)
    
    #bigwigs and peaks from Terri's MCF10A histone data
    wd2 = "/slipstream/galaxy/uploads/working/qc_framework/output"
    bw_hm = (dir(wd2, full.names = TRUE, pattern = "^MCF10A.+ed$") %>% dir(pattern = "_FE.bw$", full.names = TRUE))[1:4]
    np_hm = sub("_FE.bw", "_peaks_passIDR.05.narrowPeak", bw_hm)
    stopifnot(file.exists(np_hm))
    np_gr = easyLoad_narrowPeak(np_hm)
    
    #pick 5000 most significant peaks per histone mark
    olaps_hm = ssvOverlapIntervalSets(lapply(np_gr, function(x)x[order(x$pValue, decreasing = TRUE)][1:5000]))
    
    #include Andy's higher confidence Runx1 data
    bw_af = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_RUNX1_ChIP/AF-MCF10A_RUNX1_pooled/AF-MCF10A_RUNX1_pooled_FE.bw"
    
    #clean up names such that each mark is unique and consistent
    cfg_dt = data.table(file= c(bw_book, bw_hm, bw_af))
    cfg_dt[, c("cell", "mark") := file %>% basename %>% tstrsplit(., "_", keep = 1:2)]
    cfg_dt[, mark := sub("-4336BF", "", mark)]
    cfg_dt[cell != "MCF10A", mark := paste0(mark, sub("MCF10A", "", cell))]
    cfg_dt[, cell := sub("-.+", "", cell)]
    cfg_dt[cell == "AF", c("cell", "mark") := .("MCF10A", "Runx1-AF")]
    
    #uses Josh's published peak call
    load("~/R/JR_reseq/Josh_final_peaks.save")
    olaps_runx = ssvOverlapIntervalSets(my_np)
    ssvFeatureVenn(olaps_runx)
    
    #remove histone peaks at runx sites
    olaps_hm = subsetByOverlaps(olaps_hm, olaps_runx, invert = TRUE)
    mcols(olaps_hm) = NULL
    for(mc in colnames(mcols(olaps_runx))){
        mcols(olaps_hm)[[mc]] = FALSE
    }
    
    #combine and resize
    query_gr = resize(c(olaps_runx, olaps_hm), 5000, "center")
    
    #fetch raw data to establish norm_factor
    options(mc.cores = 32)
    raw_bw = stsFetchTsneInput(cfg_dt, query_gr, cap_value = Inf)$bw_dt
    dt_quant = ssvSignalBandedQuantiles(raw_bw, by_ = "mark", return_data = TRUE)
    dt_norm = dt_quant[q_range == "90-95%", .(norm_factor = 1/max(high)), by = .(facet_group)]
    colnames(dt_norm)[1] = "mark"
    cfg_dt = merge(cfg_dt, dt_norm)
    cfg_dt = cfg_dt[, .(file, cell, mark, norm_factor)]
    
    #fetch input for tsne
    tsne_input = stsFetchTsneInput(cfg_dt, query_gr, cap_value = 1)
    profile_dt = tsne_input$bw_dt
    query_gr = tsne_input$query_gr
    
    #run tsne
    tsne_dt = stsRunTsne(profile_dt, perplexity = 300)
    
    #set color_mapping
    color_mapping = safeBrew(nrow(cfg_dt))
    names(color_mapping) = cfg_dt$mark
    
    #calculate max enrichment per region
    agg_dt = profile_dt[, .(value = max(y)), .(id, cell, mark)]
    agg_dt = merge(agg_dt, tsne_dt)
    
    overlap_dt = as.data.table(tsne_input$query_gr)
    memb = as.data.table(ssvFactorizeMembTable(overlap_dt[, .(blocked, dmso, released)]))
    memb$id = overlap_dt$id
    overlap_dt = merge(overlap_dt, memb)
    overlap_dt = overlap_dt[, .(id, blocked, dmso, released, group)]
    overlap_dt = merge(tsne_dt, overlap_dt, by = "id")
    
    ref_gr =rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", feature.type = "transcript", format = "gtf")
    ref_dist = distanceToNearest(ref_gr, query_gr)
    ref_dist = subset(ref_dist, !is.na(distance))
    annotation_dt = data.table(
        gene_name = ref_gr$gene_name[queryHits(ref_dist)],
        id = query_gr$id[subjectHits(ref_dist)],
        distance = mcols(ref_dist)$distance
    )
    
    list(profile_dt, tsne_dt, query_gr, agg_dt, overlap_dt, annotation_dt)
})

profile_dt = res[[1]]
tsne_dt = res[[2]]
query_gr = res[[3]]
agg_dt = res[[4]]
overlap_dt = res[[5]]
annotation_dt = res[[6]]
remove(res)