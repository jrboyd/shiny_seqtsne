library(magrittr)
library(seqtsne)
library(seqsetvis)
library(GenomicRanges)
library(cowplot)
bfc = BiocFileCache::BiocFileCache()

rname = "Bcell_biv_wide_full"

res = bfcif(bfc,
            force_overwrite = FALSE,
            rname, function(){
                # bigwigs matching Josh's published peak call
                wd = "~/tsne_data/mm10_Bcell_tsne/"
                bw_book = dir(wd, full.names = TRUE) %>% dir(., pattern = "3K(4|27)me3.+_FE.bw$", full.names = TRUE)
                stopifnot(file.exists(bw_book))
                
                np_book = dir(wd, full.names = TRUE) %>% dir(., pattern = "Peak$", full.names = TRUE)
                stopifnot(file.exists(np_book))
                
                np_gr = easyLoad_narrowPeak(np_book)
                
                
                
                #pick 5000 most significant peaks per histone mark
                # biv_gr = lapply(unique(sub("_h.+", "", names(np_gr))), function(cl){
                #     gr = ssvOverlapIntervalSets(np_gr[grepl(cl, names(np_gr))])
                #     gr[rowSums(as.data.frame(mcols(gr))) == 2]
                # })
                # sapply(
                # olaps_hm = ssvOverlapIntervalSets(lapply(np_gr, function(x)x[order(x$pValue, decreasing = TRUE)][1:5000]))
                olaps_hm = ssvOverlapIntervalSets(np_gr)
                # ssvFeatureVenn(olaps_hm)
                #clean up names such that each mark is unique and consistent
                cfg_dt = data.table(file= bw_book)
                cfg_dt[, c("cell", "mark") := file %>% basename %>% tstrsplit(., "_", keep = c(2, 1))]
                cfg_dt[, wide_var := paste(cell, mark)]
                cfg_dt[, tall_var := "-"]
                
                
                #resize
                query_gr = resize(olaps_hm, 2000, "center")
                
                #fetch raw data to establish norm_factor
                options(mc.cores = 32)
                raw_bw = stsFetchTsneInput(cfg_dt, query_gr, cap_value = Inf)$bw_dt
                dt_quant = ssvSignalBandedQuantiles(raw_bw, by_ = c("wide_var"), return_data = TRUE)
                dt_norm = dt_quant[q_range == "90-95%", .(norm_factor = 1/max(high)), by = .(facet_group)]
                colnames(dt_norm)[1] = "wide_var"
                cfg_dt = merge(cfg_dt, dt_norm)
                # cfg_dt = cfg_dt[, .(file, tall_var, wide_var, norm_factor)]
                
                #fetch input for tsne
                tsne_input = stsFetchTsneInput(cfg_dt, query_gr, cap_value = 1)
                profile_dt = tsne_input$bw_dt
                query_gr = tsne_input$query_gr
                
                #run tsne
                tsne_dt = stsRunTsne(profile_dt, perplexity = 300)
                
                #set color_mapping
                # color_mapping = safeBrew(length(unique(cfg_dt$mark)))
                # names(color_mapping) = unique(cfg_dt$mark)
                color_mapping = col2hex(c("firebrick1", "forestgreen"))
                names(color_mapping) = c("H3K27me3", "H3K4me3")
                
                cm_alt = unique(cfg_dt$cell)
                color_mapping.alt = safeBrew(n = length(cm_alt))
                names(color_mapping.alt) = cm_alt
                
                #calculate max enrichment per region
                agg_dt = profile_dt[, .(value = max(y)), .(id, wide_var, tall_var)]
                agg_dt = merge(agg_dt, tsne_dt)
                
                overlap_dt = as.data.table(tsne_input$query_gr)
                cn = colnames(mcols(tsne_input$query_gr))
                cn = cn[cn != "id"]
                memb = as.data.table(ssvFactorizeMembTable(overlap_dt[, cn, with = FALSE]))
                memb$id = overlap_dt$id
                overlap_dt = merge(overlap_dt, memb, by = "id")
                overlap_dt = overlap_dt[, c("id", cn, "group"), with = FALSE]
                overlap_dt = merge(tsne_dt, overlap_dt, by = "id")
                
                ref_gr =rtracklayer::import.gff(
                    "~/gencode.vM20.annotation.gtf.gz", 
                    feature.type = "transcript", format = "gtf")
                ref_dist = distanceToNearest(ref_gr, query_gr)
                ref_dist = subset(ref_dist, !is.na(distance))
                annotation_dt = data.table(
                    gene_name = ref_gr$gene_name[queryHits(ref_dist)],
                    id = query_gr$id[subjectHits(ref_dist)],
                    distance = mcols(ref_dist)$distance
                )
                
                # color_mapping = safeBrew(length(unique(sub(".+_", "", cfg_dt$mark))))
                # color_mapping = rep(color_mapping, length(unique(sub("_.+", "", cfg_dt$mark))))
                # names(color_mapping) = cfg_dt$mark
                
                res = list(profile_dt = profile_dt, 
                     tsne_dt = tsne_dt, 
                     query_gr = query_gr, 
                     agg_dt = agg_dt, 
                     overlap_dt = overlap_dt, 
                     annotation_dt = annotation_dt, 
                     config_dt = cfg_dt, 
                     color_mapping = color_mapping,
                     color_mapping.alt = color_mapping.alt)
                res
            })

if(FALSE){
    res = bfcif(bfc,
                force_overwrite = TRUE,
                rname, function(){
                    res
                })    
}


profile_dt = res[[1]]
tsne_dt = res[[2]]
query_gr = res[[3]]
agg_dt = res[[4]]
overlap_dt = res[[5]]
annotation_dt = res[[6]]
config_dt = res[[7]]
color_mapping = res[[8]]
# color_mapping.alt = res[[9]]

# cm_alt = unique(config_dt$cell)
# color_mapping.alt = safeBrew(n = length(cm_alt))
# names(color_mapping.alt) = cm_alt
# 
# cm_full = unique(paste(config_dt$cell,config_dt$mark))
# color_mapping.full = safeBrew(n = length(cm_full))
# names(color_mapping.full) = cm_full

remove(res)

