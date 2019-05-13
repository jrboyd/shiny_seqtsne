library(magrittr)
library(seqtsne)
library(seqsetvis)
library(GenomicRanges)
library(cowplot)
bfc = BiocFileCache::BiocFileCache()

# rname = "MSC_biv_wide_full"
rname = "MSC_biv2_gene_driven"

res = bfcif(bfc,
            force_overwrite = FALSE,
            rname, function(){
                # bigwigs matching Josh's published peak call
                wd = "~/tsne_data/mm10_bmsc/"
                # bw_book = c(
                #     dir(wd, pattern = "3k(4|27)me3.+(d00|d21).+fault_FE.bw$", full.names = TRUE),
                #     dir(wd, pattern = "3k9(me3|ac).+(d00|d21).+fault_FE.bw$", full.names = TRUE)
                # )
                
                bw_book = c(
                    dir(wd, pattern = "(k4me3|k27me3).+(d00|d21).+fault_FE.bw$", full.names = TRUE)
                )
                            
                stopifnot(file.exists(bw_book))
                bb_book = sub("_default_FE.bw", "_default_peaks_narrow.bb", bw_book)
                stopifnot(file.exists(bb_book))
                
                
                np_gr = lapply(bb_book, function(bb){
                    cmd = paste("bigBedToBed", bb, "/dev/stdout")
                    dt = data.table(str =  system(cmd, intern = TRUE))
                    dt = dt[, tstrsplit(str, "\t")]
                    dt[, V4 := sub(".+_peak_", "peak_", V4)]
                    colnames(dt) = c("seqnames", "start", "end", "id", "width", 
                                     "strand", "signalValue", "pValue", 
                                     "qValue", "relSummit")
                    GRanges(dt)
                })
                names(np_gr) = 
                    basename(bb_book) %>% 
                    strsplit(., "_") %>% 
                    sapply(., function(x)
                        paste(x[2], x[3], x[1], sep = "_"))
                np_gr
                gr = reduce(unlist(GRangesList(np_gr)))
                gr = subset(gr, seqnames != "chrM")
                
                ref_gr = rtracklayer::import.gff("~/gencode.vM20.annotation.gtf.gz", 
                                                 feature.type = "transcript", format = "gtf")
                pr_gr = promoters(ref_gr, 1, 1)
                hit_gr = subsetByOverlaps(pr_gr, gr)
                hit_gr = split(hit_gr, hit_gr$gene_name)
                
                qgr = resize(unlist(reduce(resize(hit_gr, 200, fix = "center"))), 3000, fix = "center")
                dt_id = data.table(gene_name = names(qgr))
                dt_id[, n := seq(.N), gene_name]
                dt_id[, id := paste0(gene_name, "_", n)]
                
                stopifnot(all(dt_id$gene_name == names(qgr)))
                
                
                
                qgr$id = dt_id$id
                
                #check overlap rate
                sapply(np_gr, function(gr)length(subsetByOverlaps(gr, qgr)) / length(gr)) %>% round(., 3) *100
                
                query_gr = qgr
                
                #pick 5000 most significant peaks per histone mark
                # biv_gr = lapply(unique(sub("_h.+", "", names(np_gr))), function(cl){
                #     gr = ssvOverlapIntervalSets(np_gr[grepl(cl, names(np_gr))])
                #     gr[rowSums(as.data.frame(mcols(gr))) == 2]
                # })
                # sapply(
                # olaps_hm = ssvOverlapIntervalSets(lapply(np_gr, function(x)x[order(x$pValue, decreasing = TRUE)][1:5000]))
                # olaps_hm = ssvOverlapIntervalSets(np_gr)
                # ssvFeatureVenn(olaps_hm)
                #clean up names such that each mark is unique and consistent
                cfg_dt = data.table(file= bw_book)
                cfg_dt[, c("cell", "mark") := file %>% basename %>% tstrsplit(., "_", keep = c(3, 1))]
                cfg_dt[, wide_var := paste(cell, mark)]
                cfg_dt[, tall_var := "-"]
                
                
                #resize
                # query_gr = resize(olaps_hm, 3000, "center")
                
                #fetch raw data to establish norm_factor
                options(mc.cores = 32)
                raw_bw = stsFetchTsneInput(cfg_dt, query_gr, cap_value = Inf)$bw_dt
                dt_quant = ssvSignalBandedQuantiles(raw_bw, by_ = c("wide_var"), return_data = TRUE)
                dt_norm = dt_quant[q_range == "90-95%", .(norm_factor = 1/max(high)), by = .(facet_group)]
                colnames(dt_norm)[1] = "wide_var"
                cfg_dt = merge(cfg_dt, dt_norm)
                # cfg_dt = cfg_dt[, .(file, tall_var, wide_var, norm_factor)]
                cfg_dt
                
                cfg_dt$wide_var = factor(cfg_dt$wide_var, levels = cfg_dt[order(mark)]$wide_var)
                ggplot(cfg_dt, aes(x = wide_var, y = norm_factor, fill = mark)) + 
                    geom_bar(stat = "identity") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                max_norm = 1/20
                for(m in unique(cfg_dt$mark)){
                    penalty = cfg_dt[mark == m, max(norm_factor)] / max_norm
                    if(penalty > 1){
                        cfg_dt[mark == m, norm_factor := norm_factor/penalty]
                    }
                }
                
                # cfg_dt[, norm_factor := mean(norm_factor), by = .(mark)]
                ggplot(cfg_dt, aes(x = wide_var, y = norm_factor, fill = mark)) + 
                    geom_bar(stat = "identity") + 
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                #fetch input for tsne
                tsne_input = stsFetchTsneInput(cfg_dt, query_gr, cap_value = 1)
                profile_dt = tsne_input$bw_dt
                query_gr = tsne_input$query_gr
                
                #run tsne
                tsne_dt = stsRunTsne(profile_dt, perplexity = 300)
                
                #set color_mapping
                # color_mapping = safeBrew(length(unique(cfg_dt$mark)))
                # names(color_mapping) = unique(cfg_dt$mark)
                # color_mapping = col2hex(c("firebrick1", "forestgreen"))
                # names(color_mapping) = c("h3k27me3", "h3k4me3")
                color_mapping = NULL
                #calculate max enrichment per region
                agg_dt = profile_dt[, .(value = max(y)), .(id, wide_var, tall_var)]
                agg_dt = merge(agg_dt, tsne_dt)
                
                # overlap_dt = as.data.table(tsne_input$query_gr)
                # cn = colnames(mcols(tsne_input$query_gr))
                # cn = cn[cn != "id"]
                # memb = as.data.table(ssvFactorizeMembTable(overlap_dt[, cn, with = FALSE]))
                # memb$id = overlap_dt$id
                # overlap_dt = merge(overlap_dt, memb, by = "id")
                # overlap_dt = overlap_dt[, c("id", cn, "group"), with = FALSE]
                # overlap_dt = merge(tsne_dt, overlap_dt, by = "id")
                
                # ref_gr =rtracklayer::import.gff(
                #     "~/gencode.vM20.annotation.gtf.gz", 
                #     feature.type = "transcript", format = "gtf")
                # ref_dist = distanceToNearest(ref_gr, query_gr)
                # ref_dist = subset(ref_dist, !is.na(distance))
                # annotation_dt = data.table(
                #     gene_name = ref_gr$gene_name[queryHits(ref_dist)],
                #     id = query_gr$id[subjectHits(ref_dist)],
                #     distance = mcols(ref_dist)$distance
                # )
                
                overlap_dt = NULL
                annotation_dt = dt_id
                annotation_dt$distance = 0
                
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
                           color_mapping = color_mapping)
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
cm_alt = unique(config_dt$cell)
color_mapping.alt = safeBrew(n = length(cm_alt))
names(color_mapping.alt) = cm_alt
remove(res)

