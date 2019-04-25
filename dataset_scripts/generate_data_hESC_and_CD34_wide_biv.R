library(magrittr)
library(seqtsne)
library(seqsetvis)
library(GenomicRanges)
library(cowplot)
bfc = BiocFileCache::BiocFileCache()

res = bfcif(bfc,
            force_overwrite = FALSE,
            "hESC_CD34_seqtsne_wide_biv", function(){
                # bigwigs matching Josh's published peak call
                wd = "/slipstream/galaxy/uploads/working/qc_framework/output_hESC_court"
                bw_book = dir(wd, full.names = TRUE, pattern = "ed$") %>% dir(pattern = "FE.bw$", full.names = TRUE)
                np_book = bw_book %>%
                    sub("_H3K27me3_pooled_FE.bw", "_H3K27me3_pooled_peaks.broadPeak", .) %>%
                    sub("_H3K4me3_pooled_FE.bw", "_H3K4me3_pooled_peaks.narrowPeak", .)
                stopifnot(file.exists(np_book))

                #bigwigs and peaks from Terri's MCF10A histone data
                wd2 = "/slipstream/galaxy/uploads/working/qc_framework/output_waldron_bivalency"
                bw_hm = (dir(wd2, full.names = TRUE, pattern = "^CD.+ed$") %>% dir(pattern = "_FE.bw$", full.names = TRUE))
                np_hm = bw_hm %>%
                    sub("_H3K27me3_pooled_FE.bw", "_H3K27me3_pooled_peaks.broadPeak", .) %>%
                    sub("_H3K4me3_pooled_FE.bw", "_H3K4me3_pooled_peaks.narrowPeak", .)
                stopifnot(file.exists(np_hm))

                np_f = c(np_hm, np_book)
                np_f[grepl("K4", np_f)]
                np_gr = c(easyLoad_narrowPeak(np_f[grepl("K4", np_f)]), easyLoad_broadPeak(np_f[grepl("K27", np_f)]))

                #pick 5000 most significant peaks per histone mark
                biv_gr = lapply(unique(sub("_.+", "_", names(np_gr))), function(cl){
                    gr = ssvOverlapIntervalSets(np_gr[grepl(cl, names(np_gr))])
                    gr[rowSums(as.data.frame(mcols(gr))) == 2]
                })
                # sapply(
                # olaps_hm = ssvOverlapIntervalSets(lapply(np_gr, function(x)x[order(x$pValue, decreasing = TRUE)][1:5000]))
                olaps_hm = ssvOverlapIntervalSets(biv_gr)
                #clean up names such that each mark is unique and consistent
                cfg_dt = data.table(file= c(bw_book, bw_hm))
                cfg_dt[, c("cell", "mark") := file %>% basename %>% tstrsplit(., "_", keep = 1:2)]
                cfg_dt[, wide_var := paste(cell, mark, sep = "_")]
                cfg_dt[, tall_var := "-"]


                #resize
                query_gr = resize(olaps_hm, 1000, "center")

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

                #calculate max enrichment per region
                agg_dt = profile_dt[, .(value = max(y)), .(id, cell, mark)]
                agg_dt = merge(agg_dt, tsne_dt)

                overlap_dt = as.data.table(tsne_input$query_gr)
                cn = colnames(mcols(tsne_input$query_gr))
                cn = cn[cn != "id"]
                memb = as.data.table(ssvFactorizeMembTable(overlap_dt[, cn, with = FALSE]))
                memb$id = overlap_dt$id
                overlap_dt = merge(overlap_dt, memb, by = "id")
                overlap_dt = overlap_dt[, c("id", cn, "group"), with = FALSE]
                overlap_dt = merge(tsne_dt, overlap_dt, by = "id")

                ref_gr =rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", feature.type = "transcript", format = "gtf")
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

                list(profile_dt, tsne_dt, query_gr, agg_dt, overlap_dt, annotation_dt, cfg_dt, color_mapping)
                # res
            })

profile_dt = res[[1]]
tsne_dt = res[[2]]
query_gr = res[[3]]
agg_dt = res[[4]]
overlap_dt = res[[5]]
annotation_dt = res[[6]]
config_dt = res[[7]]
color_mapping = res[[8]]
remove(res)

