library(magrittr)
library(seqtsne)
library(seqsetvis)
library(GenomicRanges)
library(cowplot)
bfc = BiocFileCache::BiocFileCache()

rname = "seqtsne_example"

res = bfcif(bfc,
            force_overwrite = FALSE,
            # force_overwrite = TRUE,
            rname, function(){
                ref_gr =rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", feature.type = "transcript", format = "gtf")
                ref_dist = distanceToNearest(ref_gr, seqtsne::query_gr)
                ref_dist = subset(ref_dist, !is.na(distance))
                annotation_dt = data.table(
                    gene_name = ref_gr$gene_name[queryHits(ref_dist)],
                    id = seqtsne::query_gr$id[subjectHits(ref_dist)],
                    distance = mcols(ref_dist)$distance
                )
                
                # color_mapping = safeBrew(length(unique(sub(".+_", "", cfg_dt$mark))))
                # color_mapping = rep(color_mapping, length(unique(sub("_.+", "", cfg_dt$mark))))
                # names(color_mapping) = cfg_dt$mark
                
                prof_dt = seqtsne::profile_dt[, .(id, tall_var = "-", wide_var = paste(tall_var, wide_var), cell = tall_var, mark = wide_var, x, y)]
                tsne_res = seqtsne::stsRunTsne(prof_dt, perplexity = 5)
                # cfg_dt = unique(prof_dt[, .(wide_var, cell, mark, tall_var)])
                # cfg_dt$norm_factor = 1
                # 
                # src_dir = "/slipstream/galaxy/uploads/working/qc_framework/output_hESC_court"
                # cfg_dt[, file := file.path(src_dir, paste(cell, mark, "pooled", sep = "_"))]
                # cfg_dt[, file := file.path(file, paste(cell, mark, "pooled_FE.bw", sep = "_"))]
                # 
                cfg_dt = seqtsne::ex_cfg_dt
                # colnames(cfg_dt)[1] = "files"
                cfg_dt[, file := bw_files]
                cfg_dt = merge(cfg_dt, seqtsne::ex_cfg_dt.bam[, .(bam_files, tall_var, wide_var)])
                cfg_dt[, cell := tall_var]
                cfg_dt[, mark := wide_var]
                cfg_dt[, tall_var := "-"]
                cfg_dt[, wide_var := paste(cell, mark)]
                
                agg_dt = prof_dt[, .(value = max(y)), .(tall_var, wide_var, cell, mark, id)]
                agg_dt = merge(agg_dt, tsne_res, by = c('id', "tall_var"))
                
                res = list(profile_dt = prof_dt, 
                           tsne_dt = tsne_res, 
                           query_gr = seqtsne::query_gr, 
                           agg_dt = agg_dt, 
                           overlap_dt = NULL,#overlap_dt, 
                           annotation_dt = annotation_dt, 
                           config_dt = cfg_dt)
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
# overlap_dt = res[[5]]
annotation_dt = res[[6]]
config_dt = res[[7]]
# color_mapping = res[[8]]
remove(res)

