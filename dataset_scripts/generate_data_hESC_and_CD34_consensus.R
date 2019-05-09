library(BiocFileCache)
bfc = BiocFileCache::BiocFileCache()

# rname = "waldron_tsne_wide"
rname = "waldron_tsne_gene_driven"

bfc = BiocFileCache::BiocFileCache()
load(bfcrpath(bfc, rname))

profile_dt = res[[1]]
tsne_dt = res[[2]]
query_gr = res[[3]]
agg_dt = res[[4]]
# overlap_dt = res[[5]]
annotation_dt = res[[6]]
config_dt = res[[7]]
color_mapping = c("H3K27me3" = "firebrick", "H3K4me3" = "forestgreen")#res[[8]]
color_mapping.alt = res[[9]]
remove(res)

agg_dt = merge(agg_dt, tsne_dt)
agg_dt[, value := y]
config_dt[, c("cell", "mark") := tstrsplit(basename(file), "_", keep = 1:2)]
config_dt[, bam_file := sub("_FE.bw", ".bam", file)]
# browser()

# sapply(config_dt$bam_file, seqsetvis:::getReadLength, query_gr = qgr)
