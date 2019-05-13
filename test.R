library(ggplot2)
library(seqsetvis)

qid = annotation_dt[gene_name %in% c("GATA6", "SUPT3H", "RUNX1", "RUNX2")]$id

pdt = profile_dt[id %in% qid]
ggplot(pdt, aes(x = x, y = y, color = mark)) + facet_grid("cell~id") + geom_path()


qgr = subset(query_gr, id %in% qid)
names(qgr) = qgr$id
qgr = resize(qgr, 30000, fix = "center")
bw_dt = ssvFetchBigwig(config_dt, qgr, return_data.table = TRUE, win_method = "summary", win_size = 60)
# ggplot(bw_dt, aes(x = x, y = y+20, color = mark, group = sample)) + facet_grid("cell~id") + geom_path() + scale_y_log10(labels = function(x)x-20)
pdt2 = bw_dt[, .(y = mean(y)), by = .(x, cell, mark, id)]
ggplot(pdt2, aes(x = x, y = y+20, color = mark)) + facet_grid("cell~id") + geom_path() + scale_y_log10(labels = function(x)x-20)

for(m in as.data.table(qgr)[, paste(seqnames, start, end)]) message(m)
