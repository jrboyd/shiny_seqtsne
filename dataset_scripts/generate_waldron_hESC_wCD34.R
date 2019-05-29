load("~/R/waldron/waldron_tsne_v2.save")
load("~/R/waldron/waldron_tsne2_wCD34_v1.save")
library(data.table)
prof_dt_wCD34[, wide_var := paste(tall_var, wide_var)]
prof_dt_wCD34[, tall_var := "."]
profile_dt = prof_dt_wCD34[, .(tall_var, wide_var, id, x, y, cell = tall_var, mark = wide_var)]
tsne_dt = tsne_dt2
tsne_dt[, tall_var := "."]
query_gr = query_gr
agg_dt = prof_dt_wCD34[, .(value = max(y[abs(x) < .2])), .(tall_var, wide_var, id)]
# agg_dt[, value := y]
agg_dt = merge(agg_dt, tsne_dt[, .(id, tx, ty)], by = "id")
# overlap_dt = res[[5]]
annotation_dt = data.table(id = query_gr$id, gene_name = query_gr$gene_name, distance = 0)
config_dt = qdt[, .(file = file.path("/slipstream/galaxy/uploads/working/qc_framework", file), tall_var = ".", wide_var = paste(cell, mark), cell, mark)]
# color_mapping = res[[8]]
# remove(res)
