load("~/R/waldron/waldron_tsne_v2.save")

profile_dt = prof_dt[, .(tall_var, wide_var, id, x, y, cell = tall_var, mark = wide_var)]
tsne_dt = tsne_res
query_gr = query_gr
agg_dt = prof_dt[, .(value = max(y[abs(x) < .2])), .(tall_var, wide_var, id)]
# agg_dt[, value := y]
agg_dt = merge(agg_dt, tsne_dt[, .(id, tx, ty)], by = "id")
# overlap_dt = res[[5]]
annotation_dt = data.table(id = query_gr$id, gene_name = query_gr$gene_name, distance = 0)
config_dt = qdt[, .(file = file.path("/slipstream/galaxy/uploads/working/qc_framework", file), tall_var = "hESC", wide_var = mark, cell, mark)]
# color_mapping = res[[8]]
# remove(res)
