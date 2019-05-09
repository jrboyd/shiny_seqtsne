library(seqtsne)
library(magrittr)
library(parallel)
library(GenomicRanges)
# source("dataset_scripts/generate_example.R")
source("dataset_scripts/generate_data_hESC_and_CD34_consensus.R")


query_gr
names(query_gr) = query_gr$id

# qgr = subset(query_gr, group == "bivalent")[1:3]
qgr = query_gr[1:3]

options(mc.cores = 10)

rl = seqsetvis:::getReadLength(config_dt$bam_file[1], query_gr = qgr)

shift_dt = mclapply(config_dt$bam_file, function(f){
    shift_dt = seqsetvis::crossCorrByRle(flip_strand = FALSE,
                                         f,
                                         query_gr = resize(qgr, 500, fix = "center"),
                                         fragment_sizes = seq(100, 500, by = 5),
                                         read_length = rl
    )
    shift_dt$bam_file = f
    shift_dt = merge(shift_dt, config_dt[, .(bam_file, cell, mark)], by = "bam_file")
    shift_dt$file = NULL
    shift_dt
}) %>% rbindlist


# sapply(config_dt$bam_file, seqsetvis:::getReadLength, query_gr = qgr)
# ggplot(shift_dt, aes(
#     x = shift,
#     y = correlation,
#     color = group,
#     group = paste(cell, group)
# )) + geom_path() + facet_grid(mark+cell~id)
shift_max = shift_dt[, .(fragLen = shift[which.max(correlation)], correlation = max(correlation)), by = .(id, cell, mark)]

ggplot() + 
    geom_path(data = shift_dt[mark == "H3K4me3" & id == "1" & grepl("HUES", cell)], 
              aes(
                  x = shift,
                  y = correlation,
                  group = cell
              )) + 
    geom_point(data = shift_max[mark == "H3K4me3" & id == "1" & grepl("HUES", cell)], 
               aes(x = fragLen, y = correlation)) + 
    coord_cartesian(ylim = c(0, 1)) +
    facet_grid(mark+cell~id)

step = 10
prof_dt = stsFetchTsneInput(config_dt[, .(file = bam_file, wide_var = cell, tall_var = mark, cell, mark)], 
                            qwin = step/2, 
                            qmet = "sample",
                            qgr = resize(qgr, 1000, fix = "center"), 
                            cap_value = Inf, high_on_right = FALSE,
                            fetch_FUN = fetch_bam_stranded_dt)$bw_dt
# prof_dt[, strand := ifelse(strand == "+", "-", "+")]

ggplot(prof_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_grid(mark+cell~id)


qshift = seq(0, 500/2, step/2)
library(data.table)
off_dt = dcast(prof_dt[, .(cell, mark, id, x, y, strand)], cell+mark+id+x~strand, value.var = "y")
off_dt[, delta := `+` - `-`]
ggplot(off_dt, aes(x = x, y = delta, color = cell)) + geom_path() + facet_grid(mark~id)

all_shift = mclapply(qshift, function(s){
    off_dt = dcast(prof_dt[, .(cell, mark, id, x = ifelse(strand == "+", x + s, x - s), y, strand)], 
                   cell+mark+id+x~strand, value.var = "y")
    off_dt[, delta := `+` - `-`]
    # off_dt[!is.na(delta)]
    off_dt$shift = s
    off_dt
}) %>% rbindlist()




ggplot() +
    geom_raster(data = all_shift[mark == "H3K4me3" & id == "1" & grepl("HUES", cell)], 
                aes(x = x, y = shift * 2, fill = delta)) +
    geom_segment(data = shift_max[mark == "H3K4me3" & id == "1" & grepl("HUES", cell)],
                 aes(
                     x = min(all_shift$x),
                     xend = max(all_shift$x),
                     y = 150,
                     yend = 150
                 ),
                 color = "red") +
    facet_grid(cell + mark ~ id) +
    annotate("segment", x = 200, xend = 200, y = 0, yend = 500) +
    scale_fill_viridis_c()

library(seqsetvis)
bam_dt = ssvFetchBam(
    config_dt[, .(file = bam_file, cell, mark)],
    qgr = qgr,
    target_strand = "both",
    fragLens = 150,
    return_data.table = TRUE,
    win_size = 10
)
ggplot(bam_dt[mark == "H3K4me3" & id == "1" & grepl("HUES", cell)], 
       aes(x = x, y = y, color = strand)) + 
    geom_path() + 
    facet_grid(cell + mark ~ id)
