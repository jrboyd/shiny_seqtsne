source("generate_data.R")
ggplot(agg_dt, aes(x = tx, y = ty, color = value)) + facet_wrap("mark") + geom_point(size = .2)
ggsave("agg_plot.pdf", width = 10, height = 10)

#seperate summaries per mark
summary_plots = lapply(cfg_dt$mark, function(m){
    seqtsne::stsPlotSummaryProfiles(profile_dt[mark == m], tsne_res, 8, 
                                    line_color_mapping = color_mapping, ylim = c(0,1)) + 
        labs(title = m) + theme_classic()
})
pg = cowplot::plot_grid(plotlist = summary_plots)
ggsave("summary_plots.pdf", pg, width = 10, height = 10)

#runx1 peak call defined groups
ggplot(anno_dt, aes(x = tx, y = ty)) + 
    geom_point(alpha = 1, shape = 16, size = .4, color = "gray") + theme_void() +
    coord_cartesian(expand = FALSE)
ggsave("olaps_bg.png", width = 5, height = 5)
dt = data.table(file = "olaps_bg.png")
ggplot() + 
    geom_image.rect(data = dt, aes(image = "olaps_bg.png", xmin = -.5, xmax = .5, ymin = -.5, ymax = .5)) +
    geom_point(data = anno_dt, 
               aes(x = tx, y = ty, color = group), 
               alpha = 1, shape = 16, size = .4) + 
    guides(color = "none") +
    facet_wrap("group")
ggsave("olaps_plot.pdf", width = 10, height = 10)

#use cell to code runx overlap groups
group_dt = merge(profile_dt, anno_dt[, .(id, group)])
group_dt[, cell := group]
tgroup_dt = merge(tsne_res, anno_dt[, .(id, group)], by = "id")
tgroup_dt[, cell := group]

x_points = y_points = 12

#summaries grouped by mark and then overlaps
summary_plots2 = lapply(cfg_dt$mark, function(m){
    img_res = seqtsne::stsPlotSummaryProfiles(group_dt[mark == m], tgroup_dt, x_points, force_rewrite = TRUE, min_size = 0,
                                              line_color_mapping = color_mapping, ylim = c(0,1),
                                              facet_by = "cell", return_data = TRUE)  
    img_res[, img_size := N / max(N), by = .(cell)]
    xrng = c(-.5, .5)
    yrng = c(-.5, .5)
    
    xspc = diff(xrng) / x_points / 2
    yspc = diff(yrng) / y_points / 2
    
    img_res[, xmin := tx - xspc * img_size]
    img_res[, xmax := tx + xspc * img_size]
    img_res[, ymin := ty - yspc * img_size]
    img_res[, ymax := ty + yspc * img_size]
    
    # plot_tsne_img_byCell(img_res, 8)    
    ggplot() + geom_image.rect(data = img_res[img_size > .3], aes(image = png_file, 
                                                                  xmin = xmin, xmax = xmax, 
                                                                  ymin = ymin, ymax = ymax)) + facet_wrap("cell") + 
        coord_cartesian(xlim = c(-.5, .5), ylim = c(-.5, .5)) +
        labs(title = m) + theme_classic() + 
        theme(strip.text = element_text(size = 6), 
              axis.text = element_text(size = 6))
        
})
pg = cowplot::plot_grid(plotlist = summary_plots2)
ggsave("summary_plots2.pdf", pg, width = 10, height = 10)

#summaries grouped by overlaps and then mark
grp = levels(group_dt$group)[1]
m = cfg_dt$mark[1]
summary_plots3 = lapply(as.character(unique(group_dt$group)), function(grp){
    img_res = lapply(cfg_dt$mark, function(m){
        tmp = seqtsne::stsPlotSummaryProfiles(group_dt[mark == m & cell == grp], tgroup_dt, x_points, force_rewrite = TRUE, min_size = 0,
                                              line_color_mapping = color_mapping, ylim = c(0,1),
                                              facet_by = "cell", return_data = TRUE)  
        tmp$mark = m
        tmp
    }) %>% rbindlist
    img_res[, img_size := N / max(N), by = .(cell)]
    xrng = c(-.5, .5)
    yrng = c(-.5, .5)
    
    xspc = diff(xrng) / x_points / 2
    yspc = diff(yrng) / y_points / 2
    
    img_res[, xmin := tx - xspc * img_size]
    img_res[, xmax := tx + xspc * img_size]
    img_res[, ymin := ty - yspc * img_size]
    img_res[, ymax := ty + yspc * img_size]
    
    # plot_tsne_img_byCell(img_res, 8)    
    ggplot() + geom_image.rect(data = img_res[img_size > .3],
                               aes(
                                   image = png_file,
                                   xmin = xmin,
                                   xmax = xmax,
                                   ymin = ymin,
                                   ymax = ymax
                               )) +
        coord_cartesian(xlim = c(-.5, .5), ylim = c(-.5, .5)) +
        facet_wrap("mark") +
        labs(title = grp) + 
        theme_classic() + 
        theme(strip.text = element_text(size = 6),
              axis.text = element_text(size = 6))
})
pg = cowplot::plot_grid(plotlist = summary_plots3)
ggsave("summary_plots3.pdf", pg, width = 10, height = 10)
