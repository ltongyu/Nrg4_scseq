

#### add noise function####
data_adj <- function(data_input, random_seed=1){
  set.seed(random_seed)
  data_output <- data_input
  for (i in 1: dim (data_input) [1]){
    for (j in 1: dim(data_input) [2]){
      data_output [i ,j] <- data_input [i, j] + rnorm (1, 0, 0.1)
    }
  }
  return(data_output)
}



virtual_flow_plot_2 <- function (vals_f, marker_1, marker_2, cutoff_1 , cutoff_2, geno_1, geno_2){
  y <- vals_f [, c(marker_1, marker_2, 'meta_data.geno')]
  y$geno_sort <- factor (y$meta_data.geno, levels = c(geno_1, geno_2))
  colnames(y) <- c('gene1','gene2', 'geno', 'geno_sort')

  a <- y[which(y$geno %in% geno_1), ]
  b <- length(which (a$gene1 > cutoff_1 & a$gene2 > cutoff_2))
  perc_text_WT <- paste0(geno_1,':',b,'/',dim(a)[1] , '=',round(b/dim(a)[1],2))
  a <- y[which(y$geno %in% geno_2), ]
  b <- length(which (a$gene1 > cutoff_1 & a$gene2 > cutoff_2))
  perc_text_KO <- paste0(geno_2, ':',b,'/',dim(a)[1] , '=',round(b/dim(a)[1],2))
  per_label <- c(perc_text_WT , perc_text_KO)
  names(per_label) <- c(geno_1, geno_2)
  
  p <- ggplot (y, aes (x = gene1, y =gene2)) +
    stat_density2d(geom = "tile", aes (fill = stat(count)), contour = F)+
    scale_fill_gradientn(colours=colfunc(400)) +

    geom_point(size = 0.5, alpha = 0.1,shape = 16) +
    ylab(marker_2)+ xlab(marker_1) +

    geom_vline(xintercept=cutoff_1, linetype="dashed", color = "red") +
    geom_hline(yintercept=cutoff_2, linetype="dashed", color = "red") +
    facet_wrap(vars(geno_sort), labeller = labeller(geno_sort = per_label))

  return (p)
}


