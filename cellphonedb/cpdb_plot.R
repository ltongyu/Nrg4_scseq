rm(list=ls())
library(ggplot2)

means_separator = "\t"
pvalues_separator = '\t'
selected_rows = NULL
selected_columns = c("CD8|Macrophage","Macrophage|CD8")

res_list <- list()
#
for(iset in c("chow","nash")){

  if(iset == "combined"){
    all_pval = read.table(paste0("output/cpdb/pvals_stat_cd_mac_10k.txt"), header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
    all_means = read.table(paste0("output/cpdb/means_stat_cd_mac_10k.txt"), header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  }else{
    all_pval = read.table(paste0("output/cpdb/pvals_stat_cd_mac_",iset,"_10k.txt"), header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
    all_means = read.table(paste0("output/cpdb/means_stat_cd_mac_",iset,"_10k.txt"), header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  }


  intr_pairs = all_pval$interacting_pair
  all_pval  = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  sel_pval  = all_pval[,selected_columns]
  sel_means = all_means[,selected_columns]

  rownames(sel_pval) <- rownames(sel_means) <- intr_pairs
  colnames(sel_pval) <- colnames(sel_means) <- paste0(selected_columns,"|",iset)


  res_list[[iset]] <- list(pval=sel_pval,means=sel_means)
  rm(sel_pval,sel_means,all_means,all_pval,intr_pairs)
}


pval_merge_two <- merge(res_list[[1]]$pval,res_list[[2]]$pval,by="row.names",all=TRUE)
rownames(pval_merge_two) <- pval_merge_two$Row.names

means_merge_two <- merge(res_list[[1]]$means,res_list[[2]]$means,by="row.names",all=TRUE)
rownames(means_merge_two) <- means_merge_two$Row.names

rownames(means_merge_two) <- means_merge_two[,1]
rownames(pval_merge_two) <- pval_merge_two[,1]

identical(rownames(pval_merge_two),rownames(means_merge_two))
identical(colnames(pval_merge_two),colnames(means_merge_two))


pval_merge_two        <- pval_merge_two[,-1]
means_merge_two       <- means_merge_two[,-1]


idirect <- "cd2mac"
if(idirect=="cd2mac"){
  pval_sub  <- pval_merge_two[,c(1,3)]
  means_sub <- means_merge_two[,c(1,3)]
}else{
  pval_sub  <- pval_merge_two[,c(2,4)]
  means_sub <- means_merge_two[,c(2,4)]
}

rm(pval_merge_two,means_merge_two)


selected_rows         <- rownames(pval_sub)
new_selected_columns  <- colnames(pval_sub)


kept_idx        <- which(apply(pval_sub,1,function(x){sum(x<0.05,na.rm=TRUE)})!=0)
sel_pval_kept   <- pval_sub[kept_idx,]
sel_means_kept  <- means_sub[kept_idx,]
selected_rows_kept <- selected_rows[kept_idx]

df_names = expand.grid(selected_rows_kept, new_selected_columns)
# pval = unlist(sel_pval)
pval = unlist(sel_pval_kept)
pval[pval==0] = 0.0009
plot.data = cbind(df_names,pval)
pr = unlist(as.data.frame(sel_means_kept))
pr[pr==0] = 1
plot.data = cbind(plot.data,log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

plot.data$source  <- sapply(strsplit(as.character(plot.data$clusters),split="\\|"),"[[",1)
plot.data$target  <- sapply(strsplit(as.character(plot.data$clusters),split="\\|"),"[[",2)
plot.data$Geno    <- sapply(strsplit(as.character(plot.data$clusters),split="\\|"),"[[",3)

plot.data$shared_clusters <- paste0(plot.data$source ,"|",plot.data$target )

color.heatmap = "YlGnBu"
direction = -1
n.colors = 9
if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
        RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
        (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
}else{
    color.use <- color.heatmap
}

if (direction == -1) {
    color.use <- rev(color.use)
}

# my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
my_palette <-  colorRampPalette(color.use)(99)


g1 <- ggplot(plot.data,aes(x=clusters,y=pair)) +
geom_point(aes(size=-log10(pvalue),color=mean)) +
scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text= element_blank(),
      # axis.text=element_text(size=14, colour = "black"),
      # axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
      # axis.text.y = element_text(size=12, colour = "black"),
      axis.title=element_blank(),
      panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
      legend.position="none")

