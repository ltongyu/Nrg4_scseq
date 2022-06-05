## Passing cluster, umap, color to velocyto
## ==================================================

## Data is quite noisy

rm(list=ls())
library(Seurat)

load(paste0("data/macro_res0.1.Rda"))

subdata <- subset(sub_cluster, idents=c("0","1","2","3","4"))

load(paste0("data/macro_anno.Rda"))

new_labels <- Idents(sub_cluster_anno)[rownames(subdata@meta.data)]
subdata@meta.data$anno <- new_labels[rownames(subdata@meta.data)]
macs_obj <- subdata

p1 	<- DimPlot(object = macs_obj , reduction = "umap", label = TRUE) 
p2 	<- DimPlot(object = macs_obj , reduction = "umap",group.by="anno" ,label = TRUE) 

p3 <- p2 + NoLegend()

library(ggplot2)


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


color_list2 <- ggplotColours(n=3)

## need to add the cell names to make sure the correct order
clust_idx <- as.numeric(macs_obj@meta.data$anno)
seurat_clusters <- cbind.data.frame(clusters = clust_idx, color = color_list2[clust_idx], anno=macs_obj@meta.data$anno)

rownames(seurat_clusters) <- rownames(macs_obj@meta.data)

write.csv(Cells(macs_obj), file = paste0("data/velo/chow_nash_mac_sub_cellID_obs.csv"), row.names = FALSE)
write.csv(Embeddings(macs_obj, reduction = "umap"), file = paste0("data/velo/chow_nash_mac_sub_cell_embeddings.csv"))
write.csv(seurat_clusters, file = paste0("data/velo/chow_nash_mac_sub_clusters_obs.csv"))





