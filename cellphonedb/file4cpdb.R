
rm(list=ls())
library(Seurat)
library(SeuratObject)
library(Matrix) 


load(paste0('data/new_analysis_TL_anno_JZ_detailed.Rda'))
Idents(chow_nash_jq_anno) <- chow_nash_jq_anno@meta.data$detailed_anno

chow_nash_jq_anno@meta.data$mac_anno <- chow_nash_jq_anno@meta.data$detailed_anno
chow_nash_jq_anno@meta.data$mac_anno[chow_nash_jq_anno@meta.data$detailed_anno%in% c("KC","MDM","NAM")] <- "Macrophage"

Idents(chow_nash_jq_anno) <- chow_nash_jq_anno@meta.data$mac_anno

CD8_NAM   <- subset(chow_nash_jq_anno, subset=mac_anno %in% c("CD8","Macrophage"))

rm(chow_nash_jq_anno)
cnt_mat <- CD8_NAM@assays$RNA@counts
gn <- rownames(cnt_mat)

load(paste0("data/one2one_mgi_hgnc_nrg4.rds"))

match_one2one <- extend_map_final[extend_map_final$MGI.symbol %in% gn,]
rownames(match_one2one) <- match_one2one$MGI.symbol

fil_cnt_mat <- cnt_mat[match_one2one$MGI.symbol,]
select_gn <- rownames(fil_cnt_mat)
rownames(fil_cnt_mat) <- match_one2one[select_gn,"HGNC.symbol"]

writeMM(fil_cnt_mat, file = paste0('data/cpdb/cd_mac/matrix.mtx'))
write(x = rownames(fil_cnt_mat), file = paste0('data/cpdb/cd_mac/features.tsv'))
write(x = colnames(fil_cnt_mat), file = paste0('data/cpdb/cd_mac/barcodes.tsv'))

CD8_NAM@meta.data$Cell = rownames(CD8_NAM@meta.data)
CD8_NAM@meta.data$cell_type = Idents(CD8_NAM)
df = CD8_NAM@meta.data[, c('Cell', 'cell_type')]

write.table(df, file =paste0('data/cpdb/CD8_MAC_meta_hgnc.tsv'), sep = '\t', quote = F, row.names = F)

