
# https://www.biostars.org/p/450410/

rm(list=ls())
library(Seurat)
library(SeuratObject)
library(Matrix) 

load(paste0("data/new_analysis_TL_anno.Rda"))

cnt_mat <- chow_nash_jq_anno@assays$RNA@counts
gn <- rownames(cnt_mat)

require("biomaRt")
##===================================
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# save(mouse,file="data/biomart_data/mouse_mart.rds")
# save(human,file="data/biomart_data/human_mart.rds")

# genesV2 = getLDS(attributes = c("mgi_symbol"), 
#                     filters = "mgi_symbol", 
#                     values = gn , 
#                     mart = mouse, 
#                     attributesL = c("hgnc_symbol"), 
#                     martL = human, 
#                     uniqueRows=T)

# save(genesV2,file=paste0("data/raw_mgi_hgnc_nrg4.rds"))
##===================================

load(paste0("data/biomart_data/mouse_mart.rds"))
load(paste0("data/biomart_data/human_mart.rds"))
load(paste0("data/raw_mgi_hgnc_nrg4.rds"))

# genesV2[duplicated(genesV2$MGI.symbol),]

length(genesV2$MGI.symbol)


duplicated_gn <- names(table(genesV2$MGI.symbol))[which(table(genesV2$MGI.symbol)!=1)]


one2one_df <- genesV2[-which(genesV2$MGI.symbol%in%duplicated_gn),]

duplicated_raw <- genesV2[which(genesV2$MGI.symbol%in%duplicated_gn),]
duplicated_raw_order <- duplicated_raw[order(duplicated_raw$MGI.symbol),]

duplicated_order_dashfil <- duplicated_raw_order[-grep("-",duplicated_raw_order$HGNC.symbol),]


## part I to save
##=======================
reserve_dup1 <- duplicated_order_dashfil[duplicated_order_dashfil$MGI.symbol %in% names(table(duplicated_order_dashfil$MGI.symbol))[which(table(duplicated_order_dashfil$MGI.symbol)==1)],]

duplicated_order_dashfil_v2 <- duplicated_order_dashfil[-which(duplicated_order_dashfil$MGI.symbol %in% reserve_dup1$MGI.symbol),]

duplicated_order_dashfil_v2$upper <- toupper(duplicated_order_dashfil_v2$MGI.symbol)
perfect_match_idx <- apply(duplicated_order_dashfil_v2,1,function(x){identical(as.character(x[2]),as.character(x[3]))})

## part II to save
##=======================
reserve_dup2 <- duplicated_order_dashfil_v2[perfect_match_idx,1:2]

extend_map <- rbind(one2one_df,reserve_dup1,reserve_dup2)

duplicated_order_dashfil_v3 <- duplicated_order_dashfil_v2[-perfect_match_idx,]

## not consider the one to many case
##=======================
duplicated_order_dashfil_v4 <- duplicated_order_dashfil_v3[-which(duplicated_order_dashfil_v3$MGI.symbol %in% reserve_dup2$MGI.symbol),]


## scenario I: not consider the many to one case
##===================================================
duplicated_order_dashfil_v5 <- duplicated_order_dashfil_v4[-which(duplicated_order_dashfil_v4$HGNC.symbol%in%extend_map$HGNC.symbol),]
# length(unique(duplicated_order_dashfil_v5$MGI.symbol))

duplicated_order_dashfil_v6 <- duplicated_order_dashfil_v5[order(                                duplicated_order_dashfil_v5$MGI.symbol,
                                    duplicated_order_dashfil_v5$HGNC.symbol),]

duplicated_order_dashfil_v7 <- duplicated_order_dashfil_v6[!duplicated(duplicated_order_dashfil_v6$MGI.symbol),]

reserve_dup3 <- duplicated_order_dashfil_v7[!duplicated(duplicated_order_dashfil_v7$HGNC.symbol),1:2]


length(union(reserve_dup3$HGNC.symbol,extend_map$HGNC.symbol))
length(c(reserve_dup3$HGNC.symbol,extend_map$HGNC.symbol))

all_hgnc <- c(reserve_dup3$HGNC.symbol,extend_map$HGNC.symbol)


which(c(extend_map$HGNC.symbol,reserve_dup3$HGNC.symbol)=="ADRM1")
extend_map[c(148,307),]

extend_map <- rbind(extend_map,reserve_dup3)

#     MGI.symbol HGNC.symbol
# 163     Gm9774       ADRM1
# 326      Adrm1       ADRM1


# > length(unique(extend_map$MGI.symbol))
# [1] 14687
# > length(unique(extend_map$HGNC.symbol))
# [1] 14340


duplicated_hgnc <- names(table(extend_map$HGNC.symbol))[which(table(extend_map$HGNC.symbol)!=1)]

final_duplicated <- extend_map[extend_map$HGNC.symbol %in%duplicated_hgnc, ]

final_duplicated_order <- final_duplicated[order(final_duplicated$HGNC.symbol),]
final_duplicated_order$upper <- toupper(final_duplicated_order$MGI.symbol)

reserve_dup4 <- final_duplicated_order[apply(final_duplicated_order,1,function(x){identical(as.character(x[2]),as.character(x[3]))}),1:2]


final_duplicated_order_v2 <- final_duplicated_order[!(final_duplicated_order$HGNC.symbol %in% reserve_dup4$HGNC.symbol),]

final_duplicated_order_v3 <- final_duplicated_order_v2[order(final_duplicated_order_v2$HGNC.symbol,final_duplicated_order_v2$MGI.symbol),]


reserve_dup5 <- final_duplicated_order_v3[!duplicated(final_duplicated_order_v3$HGNC.symbol),1:2]


extend_map_fil  <- extend_map[!(extend_map$HGNC.symbol %in%duplicated_hgnc),]

extend_map_final <- rbind(extend_map_fil,reserve_dup4,reserve_dup5)

length(unique(extend_map_final$MGI.symbol))
length(unique(extend_map_final$HGNC.symbol))

save(extend_map_final,file=paste0("data/one2one_mgi_hgnc_nrg4.rds"))











