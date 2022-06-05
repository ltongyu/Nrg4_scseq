rm(list=ls())

library(Matrix)
library(Seurat)

# datadir1 = "/net/mulan/disk2/jiaqiang/linlab/zmchen/11112020/"
datadir2 = "/net/mulan/disk2/jiaqiang/linlab/tongyu/"
workdir = "/net/mulan/disk2/jiaqiang/linlab/zmchen/2021/03312021/"
# load(paste0(datadir1,"/data/clean_data_dbl_ver4.rds"))

fn <- list.files(datadir2,"gz")

ty_data <- list()


library(data.table)


for(ifn in 1:length(fn)){
	dt <-  as.data.frame(fread(paste0(datadir2,fn[ifn])))
	count_file <- dt[,-1]
	rownames(count_file) <- dt[,1]
	cn <- paste0(colnames(count_file),"_",unlist(strsplit(fn[ifn],"_"))[2])
	colnames(count_file) <- cn
	ty_data[[ifn]] <- count_file
	rm(count_file)
}




nrg_all <- list()

nrg_all[[1]] <- CreateSeuratObject(counts = ty_data[[1]], project = "chow1", min.cells = 3, min.features = 200)
nrg_all[[1]]$geno 	<- "chow"
nrg_all[[1]]$sn 	<- "chow1"
# nrg_all[[1]]$batch 	<- "chow_nash"
nrg_all[[1]][["percent.mt"]] <- PercentageFeatureSet(nrg_all[[1]], pattern = "^mt-")
# nrg_all[[1]] <- subset(zc_all[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 8750 & percent.mt < 20)
nrg_all[[1]] <- NormalizeData(nrg_all[[1]], normalization.method = "LogNormalize", scale.factor = 10000)
nrg_all[[1]] <- FindVariableFeatures(nrg_all[[1]], selection.method = "vst", nfeatures = 2000)


nrg_all[[2]] <- CreateSeuratObject(counts = ty_data[[2]], project = "chow2", min.cells = 3, min.features = 200)
nrg_all[[2]]$geno 	<- "chow"
nrg_all[[2]]$sn 	<- "chow2"
# nrg_all[[2]]$batch 	<- "chow_nash"
nrg_all[[2]][["percent.mt"]] <- PercentageFeatureSet(nrg_all[[2]], pattern = "^mt-")
nrg_all[[2]] <- NormalizeData(nrg_all[[2]], normalization.method = "LogNormalize", scale.factor = 10000)
nrg_all[[2]] <- FindVariableFeatures(nrg_all[[2]], selection.method = "vst", nfeatures = 2000)


nrg_all[[3]] <- CreateSeuratObject(counts = ty_data[[3]], project = "chow3", min.cells = 3, min.features = 200)
nrg_all[[3]]$geno 	<- "chow"
nrg_all[[3]]$sn 	<- "chow3"
# nrg_all[[3]]$batch 	<- "chow_nash"
nrg_all[[3]][["percent.mt"]] <- PercentageFeatureSet(nrg_all[[3]], pattern = "^mt-")
nrg_all[[3]] <- NormalizeData(nrg_all[[3]], normalization.method = "LogNormalize", scale.factor = 10000)
nrg_all[[3]] <- FindVariableFeatures(nrg_all[[3]], selection.method = "vst", nfeatures = 2000)





nrg_all[[4]] <- CreateSeuratObject(counts = ty_data[[4]], project = "nash1", min.cells = 3, min.features = 200)
nrg_all[[4]]$geno 	<- "nash"
nrg_all[[4]]$sn 	<- "nash1"
# nrg_all[[4]]$batch 	<- "chow_nash"
nrg_all[[4]][["percent.mt"]] <- PercentageFeatureSet(nrg_all[[4]], pattern = "^mt-")
nrg_all[[4]] <- NormalizeData(nrg_all[[4]], normalization.method = "LogNormalize", scale.factor = 10000)
nrg_all[[4]] <- FindVariableFeatures(nrg_all[[4]], selection.method = "vst", nfeatures = 2000)





nrg_all[[5]] <- CreateSeuratObject(counts = ty_data[[5]], project = "nash2", min.cells = 3, min.features = 200)
nrg_all[[5]]$geno 	<- "nash"
nrg_all[[5]]$sn 	<- "nash2"
# nrg_all[[5]]$batch 	<- "chow_nash"
nrg_all[[5]][["percent.mt"]] <- PercentageFeatureSet(nrg_all[[5]], pattern = "^mt-")
nrg_all[[5]] <- NormalizeData(nrg_all[[5]], normalization.method = "LogNormalize", scale.factor = 10000)
nrg_all[[5]] <- FindVariableFeatures(nrg_all[[5]], selection.method = "vst", nfeatures = 2000)



nrg_all[[6]] <- CreateSeuratObject(counts = ty_data[[6]], project = "nash3", min.cells = 3, min.features = 200)
nrg_all[[6]]$geno 	<- "nash"
nrg_all[[6]]$sn 	<- "nash3"
# nrg_all[[6]]$batch 	<- "chow_nash"
nrg_all[[6]][["percent.mt"]] <- PercentageFeatureSet(nrg_all[[6]], pattern = "^mt-")
nrg_all[[6]] <- NormalizeData(nrg_all[[6]], normalization.method = "LogNormalize", scale.factor = 10000)
nrg_all[[6]] <- FindVariableFeatures(nrg_all[[6]], selection.method = "vst", nfeatures = 2000)



# 
nrg.anchors <- FindIntegrationAnchors(object.list = nrg_all, dims = 1:30)
nrg.combined <- IntegrateData(anchorset = nrg.anchors , dims = 1:30)

DefaultAssay(object = nrg.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
nrg.combined <- ScaleData(object = nrg.combined, verbose = FALSE)
nrg.combined <- RunPCA(object = nrg.combined, npcs = 100, verbose = FALSE)

##-----------
## ALL
ep_search <- function(x){
	df 			<- cbind(1:length(x), x)
	line 		<- df[c(1, nrow(df)),]
	proj 		<- princurve::project_to_curve(df, line)
	optpoint 	<- which.max(proj$dist_ind)-1
	return(optpoint)
}


opt_dim 	<- ep_search(nrg.combined@reductions$pca@stdev)
nrg.combined <- RunUMAP(object = nrg.combined, reduction = "pca", dims = 1:opt_dim)
nrg.combined <- FindNeighbors(object = nrg.combined, reduction = "pca", dims = 1:opt_dim)
# nrg.combined <- FindClusters(nrg.combined, resolution = 0.5)
nrg.combined <- FindClusters(nrg.combined, resolution = 0.5)


nrg.markers <- FindAllMarkers(object = nrg.combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

saveRDS(nrg.combined, file = paste0(workdir,"/output/nrg_chownash_integrated_anchor30_res05.rds"))
saveRDS(nrg.markers, file = paste0(workdir,"/output/nrg_chownash_markers_integrated_anchor30_res05.rds"))


write.csv(nrg.markers,file=paste0(workdir,"/result/nrg_chownash_markers_integrated_anchor30_res05.csv"),row.names=T)



# FeaturePlot(pbmc3k.final, features = features)

p1 <- DimPlot(object = nrg.combined, reduction = "umap", group.by = "geno")
p2 <- DimPlot(object = nrg.combined, reduction = "umap", label = TRUE)
p3 <- DimPlot(object = nrg.combined, reduction = "umap", split.by="sn",ncol=3)


png(paste0(workdir,"/figure/nrg_chownash_integrated_anchor30_res05_umap.png"),width=960,height=480)
p1 + p2
dev.off()


png(paste0(workdir,"/figure/nrg_chownash_integrated_anchor30_res05_umap_bySample.png"),width=960,height=640)
DimPlot(nrg.combined, reduction = "umap", split.by = "sn",ncol=3)+ NoLegend()
dev.off()










