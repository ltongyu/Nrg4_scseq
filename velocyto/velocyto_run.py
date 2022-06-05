import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
# %load_ext rpy2.ipython

adata_chow1 = anndata.read_loom("data/chow1/possorted_genome_bam_YAY4T.loom")
adata_chow2 = anndata.read_loom("data/chow2/possorted_genome_bam_UPQSQ.loom")
adata_chow3 = anndata.read_loom("data/chow3/possorted_genome_bam_1EFH0.loom")
adata_nash1 = anndata.read_loom("data/nash1/possorted_genome_bam_56RJO.loom")
adata_nash2 = anndata.read_loom("data/nash2/possorted_genome_bam_3BU2F.loom")
adata_nash3 = anndata.read_loom("data/nash3/possorted_genome_bam_D1A1L.loom")

sample_obs = pd.read_csv("/data/velo/chow_nash_mac_sub_cellID_obs.csv")
umap_cord = pd.read_csv("/data/velo/chow_nash_mac_sub_cell_embeddings.csv")
cell_clusters = pd.read_csv("/data/velo/chow_nash_mac_sub_clusters_obs.csv")


## chow1
##========
idx_chow1=[]
for x in adata_chow1.obs.index:
	idx_chow1.append(x.split(':')[1])

adata_chow1_v2 = adata_chow1.copy()
adata_chow1_v2.obs.index = idx_chow1
idx2_chow1 = [str1 for str1 in adata_chow1_v2.obs.index + "_chow1"] ## paste
adata_chow1_v2.obs.index = idx2_chow1
adata_filter_chow1 = adata_chow1_v2[np.isin(adata_chow1_v2.obs.index, sample_obs_chow1)]


## chow2
##========

idx_chow2=[]
for x in adata_chow2.obs.index:
	idx_chow2.append(x.split(':')[1])
adata_chow2_v2 = adata_chow2.copy()
adata_chow2_v2.obs.index = idx_chow2
idx2_chow2 = [str1 for str1 in adata_chow2_v2.obs.index + "_chow2"] ## paste
adata_chow2_v2.obs.index = idx2_chow2
adata_filter_chow2 = adata_chow2_v2[np.isin(adata_chow2_v2.obs.index, sample_obs_chow2)]



## chow3
##========

idx_chow3=[]
for x in adata_chow3.obs.index:
	idx_chow3.append(x.split(':')[1])

adata_chow3_v2 = adata_chow3.copy()
adata_chow3_v2.obs.index = idx_chow3
idx2_chow3 = [str1 for str1 in adata_chow3_v2.obs.index + "_chow3"] ## paste
adata_chow3_v2.obs.index = idx2_chow3
adata_filter_chow3 = adata_chow3_v2[np.isin(adata_chow3_v2.obs.index, sample_obs_chow3)]



## nash1
##========
idx_nash1=[]
for x in adata_nash1.obs.index:
	idx_nash1.append(x.split(':')[1])

adata_nash1_v2 = adata_nash1.copy()
adata_nash1_v2.obs.index = idx_nash1
idx2_nash1 = [str1 for str1 in adata_nash1_v2.obs.index + "_nash1"] ## paste
adata_nash1_v2.obs.index = idx2_nash1
adata_filter_nash1 = adata_nash1_v2[np.isin(adata_nash1_v2.obs.index, sample_obs_nash1)]


## nash2
##========
idx_nash2=[]
for x in adata_nash2.obs.index:
	idx_nash2.append(x.split(':')[1])

adata_nash2_v2 = adata_nash2.copy()
adata_nash2_v2.obs.index = idx_nash2
idx2_nash2 = [str1 for str1 in adata_nash2_v2.obs.index + "_nash2"] ## paste

adata_nash2_v2.obs.index = idx2_nash2
adata_filter_nash2 = adata_nash2_v2[np.isin(adata_nash2_v2.obs.index, sample_obs_nash2)]



## nash3
##========

idx_nash3=[]
for x in adata_nash3.obs.index:
	idx_nash3.append(x.split(':')[1])

adata_nash3_v2 = adata_nash3.copy()
adata_nash3_v2.obs.index = idx_nash3
idx2_nash3 = [str1 for str1 in adata_nash3_v2.obs.index + "_nash3"] ## paste

adata_nash3_v2.obs.index = idx2_nash3
adata_filter_nash3 = adata_nash3_v2[np.isin(adata_nash3_v2.obs.index, sample_obs_nash3)]


## necessary for concatenate
##================================
adata_filter_chow1.var_names_make_unique()
adata_filter_chow2.var_names_make_unique()
adata_filter_chow3.var_names_make_unique()

adata_filter_nash1.var_names_make_unique()
adata_filter_nash2.var_names_make_unique()
adata_filter_nash3.var_names_make_unique()


## Set index_unique to None to keep the current index
adata_allsample = adata_filter_chow1.concatenate(
adata_filter_chow2,
adata_filter_chow3,
adata_filter_nash1,
adata_filter_nash2,
adata_filter_nash3,index_unique=None)

## rename the column names
adata_allsample_index = pd.DataFrame(adata_allsample.obs.index)
adata_allsample_index = adata_allsample_index.rename(columns = {0:'Cell ID'})



## adding the umap information 
umap_cord = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered = adata_allsample_index.merge(umap_cord, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
adata_allsample.obsm['X_umap'] = umap_ordered.values



## add the cluster information 
clust_info = cell_clusters.rename(columns = {'Unnamed: 0':'Cell ID'})
clust_info_ordered = adata_allsample_index.merge(clust_info, on = "Cell ID")

adata_allsample.uns['Cluster_colors'] = clust_info_ordered.iloc[:,2:3]
adata_allsample.uns['Seurat_Cluster'] = clust_info_ordered.iloc[:,1:2]


scv.pp.filter_and_normalize(adata_allsample)
scv.pp.moments(adata_allsample)
scv.tl.velocity(adata_allsample, mode = "stochastic")
scv.tl.velocity_graph(adata_allsample) #this step is time consuming 



fig = scv.pl.velocity_embedding_stream(adata_allsample, basis='umap',color=adata_allsample.uns['Cluster_colors'].iloc[:,0],colorbar=False,show=False,size=10)
ff 	= fig.get_figure()

ff.savefig('figure/adata_chow_nash_macrophage_sub.png',dpi=300)


