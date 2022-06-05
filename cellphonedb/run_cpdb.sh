
## get the cpdb activated
##+=======================
bash
source activate cpdb

## for cd_mac
##+=======================
cellphonedb method statistical_analysis \
CD8_MAC_meta_hgnc.tsv \
cd_mac \
--counts-data hgnc_symbol \
--database v2.0.0 \
--iterations 10000 \
--output-path output/cpdb \
--means-result-name means_stat_cd_mac_10k \
--significant-means-result-name significant_means_stat_cd_mac_10k \
--deconvoluted-result-name deconvoluted_stat_cd_mac_10k \
--pvalues-result-name pvals_stat_cd_mac_10k \
--threads 10 
