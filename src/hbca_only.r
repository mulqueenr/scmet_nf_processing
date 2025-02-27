singularity shell \
--bind /volumes/seq/projects/metACT \
~/singularity/amethyst.sif

library(amethyst)
library(rhdf5)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(plyr)
library(future)
library(furrr)
library(purrr)
library(cowplot)
library(pheatmap)
library(plyr)
library(parallel)
dcis_cnv<-"/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/sc_bams/all_cells.scCNA.tsv"
in_dir="/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/transfer_dat/"
in_dir2="/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/transfer_dat/"
in_dir3="/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/samples/methylation_coverage/amethyst"

setwd(in_dir3)
obj<-readRDS(file="hbca_dcis.met.Rds")


cluster_by_windows<-function(obj,window_name,stepsize.=NULL,bed.=NULL,metric.="score",threads.=200,neighbors.=50){
  print(paste("Making window summaries for ",window_name))
  obj@genomeMatrices[[window_name]] <- makeWindows(obj,
                                                      stepsize = stepsize., 
                                                      type = "CG", 
                                                      metric = metric., 
                                                      bed = bed.,
                                                      threads = threads., 
                                                      index = "chr_cg", 
                                                      nmin = 2) 
  print(paste("Estimating dimensions..."))                                           
  #filter windows by cell coverage
  obj@genomeMatrices[[window_name]] <- obj@genomeMatrices[[window_name]][rowSums(!is.na(obj@genomeMatrices[[window_name]])) >= 45, ]
  est_dim<-dimEstimate(obj, genomeMatrices = c(window_name), dims = c(10), threshold = 0.95)
  print(est_dim)
  set.seed(111)
  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = est_dim, replaceNA = c(0))
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) # Optional; helps reduce coverage 
  print("Clustering on coverage regressed reduction...")

  obj <- runCluster(obj, k_phenograph = neighbors., reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
  obj <- runUmap(obj, neighbors = neighbors., dist = 0.05, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
    obj <- runTsne(obj, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print("Plotting...")

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(window_name,"_umap.pdf"))     
  return(obj)                                             
}

obj_hbca<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$sample %in% c("HBCA-19T","hbca-83l","hbca-16r","HBCA-17T"),]))
obj_hbca<-cluster_by_windows(obj_hbca,window_name="cg_100k_score",stepsize.=100000,threads.=200)

obj_hbca<-readRDS(file="hbca.met.Rds")
window_name="cg_100k_score"
obj_hbca<- runCluster(obj_hbca, k_phenograph = 500, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets

p1 <- dimFeature(obj_hbca, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
p2 <- dimFeature(obj_hbca, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
p3 <- dimFeature(obj_hbca, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
p4 <- dimFeature(obj_hbca, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")

plt<-plot_grid(p1, p2,p3, p4,ncol=2)
ggsave(plt,file=paste0(window_name,"_umap.pdf"))  