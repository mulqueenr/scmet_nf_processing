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

#forgot to add optparse to SIF, so just using ordered list of inputs
args <- commandArgs(trailingOnly = TRUE)
input_dir=args[1] #"Dir of single-cell bam files"
output_prefix=args[2] #Prefix of output
metadata=args[3] #"Input of metadata from METHYLATION_CALL csv output."
task_cpus=args[4] #"Integer number of cpus"

cpu_count=task_cpus
prefix=output_prefix
metadata_in=metadata

obj <- createObject()

#metadata MUST have a column called mcg_pct for score calculation
#metadata MUST have a column called cov to regress coverage mias
metadat<-read.table(metadata,sep=",")
colnames(metadat)<-c("cellid","mcg_cov","cov","mcg_pct")
row.names(metadat)<-metadat$cellid
obj@metadata<-metadat

head(obj@metadata)
plt<-ggplot(obj@metadata, aes(x=prefix, y = cov)) +geom_violin() + geom_jitter()
ggsave(file=paste0(prefix,"_cov_plot.pdf"),plt)

h5paths<-list.files(path=input_dir,pattern="*.h5.gz",full.names=TRUE)
cellid<-basename(h5paths) #make sure this matches with metadata
obj@h5paths <- data.frame(row.names = c(cellid), paths = h5paths)

# index files
obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = cpu_count) 


cluster_by_windows<-function(obj,window_name,stepsize.=NULL,bed.=NULL,metric.="score",threads.=100,neighbors.=50){
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

  obj <- runCluster(obj, k_phenograph = 175, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
  obj <- runUmap(obj, neighbors = neighbors., dist = 0.05, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
  print("Plotting...")

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sampleName, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(window_name,"_umap.pdf"))     
  return(obj)                                             
}

obj<-cluster_by_windows(obj,window_name="cg_100k_score",stepsize.=100000,threads.=cpu_count)


#Annotate Data
#obj@ref <- makeRef(ref="HG38",gtf="/volumes/USR2/Ryan/ref/gencode.v43.annotation.gtf") #ref stored within amethyst.sif container
gtf <- rtracklayer::readGFF("/container_ref/gencode.v43.annotation.gtf.gz")
for (i in c("gene_name", "exon_number")) {gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))}
gtf <- dplyr::mutate(gtf, location = paste0(seqid, "_", start, "_", end))
obj@ref<-gtf
protein_coding <- unique(obj@ref |> dplyr::filter(type == "gene" & gene_type == "protein_coding" & seqid != "chrM") |> dplyr::pull(gene_name))

obj@genomeMatrices[["cg_promoter"]] <- makeWindows(obj, 
                                                     genes = protein_coding,
                                                     promoter = TRUE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = cpu_count, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

obj@genomeMatrices[["cg_genebody"]] <- makeWindows(obj, 
                                                     genes = protein_coding,
                                                     promoter = FALSE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = cpu_count, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

saveRDS(obj,file=paste0(prefix,".amethyst.rds"))

