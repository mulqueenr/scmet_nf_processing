```bash
singularity shell \
--bind /volumes/seq/projects/metACT \
~/singularity/amethyst.sif
```

```R
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

#read in all sample/csv/sample.passingCellsMapMethylStats.csv files into data frame
#make a dataframe of all h5 files also <sample>\t<h5location>
dcis_cnv<-"/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/sc_bams/all_cells.scCNA.tsv"
in_dir="/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/transfer_dat/"
in_dir2="/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/transfer_dat/"
setwd(in_dir2)
#obj<-readRDS("hbca_dcis.met.Rds")
obj <- createObject()

#read in metadata
samp_dir=list.files(path=paste0(in_dir,"/report/"))
samp_dir2=list.files(path=paste0(in_dir2,"/report/"))
met1<-lapply(samp_dir,function(i){
  samp_met<-read.csv(paste0(in_dir,"/report/",i,"/csv/",i,".passingCellsMapMethylStats.csv"))
  row.names(samp_met)<-samp_met$BC
  return(samp_met)
})
metadat1<-do.call("rbind",met1)
metadat1$run<-1
metadat1<-metadat1[metadat1$CG_Cov>10000,]
metadat1<-metadat1[metadat1$sampleName %in% c("hbca-83l","hbca-16r","DCIS-41T","DCIS-66T"),]

met2<-lapply(samp_dir2,function(i){
  samp_met<-read.csv(paste0(in_dir2,"/report/",i,"/csv/",i,".passingCellsMapMethylStats.csv"))
  row.names(samp_met)<-samp_met$BC
  return(samp_met)
})
metadat2<-do.call("rbind",met2)
metadat2$run<-2
metadat2<-metadat2[metadat2$CG_Cov>10000,]
metadat2<-metadat2[metadat2$sampleName %in% c("hbca-83l" ,"hbca-16r" ,"DCIS-41T", "DCIS-66T"),]

metadat<-rbind(metadat1,metadat2)

#Add CNV metadata for DCIS
cnv_dat<-read.table(dcis_cnv,sep="\t",header=T)
cnv_dat$idx<-gsub("_","+",cnv_dat$idx)
metadat$cnv_clone<-"diploid"
metadat[match(cnv_dat$idx,row.names(metadat)),]$cnv_clone<-cnv_dat$superclones

#metadata MUST have a column called mcg_pct for score calculation
#metadata MUST have a column called cov to regress coverage mias
metadat$mcg_pct<-metadat$CG_mC_Pct
metadat$cov<-metadat$CG_Cov
obj@metadata<-metadat

head(obj@metadata)
obj@metadata<-obj@metadata[obj@metadata$CG_Cov>10000,]
plt<-ggplot(obj@metadata, aes(x=sampleName, y = CG_Cov)) +geom_violin() + geom_jitter()
ggsave(file="cov_plot.pdf",plt)

#filter cells
#obj@metadata <- obj@metadata |> dplyr::filter(cov > 100000 & cov < 40000000)
h5paths<-paste0(in_dir,"/cg_sort_cov/h5_files/",metadat1$sampleName,".",metadat1$tgmt_well,".h5")
h5paths2<-paste0(in_dir2,"/h5_files/",metadat2$sampleName,".",metadat2$tgmt_well,".h5")
obj@h5paths <- data.frame(row.names = c(rownames(metadat)), paths = c(h5paths,h5paths2))
head(obj@h5paths)

# index files
obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = 100) 

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

obj<-cluster_by_windows(obj,window_name="cg_100k_score",stepsize.=100000,threads.=200)


#rerun clustering (overcluster and then merge them visually)
window_name="cg_100k_score"
obj <- runCluster(obj, k_phenograph = 25, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
obj <- runUmap(obj, neighbors = 50, dist = 0.0001, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
print("Plotting...")
p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
p2 <- dimFeature(obj, colorBy = sampleName, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
p5 <- dimFeature(obj, colorBy = cnv_clone, reduction = "umap") + ggtitle(paste(window_name,"CNV Clones"))
plt<-plot_grid(p1, p2,p3, p4, p5, ncol=2)
ggsave(plt,file=paste0(window_name,"_umap.pdf"),height=10)    

table(obj@metadata$cluster_id)
#  1  10  11  12  13  14  15  16  17  18  19   2  20   3   4   5   6   7   8   9 
# 87  92 100 129  79 109  42  91 112  73  43 171  65 129 159 112 199  29 100 265 
#plot by feature

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
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

obj@genomeMatrices[["cg_genebody"]] <- makeWindows(obj, 
                                                     genes = protein_coding,
                                                     promoter = FALSE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 


cluster_by_windows<-function(obj,window_name=cg_genebody,est_dim=5,threads.=100,neighbors.=50){
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

#############Gene body methylation#############
cluster_by_windows(obj,est_dim=3, window_name="cg_genebody")
cluster_by_windows(obj,est_dim=3,window_name="cg_promoter")

#find cluster markers on genebodies
cluster_genebody_markers <- FindClusterMarkers(obj, 
                                               matrix = "cg_genebody", 
                                               genes = protein_coding, group_by="cluster_id",
                                               threads = 400)


cluster_genebody_markers <- cluster_genebody_markers |> dplyr::filter(p.val < 0.1)
write.table(cluster_genebody_markers,file="cluster_genebody_markers.tsv",sep="\t",col.names=T)
head(cluster_genebody_markers)

posteriori_markers<-cluster_genebody_markers |> 
dplyr::filter(p.val < 0.05)  |> 
dplyr::filter(is.finite(logFC)) |> 
dplyr::filter(!duplicated(gene)) |> 
dplyr::filter(direction=="hypomethylated") |> 
group_by(cluster_id) |> slice_min(n = 10, order_by = p.adj) 


paste_markers<-function(i){
  print(posteriori_markers[posteriori_markers$cluster_id==i,])
  print(i)
  paste(unlist(posteriori_markers[posteriori_markers$cluster_id==i,]$gene),collapse=",")}

lapply("19",paste_markers)


plt_list<-lapply(unique(posteriori_markers$cluster_id), function(i) {
  plt<-dotM(obj, genes = posteriori_markers[posteriori_markers$cluster_id == i,]$gene, groupBy = "cluster_id", matrix = "cg_genebody") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "black")) + scale_size(range = c(1, 15)) + guides(fill="none")
  return(plt)
})
plt_out<-plot_grid(plotlist=plt_list,nrow=1)
ggsave(plt_out,file="denovomarkers_genebody_dotplot_clusterid.pdf",height=10,width=length(features)*8,limitsize=F)


###########################################################################################

###### DMR Calculation ####
obj_dcis<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$sample %in% c("DCIS-41T"),]))


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
  return(obj)                                             
}
   

obj_dcis<-cluster_by_windows(obj_dcis,window_name="cg_100k_score",stepsize.=100000,threads.=200,neighbors.=400)
window_name="cg_100k_score"
p1 <- dimFeature(obj_dcis, colorBy = cnv_clone, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
ggsave(p1,file=paste0(window_name,"_umap.pdf"))  

cluster1kbwindows <- calcSmoothedWindows(obj_dcis, 
                                         type = "CG", 
                                         threads = 300,
                                         step = 1000,
                                         smooth = 3,
                                         index = "chr_cg",
                                         groupBy = "cnv_clone",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
obj_dcis@genomeMatrices[["cg_cluster_tracks"]] <- cluster1kbwindows[["pct_matrix"]]
pal=c("#E5E6E4","#CFD2CD","#A6A2A2","#847577","#6E44FF")
dmrs<-testDMR(cluster1kbwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 0, nminGroup = 0) # or use cluster1kbwindows[["sum_matrix"]] and rename
dmrs2<-filterDMR(dmrs, method = "bonferroni", filter = FALSE) #add additional columns direction column
celltype_test<-setNames(nm=unique(dmrs2$test), gsub("_c","",colnames(cluster1kbwindows[["sum_matrix"]])[grepl(pattern="_c",colnames(cluster1kbwindows[["sum_matrix"]]))]))
dmrs2$celltype<-celltype_test[dmrs2$test]

collapsed_dmrs <- collapseDMR(obj_dcis, dmrs2, maxDist = 4000, minLength = 1000, reduce = T, annotate = T) 
collapsed_dmrs$celltype<-celltype_test[collapsed_dmrs$test]
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(celltype, direction) |> dplyr::summarise(n = n()), 
       aes(y = celltype, x = n, fill = celltype)) + geom_col() + 
  facet_grid(vars(direction), scales = "free_y") + scale_fill_manual(values = makePalette(option = 7, n = 13) ) + theme_classic()
ggsave(plt,file="met_per_dmr_cnv_cluster.pdf")

top_dmrs <- collapsed_dmrs |> 
  dplyr::group_by(celltype, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  group_by(celltype, direction) |> slice_min(n = 20, order_by = total_rank) |>
  dplyr::mutate(location = paste0(chr, "_", (dmr_start - 2000), "_", (dmr_end + 2000))) |> dplyr::arrange(direction)
write.table(top_dmrs,file="cluster_dmr_fine.tsv")

plt<-heatMap(obj, matrix = "cg_celltype_tracks", regions = top_dmrs$location, nrow = 10, legend = FALSE, width = 1000, arrowOverhang = 2000)
ggsave(plt,file="dmrs_heatmap_celltype.pdf",width=50,height=50,limitsize=FALSE)


###########################################################################################

#Way to cell type is: overcluster to 1kb windows, put all clusters into IGV, color code by grouping, use both scRNA and snRNA markers to look for enrichment. CGI annotation, since CGI are univerally low methylation.

######### Set Initial Celltypes #############
celltypes<-setNames(nm=c(
  3,8,
  5,
  1,
  10,
  11,
  4,7,
  19,
  18,6,17,
  2,9,12,
  16,14,
  13,
  20,
  15),
c("fibro","fibro",
"tcell",
"unknown1", #endo
"immune", #macro/tcell
"unknown2", 
"basal","basal", #basal
"adipo",
"lumsec","lumsec","lumsec", #stromal
"lumhr","lumhr","lumhr",
"dcis_lumhr1","dcis_lumhr1", #these look weird, maybe stromal associated with DCIS?
"dcis_lumhr2",
"dcis_lumhr3",
"dcis_lumhr4"))

obj@metadata$celltype<-celltypes[obj@metadata$cluster_id]


###### DMR on cell types ####
celltype1kbwindows <- calcSmoothedWindows(obj, 
                                         type = "CG", 
                                         threads = 300,
                                         step = 1000,
                                         smooth = 3,
                                         index = "chr_cg",
                                         groupBy = "celltype",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
obj@genomeMatrices[["cg_celltype_tracks"]] <- celltype1kbwindows[["pct_matrix"]]

saveRDS(obj,file="hbca_dcis.met.Rds")
obj<-readRDS(file="hbca_dcis.met.Rds")

pal=c("#E5E6E4","#CFD2CD","#A6A2A2","#847577","#6E44FF")
dmrs<-testDMR(celltype1kbwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 0, nminGroup = 0) # or use cluster1kbwindows[["sum_matrix"]] and rename
dmrs2<-filterDMR(dmrs, method = "bonferroni", filter = FALSE) #add additional columns direction column
celltype_test<-setNames(nm=unique(dmrs2$test), gsub("_c","",colnames(celltype1kbwindows[["sum_matrix"]])[grepl(pattern="_c",colnames(celltype1kbwindows[["sum_matrix"]]))]))
dmrs2$celltype<-celltype_test[dmrs2$test]
save(dmrs2,file="hbca_dcis.celltype.dmr.Rds")

collapsed_dmrs <- collapseDMR(obj, dmrs2, maxDist = 4000, minLength = 2000, reduce = T, annotate = T) 
collapsed_dmrs$celltype<-celltype_test[collapsed_dmrs$test]
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(celltype, direction) |> dplyr::summarise(n = n()), 
       aes(y = celltype, x = n, fill = celltype)) + geom_col() + 
  facet_grid(vars(direction), scales = "free_y") + scale_fill_manual(values = makePalette(option = 7, n = 13) ) + theme_classic()
ggsave(plt,file="met_per_dmr_celltype.pdf")

top_dmrs <- collapsed_dmrs |> 
  dplyr::group_by(celltype, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  group_by(celltype, direction) |> slice_min(n = 5, order_by = total_rank) |>
  dplyr::mutate(location = paste0(chr, "_", (dmr_start - 2000), "_", (dmr_end + 2000))) |> dplyr::arrange(direction)

plt<-heatMap(obj, matrix = "cg_celltype_tracks", regions = top_dmrs$location, nrow = 10, legend = FALSE, width = 1000, arrowOverhang = 2000)
ggsave(plt,file="dmrs_heatmap_celltype.pdf",width=50,height=50,limitsize=FALSE)




#find celltype markers on genebodies
celltype_genebody_markers <- FindClusterMarkers(obj, 
                                               matrix = "cg_genebody", 
                                               genes = protein_coding, group_by="celltype",
                                               threads = 200)


celltype_genebody_markers <- celltype_genebody_markers |> dplyr::filter(p.val < 0.1)
write.table(celltype_genebody_markers,file="celltype_genebody_markers.tsv",sep="\t",col.names=T)
head(celltype_genebody_markers)

posteriori_markers<-celltype_genebody_markers |> 
dplyr::filter(p.adj < 0.05)  |> 
dplyr::filter(is.finite(logFC)) |> 
dplyr::filter(!duplicated(gene)) |> 
dplyr::filter(direction=="hypomethylated") |> 
group_by(cluster_id) |> slice_min(n = 10, order_by = p.adj) 

plt_list<-lapply(unique(posteriori_markers$cluster_id), function(i) {
  plt<-dotM(obj, genes = posteriori_markers[posteriori_markers$cluster_id == i,]$gene, groupBy = "celltype", matrix = "cg_genebody") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "black")) + scale_size(range = c(1, 15)) + guides(fill="none")
  return(plt)
})
plt_out<-plot_grid(plotlist=plt_list,nrow=1)
ggsave(plt_out,file="denovomarkers_genebody_dotplot_celltype.pdf",height=10,width=length(features)*8,limitsize=F)


#######################################

###### ChromVAR on DMR sites ####
library(BiocParallel)
register(MulticoreParam(5)) # Use 50 cores
library(chromVAR)
library(motifmatchr)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38) #locally installed
library(TFBSTools)
library(SummarizedExperiment)
library(dendextend)#local
library(reshape2)


dmr_counts<-dmrs2 |> 
  filter(padj < 0.05) |> 
  filter(direction=="hypo")

peaks <- GRanges(dmrs2)
mcols(peaks)<-NULL
peaks<-unique(peaks)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = FALSE))

motif_ix <- matchMotifs(pfm, peaks, genome = "hg38") 

hypergeo_motif<-function(cluster){
  motif_dmrs<-findOverlaps(
                          query=GRanges(dmr_counts[dmr_counts$test==cluster,]),
                          subject=motif_ix,minoverlap=min(width(motif_ix)))
  
  motif_dmr_cluster<-setNames(
    unlist(lapply(1:ncol(motif_ix@assays@data$motifMatches),
    function(i){
      sum(motif_ix@assays@data$motifMatches[unique(subjectHits(motif_dmrs)),i])
      })),
      nm=motif_ix@colData$name)

  motif_dmr_background<-setNames(
                                colSums(motif_ix@assays@data$motifMatches),
                                nm=motif_ix@colData$name)

  pval_calc<-do.call("rbind",lapply(motif_ix@colData$name, function(motif){
  q=motif_dmr_cluster[motif]
  m=motif_dmr_background[motif]
  n=nrow(motif_ix)-motif_dmr_background[motif]
  k=nrow(dmr_counts[dmr_counts$test==cluster,])
  pval<-phyper(q=q,m=m,n=n,k=k,lower.tail=FALSE)
  out_calc<-c(q,m,n,k,pval)
  return(out_calc)
  }))
out_frame=cbind(data.frame(cluster=cluster,motif=motif_ix@colData$name),pval_calc)
return(out_frame)
}

motif_enrichment<-do.call("rbind",lapply(unique(dmr_counts$test),hypergeo_motif))
colnames(motif_enrichment)<-c("cluster","motif","q","m","n","k","pval")
motif_enrichment$padj<-p.adjust(motif_enrichment$pval, method = "bonferroni", n = length(motif_enrichment$pval))
motif_enrichment$celltype<-celltype_test[motif_enrichment$cluster]

#get list of TFs that are sig enriched at least once
out_heatmap<-dcast(motif_enrichment, motif ~ celltype, value.var="padj")
row.names(out_heatmap)<-out_heatmap$motif
out_heatmap<-out_heatmap[2:ncol(out_heatmap)]
out_heatmap<-out_heatmap[rowSums(out_heatmap<0.05)>1,] #filter out rows with no sig

hclust_tfs <- hclust(dist(out_heatmap), method = "complete")
hclust_cols <- hclust(dist(t(out_heatmap)), method = "complete")
plt<-pheatmap(out_heatmap,cluster_rows=hclust_tfs,cluster_cols=TRUE,fontsize = 8)
ggsave(plt,file="dmr_tfmotif_enrichment.pdf")

###########################################################################################

###########################################################################################

#snRNA markers
hbca_snmarkers=list()
hbca_snmarkers[["lumhr"]]=c("ANKRD30A","AFF3","ERBB4","TTC6","MYBPC1","NEK10","THSD4")
hbca_snmarkers[["lumsec"]]=c("COBL","GABRP","ELF5","CCL28","KRT15","KIT") #"AC011247.1",
hbca_snmarkers[["basal"]]=c("CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2") #"AC044810.2",
hbca_snmarkers[["fibro"]]=c("LAMA2","DCLK1","NEGR1","LINC02511","ANK2","KAZN","SLIT2")
hbca_snmarkers[["lymphatic"]]=c("PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1") #"AL357507.1",
hbca_snmarkers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
hbca_snmarkers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1") #"AC012409.2"
hbca_snmarkers[["myeloid"]]=c("F13A1","MRC1","RBPJ","TBXAS1","FRMD4B","CD163","RAB31")
hbca_snmarkers[["tcells"]]=c("SKAP1","ARHGAP15","PTPRC","THEMIS","IKZF1","PARP8","CD247")
hbca_snmarkers[["mast"]]=c("NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC")
hbca_snmarkers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
features<-llply(hbca_snmarkers, unlist)


#post markers for easy entry into cellxgene
tmp<-as.data.frame(cluster_genebody_markers %>% 
dplyr::filter(p.adj < 0.05)  %>%  
dplyr::filter(is.finite(logFC)) %>% 
dplyr::filter(!duplicated(gene)) %>% 
dplyr::filter(direction=="hypomethylated"))

tmp<- tmp %>% 
arrange(p.adj, logFC) %>%
group_by(cluster_id) %>% 
top_n(50, wt=p.adj) %>% 
select(cluster_id, gene)

paste_markers<-function(i){
  print(tmp[tmp$cluster_id==i,])
  print(i)
  paste(unlist(tmp[tmp$cluster_id==i,]$gene),collapse=",")}

lapply("7",paste_markers)

obj@metadata$celltype<-first_pass_clusters[obj@metadata$cluster_id]

plt <- dimFeature(obj, colorBy = celltype, reduction = "umap") + ggtitle(paste(window_name,"Coarse Clusters"))
ggsave(plt,file=paste0(window_name,"_celltype_umap.pdf"))  


clusterfine1kbwindows <- calcSmoothedWindows(obj, 
                                         type = "CG", 
                                         threads = 300,
                                         step = 1000,
                                         smooth = 3,
                                         index = "chr_cg",
                                         groupBy = "cluster_fine",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
obj@genomeMatrices[["cg_clusterfine_tracks"]] <- clusterfine1kbwindows[["pct_matrix"]]


saveRDS(obj,file="hbca_dcis.met.Rds")

pal=c("#E5E6E4","#CFD2CD","#A6A2A2","#847577","#6E44FF")
dmrs<-testDMR(clusterfine1kbwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 0, nminGroup = 0) # or use cluster1kbwindows[["sum_matrix"]] and rename
dmrs2<-filterDMR(dmrs, method = "bonferroni", filter = FALSE) #add additional columns direction column
celltype_test<-setNames(nm=unique(dmrs2$test), gsub("_c","",colnames(clusterfine1kbwindows[["sum_matrix"]])[grepl(pattern="_c",colnames(clusterfine1kbwindows[["sum_matrix"]]))]))
dmrs2$celltype<-celltype_test[dmrs2$test]
save(dmrs2,file="hbca_dcis.celltype_fine.dmr.Rds")

collapsed_dmrs <- collapseDMR(obj, dmrs2, maxDist = 4000, minLength = 1000, reduce = T, annotate = T) 
collapsed_dmrs$celltype<-celltype_test[collapsed_dmrs$test]
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(celltype, direction) |> dplyr::summarise(n = n()), 
       aes(y = celltype, x = n, fill = celltype)) + geom_col() + 
  facet_grid(vars(direction), scales = "free_y") + scale_fill_manual(values = makePalette(option = 7, n = 13) ) + theme_classic()
ggsave(plt,file="met_per_dmr_celltype.pdf")

top_dmrs <- collapsed_dmrs |> 
  dplyr::group_by(celltype, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  group_by(celltype, direction) |> slice_min(n = 20, order_by = total_rank) |>
  dplyr::mutate(location = paste0(chr, "_", (dmr_start - 2000), "_", (dmr_end + 2000))) |> dplyr::arrange(direction)
write.table(top_dmrs,file="cluster_dmr_fine.tsv")

plt<-heatMap(obj, matrix = "cg_celltype_tracks", regions = top_dmrs$location, nrow = 10, legend = FALSE, width = 1000, arrowOverhang = 2000)
ggsave(plt,file="dmrs_heatmap_celltype.pdf",width=50,height=50,limitsize=FALSE)





#plot histogram of markers 
plt_list<-lapply(names(features), function(i) {
  feat=unlist(features[[i]])
  feat<-feat[feat%in% protein_coding]
  plt<-histograM(obj, genes = feat, trackOverhang=20000,
        matrix = "cg_celltype_tracks", 
        width = 1000)
  ggsave(plt,file=paste0("marker_histograms.",i,".pdf"),width=30)
})

#plot heatmap of markers 
plt_list<-lapply(names(features), function(i) {
  feat=unlist(features[[i]])
  feat<-feat[feat%in% protein_coding]
  plt<-heatMap(obj, genes = feat,
        matrix = "cg_cluster_tracks", trackOverhang=20000,
        width = 1000)
  ggsave(plt,file=paste0("marker_heatmaps.",i,".pdf"),width=5,height=10)
})


#plot dotplot of markers
plt_list<-lapply(names(features), function(i) {
  feat=unlist(features[[i]])
  plt<-dotM(obj, genes = feat, groupBy = "celltype", matrix = "cg_genebody") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "#999999")) +
  scale_size(range = c(1, 10)) + guides(fill="none") + ylab(i)
  return(plt)
})
plt_out<-plot_grid(plotlist=plt_list,ncol=1)
ggsave(plt_out,file="marker_genebody_dotplot.pdf",width=10,height=length(features)*5,limitsize=F)

saveRDS(obj,file="hbca_dcis.met.Rds")



plt_list<-lapply(names(features), function(i) {
  plt<-dotM(obj, genes = features[[i]], groupBy = "celltype", matrix = "cg_promoter") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "#999999")) + scale_size(range = c(1, 8)) +coord_flip() + ggtitle(i)
  return(plt)
})
plt_out<-plot_grid(plotlist=plt_list,ncol=1)
ggsave(plt_out,file="marker_promoter_dotplot.pdf",width=10,height=length(features)*5,limitsize=F)


plt_list<-lapply(names(features), function(i) {
  plt<-dotM(obj, genes = features[[i]], groupBy = "celltype", matrix = "cg_genebody") + 
  scale_color_gradientn(colors =  c("#FF0082", "#dbdbdb", "#cccccc", "#999999")) + scale_size(range = c(1, 8)) +coord_flip() + ggtitle(i)
  return(plt)
})
plt_out<-plot_grid(plotlist=plt_list,ncol=1)
ggsave(plt_out,file="marker_genebody_dotplot.pdf",width=10,height=length(features)*5,limitsize=F)



#post markers for easy entry into cellxgene

tmp<-as.data.frame(cluster_genebody_markers %>% 
dplyr::filter(p.adj < 0.05)  %>%  
dplyr::filter(is.finite(logFC)) %>% 
dplyr::filter(!duplicated(gene)) %>% 
dplyr::filter(direction=="hypomethylated"))

tmp<- tmp %>% 
arrange(p.adj, logFC) %>%
group_by(cluster_id) %>% 
top_n(50, wt=p.adj) %>% 
select(cluster_id, gene)

paste_markers<-function(i){
  print(tmp[tmp$cluster_id==i,])
  print(i)
  paste(unlist(tmp[tmp$cluster_id==i,]$gene),collapse=",")}

#[1] "Macro"          "DCIS_Lum_Basal" "Basal"          "DCIS_Lum"      
#[5] "Fibro"          "LumHR"          "Stromal"        "LumSec" 
lapply("Basal",paste_markers)


plt1<-sampleComp(obj, groupBy = "sampleName", colorBy = "celltype") 
plt2<-sampleComp(obj,groupBy="cnv_clone",colorBy="celltype")
plt<-plot_grid(plt1, plt2, ncol=2)
ggsave(plt,file="sample_comp.pdf")

#DMR calculation

celltype1kbwindows <- calcSmoothedWindows(obj, 
                                         type = "CG", 
                                         threads = 300,
                                         step = 1000,
                                         smooth = 3,
                                         species = "human",
                                         index = "chr_cg",
                                         groupBy = "celltype",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
obj@genomeMatrices[["cg_celltype_tracks"]] <- celltype1kbwindows[["pct_matrix"]]



for(i in unique(dmr_genes$celltype)){
  for(gene in dmr_genes[dmr_genes$celltype==i,]$gene){
    plt1 <- dimFeature(obj, colorBy = obj@metadata$celltype, pointSize = 1)  
    plt2 <- dimFeature(obj, colorBy = gene, matrix="cg_genebody", pointSize = 1) + 
            scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + 
            ggtitle(gene)
    plt3<-histograM(obj, genes = gene, trackOverhang=20000,
            matrix = "cg_celltype_tracks", 
            width = 1000,legend=FALSE) 
    plt4<-heatMap(obj, genes = gene, trackOverhang=20000,
            matrix = "cg_celltype_tracks", 
            width = 1000,legend=FALSE) 
    plt_out<-plot_grid(plotlist=list(plt1,plt2,plt3,plt4),nrow=2)

    ggsave(plt_out,file=paste0("marker_celltype_",gene,"_",i,".pdf"),width=10)
  }
}


#met markers (snRNA and scRNA from HBCA and spot checked on IGV at the cluster level)
met_markers=list()
met_markers[["lumhr"]]=c("ANKRD30A","KRT18","AFF3","PIP","MYBPC1")
met_markers[["lumsec"]]=c("COBL","GABRP","ELF5","CCL28","KRT15","KIT","LTF") #"AC011247.1",
met_markers[["basal"]]=c("CARMN","KRT5","KRT14","KRT17","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2") #"AC044810.2",
met_markers[["fibro"]]=c("COL1A2","DCN","APOD","LUM","COL1A1","LAMA2","DCLK1","ANK2")
met_markers[["vascular_endo"]]=c("MECOM","BTNL9")
met_markers[["myeloid"]]=c("HLA-DPA1","FAM78A","IL1B","PTPRC","CXCR4")
met_markers[["tcells"]]=c("SIT1","GIMAP7","IKZF1")
met_markers[["mast"]]=c("NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC")
met_markers[["adipo"]]=c("TRHDE")
features<-llply(met_markers, unlist)

for(i in names(features)){
  for(j in unlist(features[i])){
    plt3<-histograM(obj, genes = j, trackOverhang=20000,
            matrix = "cg_celltype_tracks", invert_met=TRUE,
            width = 1000,legend=FALSE) 
    print(paste("Saving plot for ",i,j))
    ggsave(plt3,file=paste0("marker_celltype_",i,"_",j,".pdf"),width=10)
  }
  }



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
library(rtracklayer) #local

dcis_cnv<-"/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/sc_bams/all_cells.scCNA.tsv"
in_dir="/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/transfer_dat/"
in_dir2="/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/transfer_dat/"
setwd(in_dir2)
obj<-readRDS(file="hbca_dcis.met.Rds")

for(i in 4:ncol(obj@genomeMatrices$cg_cluster_tracks)){
cluster=names(obj@genomeMatrices$cg_cluster_tracks)[i]
out_bw<-as.data.frame(obj@genomeMatrices$cg_cluster_tracks)
out_bw<-out_bw[c("chr","start","end",cluster)]
out_bw<-GRanges(out_bw[complete.cases(out_bw),])
names(out_bw@elementMetadata)<-"score"
out_bw<-out_bw[unique(findOverlaps(out_bw, type = "any", select = "first"))]
out_bw <- resize(out_bw, width=1000, fix='start') #resize to avoid 1base overlap
genome(out_bw)<-"hg38"
hg38_seq_info<-Seqinfo(genome="hg38")
seqlengths(out_bw)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_bw@seqnames,]$seqlengths
print(paste("Saving bigwig for...",cluster))
export(out_bw,con=paste0(cluster,"_cluster.bw"),format='bigWig')
}


colors=c(
"basal"="#FA62DB","lumsec"="#00BF7D","lumhr"="#001eff",
"dcis_lumhr2"="#BEFFFF","dcis_lumhr1"="#BEFFFF","dcis_lumhr3"="#BEFFFF","dcis_lumhr4"="#BEFFFF",
"fibro"="#FF00F5","immune"="#FF6700","tcell"="#8800FF",
"unknown1"="#FF99B6","unknown2"="#FF99B6","adipo"="#FF99B6")
#4deeea
library(tidyverse)
library(colorspace)
#set colors for each facet and use greyscaling for higher met levels
#https://stackoverflow.com/questions/33221794/separate-palettes-for-facets-in-ggplot-facet-grid
histograM <- function(obj,
                      genes = c("KIT"),
                      matrix = "cg_celltype_tracks",
                      colors = NULL,
                      trackOverhang = 5000,
                      arrowOverhang = 3000,
                      ncol = length(genes),
                      invert_met = TRUE,
                      legend = TRUE,
                      removeNA = TRUE,
                      width = 1000,
                      trackScale = 1.5,
                      colorMax = 100) {

  if (!is.null(colors)) {
    pal <- colors
  } else {
    pal <- rev(c("#FF0082", "#dbdbdb", "#cccccc", "#999999"))
  }

  if (is.null(obj@ref)) {
    stop("Please make sure a genome annotation file has been added to the obj@ref slot with makeRef().")
  }

  p <- vector("list", length(genes)) # empty plot list
  for (i in 1:length(genes)) {
    ref <- obj@ref |> dplyr::filter(gene_name == genes[i])
    aggregated <- obj@genomeMatrices[[matrix]]

    toplot <- aggregated[c((aggregated$chr == ref$seqid[ref$type == "gene"] &
                              aggregated$start > (ref$start[ref$type == "gene"] - trackOverhang) &
                              aggregated$end < (ref$end[ref$type == "gene"] + trackOverhang))), ]
    ngroups <- ncol(toplot) - 3
    trackHeight <- ngroups * trackScale
    toplot <- tidyr::pivot_longer(toplot, cols = c(4:ncol(toplot)), names_to = "group", values_to = "pct_m") |> dplyr::rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))
    if (removeNA) {
      toplot <- toplot |> dplyr::filter(!is.na(pct_m))
    }
    if (invert_met){
      toplot$pct_m<-toplot$pct_m-100
      colorMax=-100
    }
    #initialize first row
    p[[i]] <- vector("list", length(colors)) # empty plot list
    
    for (j in names(colors)){
    toplot_sub<-toplot[toplot$group==j,]
    p[[i]][[j]] <- ggplot2::ggplot() + 
      ggplot2::geom_col(data = toplot_sub, ggplot2::aes(x = middle, y = pct_m, fill = pct_m), width = width) + 
      scale_fill_gradient2(high="black",mid="grey",low=colors[j],midpoint=-60) + 
      theme_void() +
      theme(legend.position="none") + ylab(j)
    if(j==first(names(colors))){
    p[[i]][[j]] <- p[[i]][[j]] + 
      ggplot2::geom_rect(fill = "pink", data = ref |> dplyr::filter(type == "gene") |>
                          dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end + 1500)),
                          promoter_end = ifelse(strand == "+", (promoter_start+3000), (promoter_start-3000))),
                          ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = trackHeight*3, ymax = trackHeight)) +
      ggplot2::geom_rect(fill = "black", data = ref |> dplyr::filter(type == "exon"),
                         ggplot2::aes(xmin = start, xmax = end, ymin = trackHeight*3, ymax = trackHeight)) +
      ggplot2::geom_segment(data = ref, aes(x = ifelse(strand == "+", (min(start) - arrowOverhang), (max(end)) + arrowOverhang),
                          y = (trackHeight*1.5),
                          xend = ifelse(strand == "+", (max(end) + arrowOverhang), (min(start)) - arrowOverhang),
                          yend = (trackHeight*1.5)), arrow = arrow(length = unit(trackHeight/40, "cm"))) + 
                          xlab(paste(genes[i],toplot$chr[1],toplot$start[1],"-",toplot$end[nrow(toplot)])) +
                          theme_void() + theme(legend.position="none")
      p[[i]][[j]] <- p[[i]][[j]]+ 
      ggplot2::scale_fill_gradientn(colors = rev(pal), limits = c(colorMax,0)) + 
      theme(axis.title.y = element_blank(),axis.text.y=element_blank()) +
      ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank()) + ylab(j) + ggtitle(genes[i]) + theme(axis.title.y = element_blank(),axis.text.y=element_blank(),legend.position="none")
  }
  if(j==last(names(colors))){
    p[[i]][[j]] <- p[[i]][[j]] + theme_void() + theme(axis.title.y = element_blank(),axis.text.y=element_blank(),legend.position="none")
 #to put genome location
  }
}
  ggsave(gridExtra::grid.arrange(grobs = p[[i]], ncol = ncol),file="test.pdf",limitsize=F,height=10)

}








#all genes
all_genes <- unique(obj@ref |> dplyr::filter(seqid != "chrM") |> dplyr::pull(gene_name))
obj@genomeMatrices[["cg_allgenes"]] <- makeWindows(obj, 
                                                     genes = all_genes,
                                                     promoter = FALSE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

#find cluster markers all genes
cluster_allgenes_markers <- FindClusterMarkers(obj, 
                                               matrix = "cg_allgenes", 
                                               genes = protein_coding, group_by="cluster_id",
                                               threads = 200)


cluster_allgenes_markers <- cluster_allgenes_markers |> dplyr::filter(p.val < 0.1)
write.table(cluster_allgenes_markers,file="cluster_allgenes_markers.tsv",sep="\t",col.names=T)
