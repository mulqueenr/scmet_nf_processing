############################################################################################################################
### BigWig File Export
#' @title exportBigWig
#' @description Export a bigwig track per column of a track stored in genomeMatrices. Can then be imported to IGV or UCSC Genome Browser for interactive viewing.
#'
#' @param obj Amethyst object on which to access the named genomeMatrices object.
#' @param matrix A character string specifying the matrix for bigwig export. Defaults to "cg_cluster_tracks" in accordance with vignettes.
#' @param clusters A list of character strings specifying which clusters to export to bigwig files. If not specified, exports all.
#' @param prefix A character prefix to prepend to output bigwig file names.
#' @return The Amethyst object with a new column in the metadata containing the transferred labels.
#' @export

exportBigWig <- function(obj,
                            matrix = "cg_cluster_tracks",
                            clusters = NULL,
                             prefix = NULL) {
  # Check if the matrix is in the obj and formatted properly.
  if (!(matrix %in% names(obj@genomeMatrices))) {
    stop(paste("Matrix",matrix,"can not be found in genomeMatrices slot."))
  } else if (!(is.null(clusters) && !(clusters %in% colnames(obj@genomeMatrices[[matrix]]))) {
    stop(paste("Not all cluster names can be found in the matrix."))
  }

  #set all cluster names if cluster is NULL
  if (is.null(clusters)) {
    clusters<-colnames(obj@genomeMatrices[[matrix]])
    clusters<-clusters[!(clusters %in% c("chr","start","end"))] #remove chr start end from colnames cluster list.

  # Check if the reduction slot exists
  if (!reduction %in% names(obj@reductions)) {
    stop(paste0("The reduction ", reduction, " does not exist in the object."))
  }

  #just initializing function with simple for loop, but this can be parallelized in future.
  for(i in clusters){
  out_bw<-as.data.frame(obj@genomeMatrices[[matrix]]))
  out_bw<-out_bw[c("chr","start","end",cluster)]
  out_bw<-GRanges(out_bw[complete.cases(out_bw),])
  names(out_bw@elementMetadata)<-"score" #set value to score to be consistent with bigwig default format
  out_bw<-out_bw[unique(findOverlaps(out_bw, type = "any", select = "first"))] #remove overlapping values 
  out_bw <- resize(out_bw, width=1000, fix='start') #resize to avoid 1base overlap
  genome(out_bw)<-"hg38"
  hg38_seq_info<-Seqinfo(genome="hg38")
  seqlengths(out_bw)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_bw@seqnames,]$seqlengths
  print(paste("Saving bigwig for...",cluster))
  if(!is.null(prefix)){
    outname=paste0(prefix,cluster,"_cluster.bw")
  } else {
    outname=paste0(cluster,"_cluster.bw")
  } 
  export(out_bw,con=outname,format='bigWig')
  }
}



############################################################################################################################
### Find Markers based on Metadata Columns
#' @title FindClusterMarkers
#' @description Export a bigwig track per column of a track stored in genomeMatrices. Can then be imported to IGV or UCSC Genome Browser for interactive viewing.
#'
#' @param obj Amethyst object on which to access the named genomeMatrices object.
#' @param matrix A character string specifying the matrix for bigwig export. Defaults to "cg_cluster_tracks" in accordance with vignettes.
#' @param clusters A list of character strings specifying which clusters to export to bigwig files. If not specified, exports all.
#' @param prefix A character prefix to prepend to output bigwig file names.
#' @return The Amethyst object with a new column in the metadata containing the transferred labels.
#' @export
#' 
FindClusterMarkers<-function (obj, matrix, genes, group_by="cluster_id",threads = 1) {
    options(scipen = 3)
    genematrix <- as.matrix(obj@genomeMatrices[[matrix]])
    membership <- dplyr::select(obj@metadata, group_by)
    colnames(membership)<-"cluster_id"
    if (threads > 1) {
        future::plan(future::multicore, workers = threads)
    }
    results <- furrr::future_map(.x = genes, .f = function(gene) {
        gene_results <- list()
        for (id in unique(membership$cluster_id)) {
            members <- rownames(dplyr::filter(membership, cluster_id == 
                id))
            nonmembers <- rownames(dplyr::filter(membership, 
                cluster_id != id))
            tryCatch({
                gene_results[[id]] <- dplyr::mutate(data.frame(p.val = stats::wilcox.test(x = genematrix[gene, 
                  members], y = genematrix[gene, nonmembers])$p.value, 
                  gene = gene, cluster_id = id, mean_1 = mean(genematrix[gene, 
                    members], na.rm = TRUE), mean_2 = mean(genematrix[gene, 
                    nonmembers], na.rm = TRUE)), logFC = log2(mean_2/mean_1), 
                  direction = ifelse(mean_1 > mean_2, "hypermethylated", 
                    "hypomethylated"))
            }, error = function(e) {
                cat("Error processing gene:", gene, "and cluster:", 
                  id, "\n")
                gene_results[[id]] <- NA
            })
        }
        gene_results
    }, .progress = TRUE)
    if (threads > 1) {
        future::plan(future::sequential)
        gc()
    }
    results <- do.call(rbind, lapply(results, function(x) do.call(rbind, 
        x)))
    results <- dplyr::select(dplyr::mutate(dplyr::group_by(results, 
        cluster_id), p.adj = stats::p.adjust(p.val, method = "bonferroni")), 
        p.val, p.adj, everything())
    return(results)
}

#TF motifs enrichment
#Kmeans classification over gene bodies? for types of methylation and correlation to expression
#Gene body summaries by cutting into N chunks and running through window features

