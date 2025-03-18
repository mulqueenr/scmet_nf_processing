library(copykit)
library(BiocParallel)
library(optparse)
#library(colorRamp2)
BiocParallel::bpparam()

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="List of single-cell bam files", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="allcells", 
              help="Prefix of output"),
  make_option(c("-c", "--task_cpus"), type="integer", default=1, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cpu_count=opt$task_cpus
prefix=opt$output_prefix

register(MulticoreParam(progressbar = T, workers = 50), default = T)

genome <- "hg38"
resolution <- "500kb"
min_bincount = 10

# bindings for NSE and data
Chr <- chr <- strand <- GeneID <- NULL
reads_assigned_bins <- reads_duplicates <- reads_total <- NULL

files <-list.files(opt$input_dir, pattern = "*.bbrd.bam", full.names = TRUE, ignore.case = TRUE )
files_names <- basename(gsub(pattern = ".bbrd.bam", "", files))

hg38_grangeslist <- hg38_grangeslist
hg38_rg<-hg38_grangeslist[[paste0("hg38_",resolution)]]
hg38_rg <- as.data.frame(hg38_rg)
rg <- hg38_rg %>%
dplyr::rename(chr = "seqnames") %>%
dplyr::mutate(GeneID = 1:nrow(hg38_rg))
rg <- dplyr::filter(rg,chr != "chrY")
varbin_counts_list_all_fields <-BiocParallel::bplapply(
            files,Rsubread::featureCounts,
            ignoreDup = TRUE,
            countMultiMappingReads = FALSE,
            annot.ext = rg,
            useMetaFeatures = FALSE,
            maxMOp=100, #this is added for methylation alignments
            verbose = TRUE,
            primaryOnly=TRUE, #this is added for methylation alignments
            isPairedEnd = TRUE,
            BPPARAM = bpparam())

names(varbin_counts_list_all_fields) <- files_names
varbin_counts_list <- lapply(varbin_counts_list_all_fields,"[[",1)
varbin_counts_list <- lapply(varbin_counts_list,as.vector)
min_bc <- which(vapply(varbin_counts_list, mean, numeric(1)) < min_bincount)
if (length(min_bc) > 0) {
varbin_counts_list <- varbin_counts_list[-min_bc]
varbin_counts_list_all_fields <- varbin_counts_list_all_fields[-min_bc]
message(length(min_bc), " bam files had less than ", min_bincount," mean bincounts and were removed.")}

# LOWESS GC normalization
message("Performing GC correction.")
varbin_counts_list_gccor <-BiocParallel::bplapply(varbin_counts_list, function(x) {
    gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
    gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
    exp(log(x) - gc_cor_z$y) * median(x)
    },
BPPARAM = bpparam())

varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_gccor), 2)

# filtering low read counts where the sum of bins does not reach more than 0
good_cells <- names(varbin_counts_df[which(colSums(varbin_counts_df) != 0)])
varbin_counts_df <- varbin_counts_df[good_cells]
rg <- rg %>%dplyr::select(-strand, -GeneID)
rg_gr <- GenomicRanges::makeGRangesFromDataFrame(rg,ignore.strand = TRUE,keep.extra.columns = TRUE)
cna_obj <- CopyKit(assays = list(bincounts = varbin_counts_df),rowRanges = rg_gr)

# Adding genome and resolution information to metadata
S4Vectors::metadata(cna_obj)$genome <- genome
S4Vectors::metadata(cna_obj)$resolution <- resolution

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:01 2021
# ADDING READS METRICS TO METADATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:24 2021

varbin_reads_list <- lapply(varbin_counts_list_all_fields,"[[",4)

# saving info and removing columns from list elements
metadata_info_names <- varbin_reads_list[[1]][c(1, 2, 8, 9, 12, 14), 1]
metadata_info_names <-
c(
"reads_assigned_bins",
"reads_unmapped",
"reads_duplicates",
"reads_multimapped",
"reads_unassigned",
"reads_ambiguous"
)

varbin_reads_info <-lapply(seq_along(varbin_reads_list), function(x) {
# RSubread seems to change underlines to dot on some cases
# Have to make more complicated lapply to extract the name of the list
# and guarantee that the cell is properly named
name <- names(varbin_reads_list)[[x]]
df <- varbin_reads_list[[x]][c(1, 2, 8, 9, 12, 14), -1, drop = FALSE]
names(df) <- name
df
})

names(varbin_reads_list) <- names(varbin_counts_list_all_fields)
bam_metrics <- dplyr::bind_cols(varbin_reads_info)

# making sure metrics match varbin_counts_df
bam_metrics <- bam_metrics[good_cells]
rownames(bam_metrics) <- metadata_info_names
bam_metrics <- as.data.frame(t(bam_metrics))

# adding total
reads_tot <- rowSums(bam_metrics)
bam_metrics$sample <- rownames(bam_metrics)
bam_metrics <-
dplyr::relocate(bam_metrics, sample, .before = reads_assigned_bins)
bam_metrics <- bam_metrics %>%
dplyr::mutate(reads_total = reads_tot,percentage_duplicates = round(reads_duplicates / reads_total, 3))

# adding to metadata
SummarizedExperiment::colData(cna_obj) <-S4Vectors::DataFrame(bam_metrics)
colnames(cna_obj) <- names(varbin_counts_df)

dat<-cna_obj
dat<-runVst(dat)
dat<-runSegmentation(dat,alpha = 1e-5, merge_levels_alpha = 1e-5,gamma = 20,name = "segment_ratios")
dat  <- runMetrics(dat)
pdf(paste0(prefix,".qc_metrics.pdf"))
plotMetrics(dat, metric = c("overdispersion", 
                              "breakpoint_count",
                              "reads_total",
                              "reads_duplicates",
                              "reads_assigned_bins",
                              "percentage_duplicates"),label = "reads_total")
dev.off()

#dat <- dat[,colData(dat)$reads_assigned_bins >= 100000]
colData(dat)$sample_name<-prefix
dat<- runMetrics(dat)
dat <- findAneuploidCells(dat)
dat <- findOutliers(dat)

pdf(paste0(prefix,".subclone.heatmap.outlier.pdf"))
plotHeatmap(dat, row_split='outlier',n_threads=cpu_count) #col=colorRamp2(c(-1,0,1),c("blue","white","red"))
dev.off()

# kNN smooth profiles
dat <- knnSmooth(dat)

# Create a umap embedding 
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10

dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK-3, k_subclones=k_clones@metadata$suggestedK+3)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat, label = 'subclones')
dev.off()

pdf(paste0(prefix,".superclone.umap.pdf"))
plotUmap(dat, label = 'superclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
dat <- calcConsensus(dat)
dat <- runConsensusPhylo(dat)
dat <- runPhylo(dat, metric = 'manhattan')
dat <- calcInteger(dat,assay="segment_ratios",method="scquantum")

#col_fun = colorRamp2(c(-0.5, 0, 5), c("blue", "white", "red"))
# Plot a copy number heatmap with clustering annotation
pdf(paste0(prefix,".subclone.heatmap.segment_ratios.pdf"))
plotHeatmap(dat, label = c('reads_assigned_bins','method','sample_name'),assay="segment_ratios",order='hclust',n_threads=cpu_count)
dev.off()

pdf(paste0(prefix,".subclone.heatmap.integer.pdf"))
plotHeatmap(dat, label = c('reads_assigned_bins','method','sample_name'),assay="integer",order='hclust',n_threads=cpu_count)
dev.off()

pdf(paste0(prefix,".subclone.phylo.pdf"))
plotPhylo(dat, label = 'subclones')
dev.off()

dat <- calcInteger(dat, method = 'scquantum', assay = 'smoothed_bincounts')
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat, label = c('reads_total'),order='hclust',n_threads=cpu_count) #,col=col_fun
dev.off()

saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(as.data.frame(dat@colData),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T,row.names=T)

