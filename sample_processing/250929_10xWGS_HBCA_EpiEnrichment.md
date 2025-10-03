# Epithelial Enrichment of HBCA samples. Running 10x WGS

less /home/rmulqueen/tools/cellranger-atac-2.1.0/lib/python/tenkit/sample_index.py 

BCMHBCA67L_2h is index Set N G1
SI_NA_G1 = SI_3A_G1 = SI_P03_G1 = ["ATGCGATT", "CATATGCG", "GGATACGA", "TCCGCTAC"]

BCMHBCA102L is index Set N G2
SI_NA_G2 = SI_3A_G2 = SI_P03_G2 = ["ATAACCTA", "CGGTGAGC", "GATCTTAT", "TCCGAGCG"]

```bash

singularity shell \
--bind /data/rmulqueen/projects/kismet \
--bind /home/rmulqueen/ref/ \
--bind /home/rmulqueen:/var/log/bcl-convert \
/home/rmulqueen/singularity/amethyst.sif

source /container_src/container_bashrc

#first need to make the output dir and the log directory for bcl-convert
outdir="/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi"
seqdir="/data/rmulqueen/projects/kismet/seq/250930_VH01788_130_222FTJKNX" 
ref="/home/rmulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"

mkdir -p ${outdir}
mkdir -p ${outdir}/logs

cd ${outdir} #move to project directory
git clone https://github.com/mulqueenr/scmet_nf_processing ./tools/scmet_nf_processing #pull github repo


echo """[Header]
FileFormatVersion,2
[BCLConvert_Settings]
CreateFastqForIndexReads,0
OverrideCycles,Y50;I8;U16;Y50
[BCLConvert_Data]
Sample_ID,index
BCMHBCA67L_2h,ATGCGATT
BCMHBCA67L_2h,CATATGCG
BCMHBCA67L_2h,GGATACGA
BCMHBCA67L_2h,TCCGCTAC
BCMHBCA102L,ATAACCTA
BCMHBCA102L,CGGTGAGC
BCMHBCA102L,GATCTTAT
BCMHBCA102L,TCCGAGCG""" > 250929_10xWGS_HBCA.samplesheet.csv

bcl-convert \
--bcl-input-directory $seqdir \
--output-directory ${outdir}/fastq \
--sample-sheet 250929_10xWGS_HBCA.samplesheet.csv --force

cd ${outdir}/fastq

~/tools/cellranger-atac-2.1.0/cellranger-atac count \
--id=BCMHBCA67L_2h \
--reference=${ref} \
--fastqs=${outdir}/fastq/ \
--sample=BCMHBCA67L_2h \
--project=BCMHBCA67L_2h \
--localcores=100 \
--localmem=1000


#make splitting barcode list from whitelist
outdir="/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi"
indir="/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA67L_2h"
dna_barcodes="${indir}/outs/filtered_peak_bc_matrix/barcodes.tsv"
bam="${indir}/outs/possorted_bam.bam"
awk -F, 'OFS="\t" {print $1,$1}' ${dna_barcodes} > ${outdir}/cell_id.tsv


cd $outdir
#split to chunks of 500 cells for i/o purposes
split -l 500 --numeric-suffixes ${outdir}/cell_id.tsv ${outdir}/cell_id.split.


#run cell splitting for each 500 chunk
for i in cell_id.split* ; do
/home/rmulqueen/.local/bin/sinto filterbarcodes --bam ${bam} --cells $i -p 100 --barcodetag "CB" --outdir ./BCMHBCA67L_2h_scbams ;
done


#generate library complexity based on 10% downsample rates
#count unique chr:start sites
function proj_complexity() {
cellid="BCMHBCA67L_${1::-6}"
for i in $(seq 0.1 0.1 1.0); do
uniq_count=$(samtools view -F 3332 -s $i $1 \
| awk 'OFS="\t"{print $3,$4}' \
| sort \
| uniq -c \
| wc -l)
total_count=$(samtools view -F 3332 -s $i $1 | wc -l)
echo "${cellid},${i},${total_count},${uniq_count}"; done > ${cellid}.projected_metrics.txt
}

export -f proj_complexity
parallel -j 300 proj_complexity ::: $(ls *-1.bam)

```


```bash
~/tools/cellranger-atac-2.1.0/cellranger-atac count \
--id=BCMHBCA102L \
--reference=${ref} \
--fastqs=${outdir}/fastq/ \
--sample=BCMHBCA102L \
--project=BCMHBCA102L \
--localcores=300 \
--localmem=1000

#make splitting barcode list from whitelist
outdir="/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi"
indir="/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA102L"
dna_barcodes="${indir}/outs/filtered_peak_bc_matrix/barcodes.tsv"
bam="${indir}/outs/possorted_bam.bam"
awk -F, 'OFS="\t" {print $1,$1}' ${dna_barcodes} > ${outdir}/cell_id.tsv
cd $outdir

#split to chunks of 500 cells for i/o purposes
split -l 500 --numeric-suffixes ${outdir}/cell_id.tsv ${outdir}/cell_id.split.

#run cell splitting for each 500 chunk
for i in cell_id.split* ; do
/home/rmulqueen/.local/bin/sinto filterbarcodes --bam ${bam} --cells $i -p 100 --barcodetag "CB" --outdir ./BCMHBCA102L_scbams ;
done

#generate library complexity based on 10% downsample rates
#count unique chr:start sites
function proj_complexity() {
cellid="BCMHBCA102L_${1::-6}"
for i in $(seq 0.1 0.1 1.0); do
uniq_count=$(samtools view -F 3332 -s $i $1 \
| awk 'OFS="\t"{print $3,$4}' \
| sort \
| uniq -c \
| wc -l)
total_count=$(samtools view -F 3332 -s $i $1 | wc -l)
echo "${cellid},${i},${total_count},${uniq_count}"; done > ${cellid}.projected_metrics.txt
}

export -f proj_complexity
parallel -j 300 proj_complexity ::: $(ls *-1.bam)
```



Project complexity
```R
library(ggplot2)
library(patchwork)
library(drc)
library(parallel)
library(dplyr)

#project complexity for cells passing QC
dat_67<-do.call("rbind",lapply(list.files(path="/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA67L_2h_scbams",pattern="*projected_metrics.txt",full.names=T),function(x) read.table(x,header=F,sep=",")))
dat_67$method<-"BCMHBCA67L_2h"

dat_102<-do.call("rbind",lapply(list.files(path="/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA102L_scbams",pattern="*projected_metrics.txt",full.names=T),function(x) read.table(x,header=F,sep=",")))
dat_102$method<-"BCMHBCA102L"

dat<-rbind(dat_67,dat_102)
colnames(dat)<-c("sample","downsample_perc","total_frag","uniq_frag","method")

#filter to top 1000 cells
top_cells<- as.data.frame(dat %>% filter(as.numeric(downsample_perc)==1)  %>% slice_max(n=5000,by=method,order_by=uniq_frag))
top_cells<- top_cells$sample

dat<-dat[dat$sample %in% top_cells,]


michaelis_menten_fit<-function(x){
    mm<-dat[dat$sample==x,]
    colnames(mm)<-c("cellid","downsample_perc","S","v","method")
    model.drm <- drm(as.numeric(v) ~ as.numeric(S), data = mm, fct = MM.2())
    km_uniq <- data.frame(S = coef(model.drm)[2])
    km_uniq$v <- predict(model.drm, newdata = km_uniq)
    vmax<-as.numeric(coef(model.drm)[1])
    km<-as.numeric(coef(model.drm)[2])
    current_total_reads<-as.numeric(mm[mm$downsample_perc==1.0,]$S)
    current_uniq_reads<-as.numeric(mm[mm$downsample_perc==1.0,]$v)
    method<-mm[mm$downsample_perc==1.0,]$method
    return(c(x,vmax,km,km_uniq$v,current_total_reads,current_uniq_reads,method))
}

projdat<-as.data.frame(do.call("rbind",mclapply(mc.cores=100,unique(dat$sample),michaelis_menten_fit)))


colnames(projdat)<-c("sample",
"projected_saturated_fragments",
"projected_optimal_seq_effort",
"projected_reads_at_optimal_effort",
"current_total_reads",
"current_uniq_reads",
"method")


projdat$projected_saturated_fragments<-as.numeric(projdat$projected_saturated_fragments)
projdat$current_total_reads<-as.numeric(projdat$current_total_reads)
projdat$projected_reads_at_optimal_effort<-as.numeric(projdat$projected_reads_at_optimal_effort)

plt1<-ggplot(projdat,aes(x=method,y=log10(projected_reads_at_optimal_effort),alpha=0.5,color=method))+
geom_jitter()+
theme_minimal()+
ylim(c(0,8))+
ggtitle(paste0("50% Saturated Library Read Counts:\n","BCMHBCA67L_2h Mean: ",
round(mean(projdat[projdat$method=="BCMHBCA67L_2h",]$projected_reads_at_optimal_effort)),
"\n",
"BCMHBCA102L Mean: ",
round(mean(projdat[projdat$method=="BCMHBCA102L",]$projected_reads_at_optimal_effort))))
ggsave(plt1,file="uniq_reads_per_method.pdf")

projdat %>% group_by(method) %>% summarize(count=n(),reads=mean(projected_reads_at_optimal_effort))


```


#Run copykit
```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```
```R
library(Rsamtools)
library(copykit)
library(BiocParallel)
register(MulticoreParam(progressbar = T, workers = 100), default = T)

dat_67 <- runVarbin("/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA67L_2h_scbams",
                    resolution="1Mb",
                    remove_Y = TRUE,
                    is_paired_end = TRUE)

dat_102 <- runVarbin("/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA102L_scbams",
                    resolution="1Mb",
                    remove_Y = TRUE,
                    is_paired_end = TRUE)

dat<-cbind(dat_67, dat_102)
dat <- runMetrics(dat)
dat <- findAneuploidCells(dat)
dat <- findOutliers(dat)

#add sample from bam file name.
# pdf("sample_heatmap.log2.pdf")
# plotHeatmap(dat, assay="smoothed_bincounts",order_cells="hclust",row_split = 'outlier', n_threads = 40)
# dev.off()

# dat <- knnSmooth(dat)
# dat<- runUmap(dat)
# dat <- findSuggestedK(dat)
# dat <- findClusters(dat)

# pdf("sample_heatmap_and_umap.log2.pdf")
# plotUmap(dat, label = c('subclones'))
# plotHeatmap(dat, assay="smoothed_bincounts",order_cells="hclust",row_split='sample_name',label= c('subclones'), n_threads = 40)
# dev.off()



```