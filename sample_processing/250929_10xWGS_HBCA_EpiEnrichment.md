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

function proj_complexity() {
cellid="${1::-4}"
#project complexity bam, just mark duplicates
samtools sort -m 10G -n $1 | \
samtools fixmate -p -m - - | \
samtools sort -m 10G | \
samtools markdup --mode t -s - ${cellid}.mkdup.bam

java -jar /home/rmulqueen/tools/picard.jar \
EstimateLibraryComplexity \
MAX_OPTICAL_DUPLICATE_SET_SIZE=-1 \
TMP_DIR="." \
I=${cellid}.mkdup.bam \
O=${cellid}.complex_metrics.txt

#format a bit
grep "^Unknown" ${cellid}.complex_metrics.txt | \
awk -v cellid=${cellid} 'OFS="," {print cellid,$3,$9,$10}' > ${cellid}.projected_metrics.picard.txt
}


export -f proj_complexity
parallel -j 300 proj_complexity ::: $(ls *BCMHBCA67L.bam)

cat *.projected_metrics.picard.txt > BCMHBCA67L.projected_metrics.txt

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

#use picard tools to project library size
function proj_complexity() {
cellid="${1::-4}"
#project complexity bam, just mark duplicates
samtools sort -m 10G -n $1 | \
samtools fixmate -p -m - - | \
samtools sort -m 10G | \
samtools markdup --mode t -s - ${cellid}.mkdup.bam

java -jar /home/rmulqueen/tools/picard.jar \
EstimateLibraryComplexity \
MAX_OPTICAL_DUPLICATE_SET_SIZE=-1 \
TMP_DIR="." \
I=${cellid}.mkdup.bam \
O=${cellid}.complex_metrics.txt

#format a bit
grep "^Unknown" ${cellid}.complex_metrics.txt | \
awk -v cellid=${cellid} 'OFS="," {print cellid,$3,$9,$10}' > ${cellid}.projected_metrics.picard.txt
}


export -f proj_complexity
parallel -j 300 proj_complexity ::: $(ls *BCMHBCA102L.bam)

cat *.projected_metrics.picard.txt > BCMHBCA102L.projected_metrics.txt

```



Project complexity
```R
library(ggplot2)
library(patchwork)
library(dplyr)

#project complexity for cells passing QC
dat_67<-read.table(file="/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA67L_2h_scbams/BCMHBCA67L.projected_metrics.txt",header=F,sep=",")
dat_67$method<-"BCMHBCA67L_2h"

dat_102<-read.table(file="/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA102L_scbams/BCMHBCA102L.projected_metrics.txt",header=F,sep=",")
dat_102$method<-"BCMHBCA102L"

dat<-rbind(dat_67,dat_102)
colnames(dat)<-c("sample","uniq","perc_dup","size","method")


plt1<-ggplot(dat,aes(x=method,y=log10(uniq),alpha=0.5,color=method))+
geom_jitter()+
theme_minimal()+
ylim(c(0,8))+
ggtitle(paste0("Library Current Seq:\n","BCMHBCA67L_2h Mean: ",
round(mean(dat[dat$method=="BCMHBCA67L_2h",]$uniq)),
"\n",
"BCMHBCA102L Mean: ",
round(mean(dat[dat$method=="BCMHBCA102L",]$uniq))))

plt2<-ggplot(dat,aes(x=method,y=log10(size),alpha=0.5,color=method))+
geom_jitter()+
theme_minimal()+
ylim(c(0,8))+
ggtitle(paste0("Library Size:\n","BCMHBCA67L_2h Mean: ",
round(mean(dat[dat$method=="BCMHBCA67L_2h",]$size)),
"\n",
"BCMHBCA102L Mean: ",
round(mean(dat[dat$method=="BCMHBCA102L",]$size))))
ggsave(plt1/plt2,file="uniq_reads_per_method.pdf")

dat %>% group_by(method) %>% summarize(count=n(),reads=mean(size))


```


#Run copykit
```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```

```R
library(Rsamtools)
library(copykit)
library(BiocParallel)
library(circlize)
library(dplyr)
register(MulticoreParam(progressbar = T, workers = 100), default = T)

dat_67 <- runVarbin("/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA67L_2h_scbams",
                    resolution="500kb",
                    remove_Y = TRUE,
                    is_paired_end = TRUE)

dat_102 <- runVarbin("/data/rmulqueen/projects/kismet/data/250929_10xWGS_HBCA_Epi/BCMHBCA102L_scbams",
                    resolution="500kb",
                    remove_Y = TRUE,
                    is_paired_end = TRUE)

dat<-cbind(dat_67, dat_102)
dat <- runMetrics(dat)
dat <- findAneuploidCells(dat)
dat <- findOutliers(dat)

colData(dat)$patient<-unlist(lapply(base::strsplit(row.names(colData(dat)),"[.]"),"[",2))

#add sample from bam file name.
 pdf("sample_heatmap.log2.pdf")
 plotHeatmap(dat,order_cells="hclust",row_split = 'patient', n_threads = 40)
 dev.off()

 dat <- knnSmooth(dat,k=8)
 dat<- runUmap(dat)
 dat <- findSuggestedK(dat)
 dat <- findClusters(dat)
dat <- calcInteger(dat, method = 'fixed', ploidy_value = 2)
#dat <- calcInteger(dat, method = 'scquantum', assay = 'segment_ratios') #try on segment ratios?


pdf("sample_heatmap_and_umap.log2.pdf")
plotUmap(dat, label = 'patient')
plotUmap(dat, label = 'subclones')
plotHeatmap(dat, 
        col=colorRamp2(breaks=c(-3,-2.5,0,2.5,3),colors=c("darkblue","blue","white","red","darkred")),
        order_cells="hclust",
        label= c('patient','subclones'), n_threads = 100)
plotHeatmap(dat, assay='integer',
        order_cells="hclust",
        label= c('patient','subclones'), n_threads = 100)
dev.off()

 saveRDS(dat,"copykit.hbca.rds")



```