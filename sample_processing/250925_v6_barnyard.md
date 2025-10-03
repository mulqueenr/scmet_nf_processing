
250424 Test
updated version trying to increase efficiency, barnyard

```bash
#first need to make the output dir and the log directory for bcl-convert
outdir="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard"
mkdir -p ${outdir}
mkdir -p ${outdir}/logs

cd /data/rmulqueen/projects/kismet/ #move to project directory
git clone https://github.com/mulqueenr/scmet_nf_processing ./tools/scmet_nf_processing #pull github repo
```

# Optimized Kismet v6 Bisulfite

```bash
#index is scalebio homebrew set A9
#TATAGCGCGG /Volumes/data/rmulqueen/projects/kismet/design/10x_met_design_250417.xlsx
#do rev comp, trim first 2 for i7
#run on barnyard genome

cd /data/rmulqueen/projects/kismet/ #move to project directory
outdir="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard"
cd $outdir

nextflow /data/rmulqueen/projects/kismet/tools/scmet_nf_processing/nextflow_running/kismet_processing.groovy \
--flowcellDir /data/rmulqueen/projects/kismet/seq/250924_VH01788_129_AAG3523M5 \
--outname kismet_optimized \
--outdir $outdir \
--ref_index /home/rmulqueen/ref/hg38_bsbolt \
--i7_idx CCGCGCTA \
--sequencing_cycles "Y151;I8;U16;Y151" \
--max_cpus 100 -resume

#barnyard_bsbolt

cd /data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard/sc_bam

function count_genomes() {
cellid="WGS_${1::-9}"
bam=$1
samtools view -F 3332 $bam \
| awk '{print $3}' \
| awk '{split($1,a,"_"); print a[1]}' \
| sort | uniq -c > ${cellid}.barnyard.count.txt 
}

export -f count_genomes
parallel -j 100 count_genomes ::: $(ls *bbrd.bam)
mkdir ../barnyard_counts
mv *barnyard.count.txt ../barnyard_counts
```

```R
library(ggplot2)
setwd("/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard/barnyard_counts")

read_barnyard_count<-function(x){
sample_name<-strsplit(x,"[.]")[[1]][1]
tmp<-read.table(x)
row.names(tmp)<-tmp$V2
colnames(tmp)<-c("count","genome")
tmp<-list(tmp[tmp$genome=="GRCh38","count"],tmp[tmp$genome=="GRCm39","count"],sample_name)
return(tmp)
}

files_in<-list.files(path="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard/barnyard_counts",pattern="*barnyard.count.txt")

dat<-as.data.frame(do.call("rbind",lapply(files_in,read_barnyard_count)))
colnames(dat)<-c("GRCh38","GRCm39","sample")
dat$GRCh38<-as.numeric(dat$GRCh38)
dat$GRCm39<-as.numeric(dat$GRCm39)

dat$perc_human<-dat$GRCh38/(dat$GRCh38+dat$GRCm39)
dat$species_call<-"mix"
dat[dat$perc_human>0.9,]$species_call<-"human"
dat[dat$perc_human<0.1,]$species_call<-"mouse"

plt<-ggplot(dat,aes(x=GRCh38,y=GRCm39,color=species_call))+geom_point()+theme_minimal()
ggsave(plt,file="kismet_barnyard_readcounts.pdf")
table(dat$species_call)
#human   mix mouse 
#  454  1750   646 

```

# Regular 10x cellranger-atac (SI-NA-D5/6 for different versions)
Run in amethyst.sif

less /home/rmulqueen/tools/cellranger-atac-2.1.0/lib/python/tenkit/sample_index.py 

BY1.1 (optimized kismetv5.2) is index Set N D5
SI_NA_D5 = SI_3A_D5 = SI_P03_D5 = ["AGACGGAT", "CCTTTAGA", "GTCGACTC", "TAGACTCG"]

```bash
singularity shell \
--bind /data/rmulqueen/projects/kismet \
--bind /home/rmulqueen/ref/ \
--bind /home/rmulqueen:/var/log/bcl-convert \
/home/rmulqueen/singularity/amethyst.sif

source /container_src/container_bashrc

ref="/home/rmulqueen/ref/refdata-cellranger-arc-GRCh38-and-GRCm39-2024-A"
outdir="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard"
cd ${outdir}
```

v5.2
```bash
echo """[Header]
FileFormatVersion,2
[BCLConvert_Settings]
CreateFastqForIndexReads,0
OverrideCycles,Y151;I8;U16;Y151
[BCLConvert_Data]
Sample_ID,index
BYv52_WGS,AGACGGAT
BYv52_WGS,CCTTTAGA
BYv52_WGS,GTCGACTC
BYv52_WGS,TAGACTCG""" > 250925_kismetv51_optimized_BYv52_WGS.samplesheet.csv

bcl-convert \
--bcl-input-directory /data/rmulqueen/projects/kismet/seq/250924_VH01788_129_AAG3523M5 \
--output-directory ${outdir}/BYv52_WGS \
--sample-sheet 250925_kismetv51_optimized_BYv52_WGS.samplesheet.csv --force

cd ${outdir}/BYv52_WGS

~/tools/cellranger-atac-2.1.0/cellranger-atac count \
--id=BYv52_WGS \
--reference=${ref} \
--fastqs=${outdir}/BYv52_WGS \
--sample=BYv52_WGS \
--project=BYv52_WGS \
--localcores=300 \
--localmem=1000

#make splitting barcode list from whitelist
outdir="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard"
indir="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard/BYv52_WGS/BYv52_WGS/"
dna_barcodes="${indir}/outs/filtered_peak_bc_matrix/barcodes.tsv"
bam="${indir}/outs/possorted_bam.bam"

awk -F, 'OFS="\t" {print $1,$1}' ${dna_barcodes} > ${outdir}/cell_id.tsv

cd $outdir
#split to chunks of 500 cells for i/o purposes
split -l 500 --numeric-suffixes ${outdir}/cell_id.tsv ${outdir}/cell_id.split.

#run cell splitting for each 500 chunk
for i in cell_id.split* ; do
/home/rmulqueen/.local/bin/sinto filterbarcodes --bam ${bam} --cells $i -p 100 --barcodetag "CB" --outdir ./sc_dna_bam_v52 ;
done

#generate library complexity based on 10% downsample rates
#count unique chr:start sites
function proj_complexity() {
cellid="WGS_${1::-6}"
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

#excluding reads that meet any below conditions:
#read unmapped (0x4)
#not primary alignment (0x100)
#read is PCR or optical duplicate (0x400)
#supplementary alignment (0x800)

```

v6

BY2.1 (kismetv6) is index Set N D6
SI_NA_D6 = SI_3A_D6 = SI_P03_D6 = ["ATGCCAAA", "CCTTATCG", "GAAGTCTT", "TGCAGGGC"]


```bash
echo """[Header]
FileFormatVersion,2
[BCLConvert_Settings]
CreateFastqForIndexReads,0
OverrideCycles,Y151;I8;U16;Y151
[BCLConvert_Data]
Sample_ID,index
BYv6_WGS,ATGCCAAA
BYv6_WGS,CCTTATCG
BYv6_WGS,GAAGTCTT
BYv6_WGS,TGCAGGGC""" > 250925_kismetv51_optimized_BY6_WGS.samplesheet.csv

bcl-convert \
--bcl-input-directory /data/rmulqueen/projects/kismet/seq/250924_VH01788_129_AAG3523M5 \
--output-directory ${outdir}/BYv6_WGS \
--sample-sheet 250925_kismetv51_optimized_BY6_WGS.samplesheet.csv --force

cd ${outdir}/BYv6_WGS
~/tools/cellranger-atac-2.1.0/cellranger-atac count \
--id=BY6_WGS \
--reference=${ref} \
--fastqs=${outdir}/BYv6_WGS/ \
--sample=BY6_WGS \
--project=BY6_WGS \
--localcores=300 \
--localmem=1000

#make splitting barcode list from whitelist
outdir="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard"
indir="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard/BYv6_WGS/BY6_WGS"
dna_barcodes="${indir}/outs/filtered_peak_bc_matrix/barcodes.tsv"
bam="${indir}/outs/possorted_bam.bam"

awk -F, 'OFS="\t" {print $1,$1}' ${dna_barcodes} > ${outdir}/cell_id.tsv

cd $outdir
#split to chunks of 500 cells for i/o purposes
split -l 500 --numeric-suffixes ${outdir}/cell_id.tsv ${outdir}/cell_id.split.

#run cell splitting for each 500 chunk
for i in cell_id.split* ; do
/home/rmulqueen/.local/bin/sinto filterbarcodes --bam ${bam} --cells $i -p 100 --barcodetag "CB" --outdir ./sc_dna_bam_v6 ;
done

#generate library complexity based on 10% downsample rates
#count unique chr:start sites
function proj_complexity() {
cellid="WGS_${1::-6}"
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
dat_kismet<-do.call("rbind",lapply(list.files(path="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard/reports/projected_size",pattern="*projected_metrics.txt",full.names=T),function(x) read.table(x,header=F,sep=",")))
dat_kismet$method<-"kismet"

dat_wgs<-do.call("rbind",lapply(list.files(path="/data/rmulqueen/projects/kismet/data/250925_kismetv6_barnyard/sc_dna_bam_v6/",pattern="*projected_metrics.txt",full.names=T),function(x) read.table(x,header=F,sep=",")))
dat_wgs$method<-"wgs_v6"

dat<-rbind(dat_kismet,dat_wgs)
colnames(dat)<-c("sample","downsample_perc","total_frag","uniq_frag","method")

#filter to kismet top 1000 cells
top_cells<- as.data.frame(dat %>% filter(as.numeric(downsample_perc)==1)  %>% slice_max(n=1000,by=method,order_by=uniq_frag))
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

plt1<-ggplot(projdat,aes(x=method,y=log10(projected_reads_at_optimal_effort*2),alpha=0.5,color=method))+
geom_jitter()+
theme_minimal()+
ylim(c(0,7))+
ggtitle(paste0("50% Saturated Library Read Counts:\n","Kismet Mean: ",
round(mean(projdat[projdat$method=="kismet",]$projected_reads_at_optimal_effort)*2),
"\nWGS Mean: ",
round(mean(projdat[projdat$method=="wgs_v6",]$projected_reads_at_optimal_effort)*2)))
ggsave(plt1,file="uniq_reads_per_method.pdf")

plt2<-ggplot(projdat,aes(x=method,y=log10(projected_saturated_fragments*2),alpha=0.5,color=method))+
geom_jitter()+
theme_minimal()+
ylim(c(0,7))+
ggtitle(paste0("100% Saturated Library Read Counts:\n","Kismet Mean: ",
round(mean(projdat[projdat$method=="kismet",]$projected_saturated_fragments)*2),
"\nWGS Mean: ",
round(mean(projdat[projdat$method=="wgs_v6",]$projected_saturated_fragments)*2)))
ggsave(plt2,file="uniq_reads_per_method_sat.pdf")

projdat %>% group_by(method) %>% summarize(count=n(),reads=mean(projected_saturated_fragments))


```

