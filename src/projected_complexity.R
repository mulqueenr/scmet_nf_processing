library(ggplot2)
library(patchwork)
library(drc)
library(parallel)

#project complexity for cells passing QC
rm02_compl<-read.table("/home/rmulqueen/projects/kismet/data/250318_kismetv4_2_231/cellranger_dna/sc_bams/RM02.projected_complexity.tsv",header=F,sep=",")
colnames(rm02_compl)<-c("sample","downsample_perc","total_frag","uniq_frag")
rm02_compl$method="10xWGS"

rm03_compl<-read.table("/home/rmulqueen/projects/kismet/data/250318_kismetv4_2_231/reports/RM03.projected_complexity.tsv",header=F,sep=",")
colnames(rm03_compl)<-c("sample","downsample_perc","total_frag","uniq_frag")
rm03_compl$method="kismet_10N"

rm04_compl<-read.table("/home/rmulqueen/projects/kismet/data/250318_kismetv4_2_231/reports/RM04.projected_complexity.tsv",header=F,sep=",")
colnames(rm04_compl)<-c("sample","downsample_perc","total_frag","uniq_frag")
rm04_compl$method="kismet_10H"

compl<-rbind(rm02_compl,rm03_compl,rm04_compl)


michaelis_menten_fit<-function(x){
    mm<-compl[compl$sample==x,]
    colnames(mm)<-c("sample","downsample_perc","S","v","method")
    model.drm <- drm(v ~ S, data = mm, fct = MM.2())
    km_uniq <- data.frame(S = coef(model.drm)[2])
    km_uniq$v <- predict(model.drm, newdata = km_uniq)
    vmax<-as.numeric(coef(model.drm)[1])
    km<-as.numeric(coef(model.drm)[2])
    current_reads<-as.numeric(mm[mm$downsample_perc==1.0,]$S)
    method<-mm[mm$downsample_perc==1.0,]$method
    return(c(x,vmax,km,km_uniq$v,current_reads,method))
}

projdat<-as.data.frame(do.call("rbind",mclapply(mc.cores=100,unique(compl$sample),michaelis_menten_fit)))

colnames(projdat)<-c("sample",
"projected_total_fragments",
"projected_optimal_seq_effort",
"projected_reads_at_optimal_effort",
"current_reads",
"method")

#filter cells with less than 10000 current reads
projdat$current_reads<-as.numeric(projdat$current_reads)

projdat<-projdat[projdat$current_reads>10000,]
proj_total=mean(as.numeric(projdat$projected_total_fragments),na.rm=T)
proj_optimal=mean(as.numeric(projdat$projected_reads_at_optimal_effort),na.rm=T)
proj_optimal_effort=mean(as.numeric(projdat$projected_optimal_seq_effort),na.rm=T)
proj_current_reads=mean(as.numeric(projdat$current_reads),na.rm=T)

optimal=setNames(nm=unique(projdat$method),unlist(lapply(unique(projdat$method),function(x) mean(as.numeric(proj[projdat$method==x,]$projected_reads_at_optimal_effort),na.rm=T))))
optimal_effort=setNames(nm=unique(projdat$method),unlist(lapply(unique(projdat$method),function(x) mean(as.numeric(proj[projdat$method==x,]$projected_optimal_seq_effort),na.rm=T))))
total=setNames(nm=unique(projdat$method),unlist(lapply(unique(projdat$method),function(x) mean(as.numeric(proj[projdat$method==x,]$projected_total_fragments),na.rm=T))))
current_reads=setNames(nm=unique(projdat$method),unlist(lapply(unique(projdat$method),function(x) mean(as.numeric(proj[projdat$method==x,]$current_reads),na.rm=T))))


plt1<-ggplot(projdat,aes(x=method,y=log10(as.numeric(current_reads)),alpha=0.5,color=method))+
geom_jitter()+
theme_minimal()+
ylim(c(0,7))+
ggtitle(paste0("Current sequencing effort:\n","Mean: ",round(proj_current_reads)))

plt2<-ggplot(projdat,aes(x=method,y=log10(as.numeric(projected_reads_at_optimal_effort)),alpha=0.5,color=method))+
geom_jitter()+
theme_minimal()+
ylim(c(0,7))+
ggtitle(paste0("Unique Fragments Per Cell \nPassing Filter at \nOptimal Sequencing Effort\n","Mean: ",round(proj_optimal)))

plt3<-ggplot(projdat,aes(x=method,y=log10(as.numeric(projected_total_fragments)),alpha=0.5,color=method))+
geom_jitter()+
theme_minimal()+
ylim(c(0,7))+
ggtitle(paste0("Unique Fragments Per Cell \nPassing Filter at \nFull Saturation\n","Mean: ",round(proj_total)))

plt<-plt1+plt2+plt3+plot_layout(widths = c(0.5,0.5,0.5))
ggsave(plt,file="250320_kismet_library_complexity.pdf",width=20)

dat<-read.table("/home/rmulqueen/projects/kismet/data/250318_kismetv4_2_231/metadata.csv",sep=",",header=F)
projdat_rm04<-projdat[startsWith(projdat$sample,prefix="RM04"),]
colnames(dat)<-c("sample","metcg","totcg","percmet")
dat<-merge(dat,projdat_rm04,by="sample")