# scmet_nf_processing
Nextflow processing and singularity container for single-cell methylation processing

250127 Test
initial run on 231
```bash
#first need to make the output dir and the log directory for bcl-convert
outdir="/home/rmulqueen/projects/kismet/data/250130_10xmet_231_nf"
mkdir -p ${outdir}
mkdir -p ${outdir}/logs

cd /home/rmulqueen/projects/kismet/ #move to project directory
git clone https://github.com/mulqueenr/scmet_nf_processing ./tools/scmet_nf_processing #pull github repo

nextflow ./tools/scmet_nf_processing/nextflow_running/kismet_processing.groovy \
--flowcellDir /home/rmulqueen/projects/kismet/seq/250127_RM10xMET_RYExome \
--outname 250130_10xMET_231_nftest \
--outdir /home/rmulqueen/projects/kismet/data/250130_10xmet_231_nf \
-with-report \
-resume
```

250318 Test
updated version with no primer dimer, on same 231 GEMs
```bash
#first need to make the output dir and the log directory for bcl-convert
outdir="/home/rmulqueen/projects/kismet/data/250318_kismetv4_2_231"
mkdir -p ${outdir}
mkdir -p ${outdir}/logs

cd /home/rmulqueen/projects/kismet/ #move to project directory
git clone https://github.com/mulqueenr/scmet_nf_processing ./tools/scmet_nf_processing #pull github repo

#index is scalebio homebrew set A9
#TATAGCGCGG #trimming first 2 index cycles (i think that's right?)
#reverse complement then trim off last 2 bases from index sequences
#CCGCGCTATA #CCGCGCTA

nextflow ./tools/scmet_nf_processing/nextflow_running/kismet_processing.groovy \
--flowcellDir /home/rmulqueen/projects/kismet/seq/250318_VH01788_96_AAGHCK2M5 \
--outname RM03_10N \
--outdir ${outdir} \
--i7_idx CCGCGCTA \
--sequencing_cycles "Y50;I8;U16;Y50" \
--max_cpus 300 \
-with-report \
-resume

#index is scalebio homebrew set A10
#CGCAAGTATT
#AATACTTGCG #AATACTTG
nextflow ./tools/scmet_nf_processing/nextflow_running/kismet_processing.groovy \
--flowcellDir /home/rmulqueen/projects/kismet/seq/250318_VH01788_96_AAGHCK2M5 \
--outname RM04_10H \
--outdir ${outdir} \
--i7_idx AATACTTG \
--sequencing_cycles "Y50;I8;U16;Y50" \
-with-report \
-resume

#and regular 10x cellranger-atac for RM02 (SI-NA-B10)
```