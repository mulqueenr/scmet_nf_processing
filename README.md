# scmet_nf_processing
Nextflow processing and singularity container for single-cell methylation processing

250127 Test
initial run on 231
```bash
#first need to make the output dir and the log directory for bcl-convert
outdir="/home/rmulqueen/projects/kismet/data/250130_10xmet_231_nf"
mkdir -p ${outdir}
mkdir -p ${outdir}/logs

cd /data/rmulqueen/projects/kismet/ #move to project directory
rm -rf /data/rmulqueen/projects/kismet/tools/scmet_nf_processing
git clone https://github.com/mulqueenr/scmet_nf_processing /data/rmulqueen/projects/kismet/tools/scmet_nf_processing #pull github repo

nextflow ./tools/scmet_nf_processing/nextflow_running/kismet_processing.groovy \
--flowcellDir /home/rmulqueen/projects/kismet/seq/250127_RM10xMET_RYExome \
--outname 250130_10xMET_231_nftest \
--outdir /home/rmulqueen/projects/kismet/data/250130_10xmet_231_nf \
-with-report \
-resume
```
