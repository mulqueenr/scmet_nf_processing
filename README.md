# scmet_nf_processing
Nextflow processing and singularity container for single-cell methylation processing

```bash
#first need to make the output dir and the log directory for bcl-convert
outdir="/home/rmulqueen/projects/kismet/data/250130_10xmet_231_nf"
mkdir -p ${outdir}
mkdir -p ${outdir}/logs

cd /volumes/USR2/Ryan/projects/10x_MET #move to project directory
git clone https://github.com/mulqueenr/scmet_nf_processing ./tools/scmet_nf_processing #pull github repo

nextflow ./tools/scmet_nf_processing/nextflow_running/kismet_processing.groovy \
-with-report \
--flowcellDir /home/rmulqueen/projects/kismet/seq/250127_RM10xMET_RYExome \
--outname 250130_10xMET_231_nftest \
--outdir /home/rmulqueen/projects/kismet/data/250130_10xmet_231_nf
--resume
```