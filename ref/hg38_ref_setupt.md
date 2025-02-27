
Set up index with now alt contigs for bsbolt
```bash

wget "https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz"

bsbolt Index -IA \
-G  ~/projects/10x_MET/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-DB ~/projects/10x_MET/ref/hg38_bsbolt

```