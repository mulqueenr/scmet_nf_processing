"""
singularity shell \
    --bind ~/projects/10x_MET \
    --bind ~/tools/ \
~/singularity/amethyst.sif
"""
#UPDATE THIS TO ALLOW FOR HAMMING DISTANCE IN THE FUTURE

import gzip
from Bio import SeqIO
import sys
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from scipy.spatial import distance
import numpy as np
import pandas as pd
import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('--fq1',default="10xmet_231.R1.trim.fastq.gz")
parser.add_argument('--fq2',default="10xmet_231.R2.trim.fastq.gz")
parser.add_argument('--whitelist',default="/volumes/USR2/Ryan/tools/cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz")
parser.add_argument('--gem_idx',default="gem_idx.txt")
parser.add_argument('--gem_cutoff',default=5000)
parser.add_argument('--cores',default=50)
parser.add_argument('--outdir',default="./sc_bam")
args = parser.parse_args()

os.system("mkdir -p "+ args.outdir)  

#read gem counts
gem_idx=pd.read_csv(args.gem_idx,names=['counts',"idx"],sep="\t",dtype={'counts':'Int64','idx':str})
gem_idx=gem_idx.loc[gem_idx['counts']>1000]
gem_idx['corrected_idx'] = [i[1:17] for i in gem_idx['idx'].to_list()]

# read whitelist
whitelist=pd.read_csv(args.whitelist,names=["idx"])#Read whitelist
whitelist['revcomp']=[str(Seq(i).reverse_complement()) for i in whitelist['idx']]

#subset gem counts by whitelist matches (assuming highest numbers have less misscalls)
gem_idx=gem_idx[gem_idx['corrected_idx'].isin(whitelist['revcomp'])]
gem_idx = gem_idx.sort_values(by='counts',ascending=False)
gem_idx_pass=gem_idx.iloc[0:int(args.gem_cutoff)+1]['corrected_idx'].to_list()


#supply gem sequence as x
#open fastq files, correct barcode read names then out fastq 1 and 2  with new read name

def split_fq(x):
    print("Running FastQ split for: " + x)
    with gzip.open(args.fq1, "rt") as handle1, \
        gzip.open(args.fq2, "rt") as handle2, \
        open(args.outdir+"/"+'10xmet.'+x+".R1.fastq", "w") as outfile_fq1, \
        open(args.outdir+"/"+'10xmet.'+x+".R2.fastq", "w") as outfile_fq2:
            for (title1, seq1, qual1), (title2, seq2, qual2) in \
            zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2)):
                    if title1.split(":")[7][1:17] == x: #try direct match first
                        fq1="@%s\n%s\n+\n%s\n" % (title1, seq1, qual1)
                        fq2="@%s\n%s\n+\n%s\n" % (title1, seq2, qual2)
                        outfile_fq1.write(fq1)
                        outfile_fq2.write(fq2)
                    elif round(distance.hamming(Seq(title1.split(":")[7][1:17]), Seq(x))*16) <= 2:
                        title1=":".join(title1.split(":")[0:7])+":r"+x+" "+title1.split(" ")[1]
                        fq1="@%s\n%s\n+\n%s\n" % (title1, seq1, qual1)
                        fq2="@%s\n%s\n+\n%s\n" % (title1, seq2, qual2)
                        outfile_fq1.write(fq1)
                        outfile_fq2.write(fq2)
    os.system('gzip '+ args.outdir+"/"+'10xmet.'+x+".R1.fastq")
    os.system('gzip '+ args.outdir+"/"+'10xmet.'+x+".R2.fastq")

#maybe just do perfect matches to start?
p = Pool(int(args.cores))
with p:
    p.map(split_fq, gem_idx_pass)
