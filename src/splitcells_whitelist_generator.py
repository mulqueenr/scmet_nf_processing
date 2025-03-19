"""
singularity shell \
    --bind ~/projects/10x_MET \
    --bind ~/tools/ \
~/singularity/amethyst.sif
"""
#UPDATE THIS TO ALLOW FOR HAMMING DISTANCE IN THE FUTURE

from Bio import SeqIO
import sys
import os
from Bio.Seq import Seq
from scipy.spatial import distance
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--whitelist',default="/cellranger/lib/python/atac/barcodes/737K-cratac-v1.txt.gz")
parser.add_argument('--i7_idx',default="ACTGGTAGAT")
parser.add_argument('--gem_idx',default="gem_idx.txt")
parser.add_argument('--prefix',default="10xmet")
parser.add_argument('--gem_cutoff',default=5000)
parser.add_argument('--outdir',default="./sc_bam")
parser.add_argument('--sequencing_cycles',default="Y151;I10;I16;Y151")
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
gem_idx_pass=pd.DataFrame()
gem_idx_pass['index2']=gem_idx.iloc[0:int(args.gem_cutoff)+1]['corrected_idx']
gem_idx_pass['index2']=[str(Seq(i).reverse_complement()) for i in gem_idx_pass['index2']]
gem_idx_pass['index']=[args.i7_idx for i in gem_idx_pass['index2']]
gem_idx_pass['Sample_ID']=[args.prefix+"_"+i for i in gem_idx_pass['index2']]
gem_idx_pass= gem_idx_pass[['Sample_ID','index','index2']]

with open("samplesheet_gemidx.csv", "w") as f:
    f.write("""[Settings],
CreateFastqForIndexReads,0
OverrideCycles,"""+args.sequencing_cycles+"""
BarcodeMismatchesIndex1,2
BarcodeMismatchesIndex2,0
[Data],
Sample_ID,index,index2
""")

gem_idx_pass.to_csv('samplesheet_gemidx.csv',mode='a',sep=",",header=False,index=False)

