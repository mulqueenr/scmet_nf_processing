"""
Example run:
singularity shell \
    --bind ~/projects/10x_MET \
    --bind ~/tools/ \
    ~/singularity/amethyst.sif

convert_to_h5() {     
    sample_name=${1};     
    python ~/projects/10x_MET/src/premethyst_cgmap_to_h5.py --input $1 ; }
export -f convert_to_h5
parallel -j 100 convert_to_h5 ::: $(ls *CGmap.gz) > cg_cov.txt
"""

#from https://github.com/adeylab/premethyst/blob/main/premethyst_commands/calls2h5.py
import h5py
import os
import numpy as np
import pandas as pd
import argparse
import sys
import re

# to run, type in the command line:
# python /container_src/premethyst_cgmap_to_h5.py 10xmet_CAAGAAAAGTCGCCTG_S939.CGmap.gz

# Define and parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input",type=str,default="10xmet_CAAGAAAAGTCGCCTG_S939.CGmap.gz")
args = parser.parse_args()

file_path=args.input
cell_id=re.sub(".CGmap.gz","",file_path)
output_file=cell_id+".h5.gz"

def read_cgmap_gz(file_path):
    try:
        cov = np.genfromtxt(file_path, delimiter='\t', dtype=[('chr', 'S10'), ('nuc', 'S1'), ('pos', int), ('context', 'S3'), ('context2', 'S3'), ('pct', float),('c',int),('t',int)])
        cov = pd.DataFrame(cov)
        cov['chr']=[x.decode() for x in cov['chr']]
        cg_cov=pd.Series(cov['t']).sum()
        mcg_pct=(pd.Series(cov['c']).sum()/pd.Series(cov['t']).sum())*100
    except TypeError:
        cov = pd.DataFrame()
    print(cell_id+"\t"+str(cg_cov),str(mcg_pct))
    return(cov)

with h5py.File(output_file, 'w') as hdf5_file:
    cg_group = hdf5_file.create_group("CG")
    cov = pd.DataFrame()
    cov_in=read_cgmap_gz(file_path)
    cov = pd.concat([cov,cov_in],ignore_index=True)
    cov=cov.reindex(columns=['chr','pos','pct','t','c'])
    cov['chr']=cov['chr'].astype('S10')
    cov=cov.sort_values(['chr', 'pos'], ascending=[True, True])
    cov=cov.to_records(index=False)
    cg_group.create_dataset(cell_id, data=cov, compression='gzip', compression_opts=9)
