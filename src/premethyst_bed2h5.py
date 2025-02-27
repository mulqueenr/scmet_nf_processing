#from https://github.com/adeylab/premethyst/blob/main/premethyst_commands/calls2h5.py
import h5py
import os
import numpy as np
import pandas as pd
import argparse
import sys

#if len(sys.argv) < 3:
#    print("\npremethyst calls2h5 [input calls folder] [output h5 file prefix]\n")
#    sys.exit(1)

# to run, type in the command line:
# python /container_src/ [input folder path] [output prefix]

# Define and parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_folder")
parser.add_argument("output_prefix")
args = parser.parse_args()

# Input folder and output file paths
folder_path = args.input_folder
#folder_path = "/volumes/USR2/Ryan/projects/metact/240115_RMMM_scalebiotest2/cg_sort_cov/HBCA-16R.2A01/HBCA-16R.2A01.CG.chroms.sort"
#output_file = "/volumes/USR2/Ryan/projects/metact/240115_RMMM_scalebiotest2/cg_sort_cov/HBCA-16R.2A01/HBCA-16R.2A01.CG.chroms.sort/test"
output_file = args.output_prefix + ".h5"

def read_bed_gz(file_path):
    try:
        cov = np.genfromtxt(file_path, delimiter='\t', dtype=[('chr', 'S10'), ('pos', int), ('pos2', int), ('cellid', 'S30'), ('met', 'S1')])
        cov = pd.DataFrame(cov)
        cov['pct']=[100.0 if x.decode('UTF-8')=="Z" else 0.0 for x in cov['met']]
        cov['t']=[0 if x.decode('UTF-8')=="Z" else 1 for x in cov['met']]
        cov['c']=[1 if x.decode('UTF-8')=="Z" else 0 for x in cov['met']]
    except TypeError:
        cov = pd.DataFrame()
    return(cov)

with h5py.File(output_file, 'w') as hdf5_file:
    cg_group = hdf5_file.create_group("CG")
    cov = pd.DataFrame()
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if file_name.endswith(".bed.gz") and "chrY" not in file_name:
            cov_in=read_bed_gz(file_path)
            if cov_in.shape[0]>0:
                cov = pd.concat([cov,cov_in],ignore_index=True)
    cov_list = [x for _, x in cov.groupby(cov['cellid'])]
    for cov_cell in cov_list:
        if cov_cell.shape[0]>10000:
            cell_id=cov_cell['cellid'].unique()[0].decode('UTF-8')
            cov_cell=cov_cell.reindex(columns=['chr','pos','pct','t','c'])
            cov_cell['chr']=cov_cell['chr'].astype('S10')
            cov_cell=cov_cell.sort_values(['chr', 'pos'], ascending=[True, True])
            cov_cell=cov_cell.to_records(index=False)
            cg_group.create_dataset(cell_id, data=cov_cell, compression='gzip', compression_opts=9)                   



""" EXAMPLE RUN
singularity shell \
--bind $HOME \
--bind /volumes/seq/projects/metACT/ \
~/singularity/amethyst.sif

export cg_sort="/volumes/USR2/Ryan/projects/metact/240115_RMMM_scalebiotest2/cg_sort_cov"
export task_cpus=50
mkdir -p ${cg_sort}/h5_files
cd $cg_sort
find $cg_sort -maxdepth 1 -type d -printf '%f\n' | grep "[0-9]$" > cg_cov_folders.txt

premethyst() {
outname=$(basename $1)
indir=$(echo $1"/"${outname}".CG.chroms.sort/")
echo $indir" "$outname
python ~/src/premethyst_calls2h5.py $cg_sort/$indir ${cg_sort}/h5_files/${outname}
}

export -f premethyst
parallel -j ${task_cpus} -a cg_cov_folders.txt premethyst

"""