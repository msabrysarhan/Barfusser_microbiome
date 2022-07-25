#!/bin/bash

#source activate /scratch/eurac/sarhan/anaconda3/bin/conda
#conda activate /apps/busco/5.2.2

#module load busco

eval "$(conda shell.bash hook)"
source /home/user/$USER/.bashrc
conda activate /apps/busco/5.2.2


busco -i $1 -f -o $3 -m genome --auto-lineage --out_path $2 --download_path src/busco_downloads #--config /scratch/eurac/sarhan/iceman_snakemake/busco_config.ini




