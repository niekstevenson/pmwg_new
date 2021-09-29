#!/bin/bash
#SBATCH -p normal
#SBATCH -n 24
#SBATCH -t 50:00:00
module load 2020
module load R/4.0.2-intel-2020a

mkdir "$TMPDIR"/output_dir

Rscript --no-save --no-restore --verbose $HOME/pmwg_new/testPMWG_multi.R FALSE > "$TMPDIR"/output_dir/Joint_Eps0.23.txt 2>&1

cp -r "$TMPDIR"/output_dir $HOME
