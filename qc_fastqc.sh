#!/bin/bash
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --output="qc_fastqc_out.%A-%a"
#SBATCH --error="qc_fastqc_out.%A-%a"

SCRIPT=${0##*/}  # removes everything until the farthest / from the left - ie basename
. "$SCRIPT"/functions.sh    # . is shorthand for "source"

if [ ! -d "$SRA_DIR" ]; then
    # ! -d means its not a valid directory
    die "SRA_DIR needs to be specified; it is the directory with .fq files"
fi


files=( 
    _
    $(ls -1 "$SRA_DIR"/*.fq)
)

verify_module_installed FastQC
mandate_slurm_array_use

set -e 
module load FastQC/0.11.5-Java-1.8.0_92

if [ $SLURM_ARRAY_TASK_ID -gt 0 -a $SLURM_ARRAY_TASK_ID -lt ${#files[*]} ]; then
    # ${#files[*]} is length of files[] array
    echo "START $(date) : task-id=$SLURM_ARRAY_TASK_ID file=${files[$SLURM_ARRAY_TASK_ID]}"
    fastqc ${files[$SLURM_ARRAY_TASK_ID]}
    echo "END $(date) : task-id=$SLURM_ARRAY_TASK_ID file=${files[$SLURM_ARRAY_TASK_ID]}"
fi

