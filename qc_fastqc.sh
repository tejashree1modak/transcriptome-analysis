#!/bin/bash
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --output="qc_fastqc_out.%A-%a"
#SBATCH --error="qc_fastqc_out.%A-%a"

function die() {
    echo "ERROR: $*" >&2
    exit 1
}

function verify_module_installed() {
    # this function verifies that `module` is installed and that FastQC/0.11.5-... is available
    if ! module avail FastQC | grep -qw FastQC ; then
        die "FastQC module is needed"
    fi
}

function mandate_slurm_array_use() {
    if [ -z "$SLURM_ARRAY_TASK_ID" ] ; then
        # -z means the variable is empty
        die "qc_fastqc.sh needs to be run via slurm as an array task"
    fi
}

if [ ! -d "$SRA_DIR" ]; then
    # ! -d means its not a valid directory
    die "SRA_DIR needs to be specified; it is the directory with .fq files"
fi


files=( 
    _
    $(ls -1 "$SRA_DIR"/*.fq)
)

verify_module_installed

mandate_slurm_array_use

set -e 
module load FastQC/0.11.5-Java-1.8.0_92

if [ $SLURM_ARRAY_TASK_ID -gt 0 -a $SLURM_ARRAY_TASK_ID -lt ${#files[*]} ]; then
    # ${#files[*]} is length of files[] array
    echo "START $(date) : task-id=$SLURM_ARRAY_TASK_ID file=${files[$SLURM_ARRAY_TASK_ID]}"
    fastqc ${files[$SLURM_ARRAY_TASK_ID]}
    echo "END $(date) : task-id=$SLURM_ARRAY_TASK_ID file=${files[$SLURM_ARRAY_TASK_ID]}"
fi

