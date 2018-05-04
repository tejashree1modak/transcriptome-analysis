#!/bin/bash
#SBATCH --job-name="hisat"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="align_hisat_out.%A-%a"
#SBATCH --error="align_hisat_out.%A-%a"

# This script requires SRA_DIR, OUT_DIR and GENOME_FILE to be set

SCRIPT=${0##*/}  # removes everything until the farthest / from the left - ie basename
DIR="${0%/*}"      # shorthand for dirname
. "${DIR}"/functions.sh    # . is shorthand for "source"

verify_module_installed HISAT2
verify_module_installed SAMtools
mandate_slurm_array_use

VALID_STEPS=( gnome_index hisat )
if [ -z "$RUN_STEP" ] ; then
    die "RUN_STEP needs to be set to one or more (':' separated) of ${VALID_STEPS[*]}"
fi

module load HISAT2/2.0.4-foss-2016b   
module load SAMtools/1.3.1-foss-2016b

if [ ! -d "$SRA_DIR" -o ! -d "$OUT_DIR" ]; then
    # ! -d means its not a valid directory
    die "both SRA_DIR and OUT_DIR needs to be specified"
fi

FILES=(
    _
    $( get_paired_sra_files "$SRA_DIR" "fq" )
)

[ ! -f "$GENOME_FILE"  ] && die "GENOME_FILE needs to be set"
GENOME_INDEX="${GENOME_FILE%.*}"


#-------------------------------------------------------------------------------------------

if echo "${RUN_STEP}" | grep -qw "genome_index" ; then

    (
        cd "${GENOME_FILE%/*}"      # shorthand for dirname
        run_with_logging "gnome_index" hisat2-build -f ${GENOME_FILE}  "${GENOME_INDEX}"
    ) # creating a subshell with ( ) means we don't change current directory
fi

# ------------------------------------------------------------------------------------------
# hisat2 align paired end reads to reference genome
if echo "${RUN_STEP}" | grep -qw "hisat" ; then
    left_file="$SRA_DIR/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fq"
    right_file=$(echo ${left_file}|sed s/_1/_2/)
    file=$(basename ${left_file%_1.fq}) 

    [ ! -f "$left_file" -o -f "$right_file" ] && die "Either '$left_file' or '$right_file' not found"
    run_with_logging "hisat" hisat2 --dta -x "${GENOME_INDEX}" -1 "${left_file}" -2 "${right_file}" \
                    -S "${OUT_DIR}/${file}.sam"
fi
