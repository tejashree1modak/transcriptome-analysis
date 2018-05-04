#!/bin/bash
#SBATCH --job-name="bbtools"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --output="qc_bbtools_out.%A-%a"
#SBATCH --error="qc_bbtools_out.%A-%a"

SCRIPT=${0##*/}  # removes everything until the farthest / from the left - ie basename
. "$SCRIPT"/functions.sh    # . is shorthand for "source"

verify_module_installed BBMAP
mandate_slurm_array_use

VALID_STEPS=( stat adapter_trim quality_filter force_trim quality_filter )

if [ -z "$RUN_STEP" ] ; then
    die "RUN_STEP needs to be set to one or more (':' separated) of ${VALID_STEPS[*]}"
fi

if [ ! -d "$SRA_DIR" -o ! -d "$OUT_DIR" ]; then
    # ! -d means its not a valid directory
    die "both SRA_DIR and OUT_DIR needs to be specified"
fi

FILES=(
    _
    $( get_paired_sra_files "$SRA_DIR" "fastq" )
)


# This script processes SRA PE end reads with BBtools to find adaptor sequences
#  with BBmerge, and then uses these for adaptor trimming and quality trimming with bbduk.sh.

module load BBMap/37.36-foss-2016b-Java-1.8.0_131

bbmap_dir=$(echo $PATH | tr ':' '\n' | grep -w "BBMap/37.36-foss-2016b-Java-1.8.0_131")
if [ ! -f "${bbmap_dir}/resources/adapters.fa" ] ; then
    die "${bbmap_dir}/resources/adapters.fa not found in BBMap directory"
fi

set -e 

if echo "${RUN_STEP}" | grep -qw "stat" ; then
    left_file="$SRA_DIR/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fastq"
    left_file_basename=$(basename $left_file)
    right_file=$(echo ${left_file}|sed s/_1/_2/)
    file=${left_file_basename%_1.fastq} 

    # [ condition ] && command is a short hand for
    # if [ condition ] ; then
    #  command
    # fi
    [ ! -f "$left_file" -o -f "$right_file" ] && die "Either '$left_file' or '$right_file' not found"

    run_with_logging "stat" bbduk.sh in1="${left_file}" in2="${right_file}" k=23 \
            ref="${bbmap_dir}"/resources/adapters.fa \
            stats="${OUT_DIR}/${file}.stat" out="${OUT_DIR}/${file}.out"
fi

# ------------------------------------------------------------------------------------------
if echo "${RUN_STEP}" | grep -qw "force_trim" ; then
    left_file="$SRA_DIR/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fastq"
    left_file_basename=$(basename $left_file)
    right_file=$(echo ${left_file}|sed s/_1/_2/)

    [ ! -f "$left_file" -o -f "$right_file" ] && die "Either '$left_file' or '$right_file' not found"

    run_with_logging "force_trim" bbduk.sh in1=${left_file}  out1=${OUT_DIR}/${left_file_basename%.fastq}.ftm \
             in2=${right_file} out2=${OUT_DIR}/${left_file_basename%_1.fastq}_2.ftm ftm=5
fi

# ------------------------------------------------------------------------------------------
#
## Trimming of adaptors found in the previous command

if echo "${RUN_STEP}" | grep -qw "adapter_trim" ; then
    left_file="$SRA_DIR/${FILES[$SLURM_ARRAY_TASK_ID]}_1.fastq"
    right_file=$(echo ${left_file}|sed s/_1/_2/)

    [ ! -f "$left_file" -o -f "$right_file" ] && die "Either '$left_file' or '$right_file' not found"
    
    # ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
    # flag -tbo specifies to also trim adaptors based on pir overlap detection using BBMerge 
    # which does not require known adapter sequences)
    # flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)
    run_with_logging "adapter_trim" bbduk.sh in1=${left_file}  out1=${left_file%.fastq}.adp \
         in2=${right_file} out2=${right_file%.fastq}.adp \
         ref="${bbmap_dir}/resources/adapters.fa" \
         ktrim=r k=23 mink=11 hdist=1 tpe tbo
fi


# ------------------------------------------------------------------------------------------

##quality trimming, of both the left and the right sides to get rid of reads that are less than quality 10
if echo "${RUN_STEP}" | grep -qw "quality_trim" ; then
    left_file="$SRA_DIR/${FILES[$SLURM_ARRAY_TASK_ID]}_1.adp"
    right_file=$(echo ${left_file}|sed s/_1/_2/)
    [ ! -f "$left_file" -o -f "$right_file" ] && die "Either '$left_file' or '$right_file' not found"

    run_with_logging "quality_trim" bbduk.sh in1=${left_file}  out1=${left_file%.adp}.trim \
             in2=${right_file} out2=${right_file%.adp}.trim qtrim=rl trimq=20
fi

# ------------------------------------------------------------------------------------------

if echo "${RUN_STEP}" | grep -qw "quality_filter" ; then
    #suffix=ftl
    suffix=fq
    left_file="$SRA_DIR/${FILES[$SLURM_ARRAY_TASK_ID]}_1.adp"
    right_file=$(echo ${left_file}|sed s/_1/_2/)
    [ ! -f "$left_file" -o -f "$right_file" ] && die "Either '$left_file' or '$right_file' not found"
    
    run_with_logging "quality_filter" bbduk.sh in1=${left_file}  out1=${left_file%.adp}.${suffix} \
         in2=${right_file} out2=${right_file%.adp}.${suffix} maq=10
fi

# ------------------------------------------------------------------------------------------

if echo "${RUN_STEP}" | grep -qw "force_trim_right" ; then
    in_suffix=adp
    out_suffix=fq
    left_file="$SRA_DIR/${FILES[$SLURM_ARRAY_TASK_ID]}_1.${in_suffix}"
    right_file=$(echo ${left_file}|sed s/_1/_2/)
    [ ! -f "$left_file" -o -f "$right_file" ] && die "Either '$left_file' or '$right_file' not found"

    run_with_logging bbduk.sh in1=${left_file}  out1=${left_file%.${in_suffix}}.${out_suffix} \
         in2=${right_file} out2=${right_file%.${in_suffix}}.${out_suffix} ftr=115
fi

# ------------------------------------------------------------------------------------------

