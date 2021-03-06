#!/bin/bash
#SBATCH --job-name="stringtie"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="stringtie_out.%A-%a"
#SBATCH --error="stringtie_out.%A-%a"

#This script takes bam files from HISAT (processed by SAMtools) and performs StringTie assembly and quantification and converts
# data into a format that is readable as count tables for DESeq2 usage

SCRIPT="${0##*/}"  # removes everything until the farthest / from the left - ie basename
DIR="${0%/*}"      # shorthand for dirname
. "${DIR}"/functions.sh    # . is shorthand for "source"

verify_module_installed StringTie
verify_module_installed gffcompare

VALID_STEPS=( assembly merge compare reestimate deseq )
if [ -z "$RUN_STEP" ] ; then
    die "RUN_STEP needs to be set to one or more (':' separated) of ${VALID_STEPS[*]}"
fi

module load StringTie/1.3.3b-foss-2016b
module load gffcompare/0.10.1-foss-2016b

# ! -d means its not a valid directory
[ ! -d "$BAM_DIR" -o ! -d "$OUT_DIR" ] && die "both BAM_DIR and OUT_DIR needs to be specified"
[ ! -f "$GENOME_GFF_FILE" ] && die "GENOME_GFF_FILE needs to be set"

FILES=( 
    _ 
    $( ls -1 "$BAM_DIR"/*.bam )
    )

MERGELIST="${OUT_DIR}/mergelist.txt"
MERGED_GTF="${OUT_DIR}/stringtie_merged.gtf"


## ------------------------------------------------------------------------------------------
# StringTie to assemble transcripts for each sample with the GFF3 annotation file
if echo "${RUN_STEP}" | grep -qw "assembly" ; then
    mandate_slurm_array_use
    i=${FILES[$SLURM_ARRAY_TASK_ID]}
    base=$(basename $i)
    # command structure: $ stringtie <options> -G <reference.gtf or .gff> -o outputname.gtf -l prefix_for_transcripts input_filename.bam
    # -o specifies the output name
    # -G specifies you are aligning with an option GFF or GTF file as well to perform novel transcript discovery 
    # -l Sets <label> as the prefix for the name of the output transcripts. Default: STRG
    # don't use -e here if you want it to assemble any novel transcripts
    run_with_logging "assembly" stringtie -G ${GENOME_GFF_FILE} \
      -o ${OUT_DIR}/${base%.bam}.gtf -l ${base%.bam} -p 10 ${i}
fi

#----------------------------------------------------------------------------------------------------	
# StringTie Merge, will merge all GFF files and assemble transcripts into a non-redundant set of transcripts,
# after which re-run StringTie with -e
# Run StringTie merge, merge transcripts from all samples (across all experiments, not just for a single experiment)
if echo "${RUN_STEP}" | grep -qw "merge" ; then
    [ -n "$SLURM_ARRAY_TASK_ID" ] && die "'merge' step cannot be run with as an array job"

    # generate the merge list
    ls -1 "${OUT_DIR}"/*.gtf > "${OUT_DIR}"/mergelist.txt
    run_with_logging "merge" stringtie --merge -G ${GENOME_GFF_FILE} -o "$MERGED_GTF" "$MERGELIST"
fi
#----------------------------------------------------------------------------
#Re-estimate transcript abundance after merge step
# -e creates more accurate abundance estimations with input transcripts, needed when converting to DESeq2 tables

if echo "${RUN_STEP}" | grep -qw "reestimate" ; then
    mandate_slurm_array_use
    i="${FILES[$SLURM_ARRAY_TASK_ID]}.bam"
    base=$(basename $i)
    run_with_logging "reestimate" stringtie -e -G "$MERGED_GTF" \
              -o "${OUT_DIR}"/${base%.bam}.merge.gtf -p 10 ${i}
fi    
#-----------------------------------------------------------------------------

# Protocol to generate count matrices for genes and transcripts for import into DESeq2 using (prepDE.py) to
# extract this read count information directly from the files generated by StringTie (run with the -e parameter).
# generates two CSV files containing the count matrices for genes and transcripts, Given a list of GTFs,
# which were re-estimated upon merging create sample_list.txt
# Generate count matrices using prepDE.py, prep_DE.py accepts a .txt file listing sample IDs and GTFs paths 

if echo "${RUN_STEP}" | grep -qw "deseq" ; then
    [ -n "$SLURM_ARRAY_TASK_ID" ] && die "'merge' step cannot be run with as an array job"
    :> "${OUT_DIR}"/sample_list.txt #to avoid concatenating file content if run twice this command will 
                       #make a new file each time you run the code.
    
    for i in $(ls -1 "${OUT_DIR}"/*.merge.gtf); do
        base=$(basename $i)
    	echo "${base%.merge.gtf} ${i}" >> "${OUT_DIR}"/sample_list.txt
    done
    (
        cd "${OUT_DIR}"
        python "${DIR}"/prepDE.py -i "${OUT_DIR}"/sample_list.txt
    )
fi    	
