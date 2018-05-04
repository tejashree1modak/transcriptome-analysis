# --------------------------------------------------------------------------------
# This is a collection of useful functions common to all the scripts
# --------------------------------------------------------------------------------

function die() {
    echo "ERROR: $*" >&2
    exit 1
}

function verify_module_installed() {
    # this function verifies that `module` is installed and that
    # e.g. FastQC/0.11.5-... is available
    # The module name is passed as $1 - i.e first argument

    if ! module avail "$1" | grep -qw "$1" ; then
        die "'$1' module is needed"
    fi
}

function mandate_slurm_array_use() {
    if [ -z "$SLURM_ARRAY_TASK_ID" ] ; then
        # -z means the variable is empty
        die "qc_fastqc.sh needs to be run via slurm as an array task"
    fi
}

function get_paired_sra_files() {
    # $1 will be the directory containing _1.fastq and _2.fastq SRA files
    # only those that have both _1 and _2 will be selected

    local file file2
    for file in $(ls -1 "$1"/*_1.fastq); do 
        file2=$(echo ${file}|sed s/_1/_2/)
        if [ -f "${file2}" ] ; then
            # if _2 exists, then echo the basename of _1 file
            # without the trailing _1.fastq
            file="${file##*/}"
            echo "${file%_1.fastq}"
        fi
    done
}

function run_with_logging() {
    # Runs the command while logging the start and end times
    # $1 - is the message that will be echo'd
    # $2, $3 .... is the command that will be run
    # e.g.
    # run_with_logging foo ls -l /tmp
    #   START foo Fri May  4 17:50:50 EDT 2018: ls -l /tmp             <-- the starting log
    #   lrwxr-xr-x@ 1 root  wheel  11 Jul 15  2017 /tmp -> private/tmp <-- running the command
    #   DONE foo Fri May  4 17:50:50 EDT 2018                          <-- the end log
    local msg="$1"
    shift # now $1 is gone from the argument list, $* is the command to be run

    echo -n "START $msg $(date): "
    echo "$*"

    if [ -z "$DEBUG" ] ; then
        $*      # this runs the command
    fi
    echo "DONE $msg $(date)"
}
