#!/bin/bash
###############################################################################
#
# This script emulates the submission of an AutoComoChem job to an HPC worker.
# It doesn't do any computational work: it only creates a result folder 
# with the key content neeed for the fitness provider to percieve a succesfully
# completed calculation and retrieve the appropriate data.
# 
# Usage: run this script with no argument or with -h to get the documentation.
#
# Crated: Marco Foscato - October 2023
#
###############################################################################

#
# Function printing the usage instructions
#

function printUsage() {
cat <<EOF 

 This is the emulator of the submission of AutoCompChem jobs. This script 
 doesn not run any computation: it only creates a fake but properly formatted 
 set of output files.

 Usage
 ===== 
 
 <command> -j <jd> [-f <file>[,<files>]] [-cnt,--help]

 Where:

 <command>    is the command you have already used to get this message, and can
              be any link or alias to execute script '$0'.

 Options (if omitted, default values will be used):

 -h (--help)  prints this help message.

 -j FILE      specifies that the job details and parameters are to be taken 
              from file FILE.

 -f FILE[,FILE] specifies pathname of the input file (or many files - comma 
              separated list) to be made available to the job, i.e., 
              copied to the work space.

 --jobname NAME sets the name of this job. By default, we try to set this name
              using either the value of the -j option, or the -f option. If
              this is not possible (e.g., no -j and multiple files given to -f)
              the job will be gven a unique name starting with "accJob-".

 --tcl FILE   requires the generation of a tack-completed label/flag file (i.e.,
              the TCL-file) named FILE. The TCL-file is used to signal that the
              job has completed, and to provide a list of files that you might
              want to recover or transfer somethere else for further processing.
              These files are specified, one by line, in the content of the 
              TCL-file. Use the --tclREGEX to control what files will be added 
              to the TCL-file.

EOF
}


###############################################################################
#
# Main
#
##############################################################################

if [ $# -lt 2 ]
then
    echo " "
    echo " ERROR! You must give me at least the pathname to the AutoCompChem job details/parameters (e.g., -j <pathname>)."
    printUsage
    exit 1;
fi

listInpFiles=""
inpFiles=()
jdFile=""
accArgs="" #NB: we threat them as a single string at this stage
cleanupRegexString=""
populateTCLRegexString=""
jobName="accJob-$(date +%s)"

# Parse application specific command line options
args=("$@")
for ((i=0; i<$#; i++))
do
    arg="${args[$i]}"
    case "$arg" in
        "-j") jdFile=${args[$i+1]};;
        "-f") listInpFiles=${args[$i+1]};;
        "--tcl") tclFile=${args[$i+1]};;
        "--jobname") jobName=${args[$i+1]};;
        *);;
    esac
done

# Process list of input filenames
for inpFile in $(echo "$listInpFiles" | sed 's/,/ /g')
do
    inpFiles+=("$inpFile")
done
echo "Found ${#inpFiles[@]} links to input files."

# Try to get a sensible job name based on input files, or give a unique name
submitDir="$(pwd -P)"
if [ 1 -eq ${#inpFiles[@]} ]
then
    jobName="$(basename "${inpFiles[0]}")"
    jobName="${jobName%.*}"
    submitDir="$(cd "$(dirname "${inpFiles[0]}")" >/dev/null 2>&1 ; pwd -P )"
elif [ 0 -eq ${#inpFiles[@]} ] && [ ! -z "$jdFile" ]
then
    jobName="$(basename "$jdFile")"
    jobName="${jobName%.*}"
    submitDir="$(cd "$(dirname "$jdFile")" >/dev/null 2>&1 ; pwd -P )"
fi
echo "Job name set to '$jobName'"

emulatorDir="$(dirname "$0")"

# Place the stored results in the output folder
case "$jdFile" in
    "P4_DFTXZ.jd")
        echo "Emulating pre-opt DFT job"
        cp "$emulatorDir/data/data_preDFT/Mol1_C_last-DFT.sdf" "$emulatorDir/data/data_preDFT/Mol1_X_last-DFT.sdf" "$emulatorDir/data/data_preDFT/Mol1_Z_last-DFT.sdf" "$emulatorDir/data/data_preDFT/Mol1_DFT.log" .
        cp "$emulatorDir/data/data_preDFT/tc_1702485380485" "$tclFile"
        ;;
   
    "sp_DFT-DZ_all-6.jd")
        echo "Emulating SP DFT job"
        cp  "$emulatorDir/data/data_SP-DFT/Mol1_DFT.log" "$emulatorDir/data/data_SP-DFT/"*last-DFT.sdf .
        cp  "$emulatorDir/data/data_SP-DFT/tc_1702636735466" "$tclFile"
        ;;

    "P4_nativeXTB_all-8.jd")
        echo "Emulating xTB job"
        cp  "$emulatorDir/data/data_xTB/Mol1_XTB.log" "$emulatorDir/data/data_xTB/"*last-xTB.sdf .
        cp  "$emulatorDir/data/data_xTB/tc_1703261267180" "$tclFile"
        ;;

    *)
        echo "ERROR! Unrecognized job. Expecting P4_DFTXZ.jd or sp_DFT-DZ_all-6.jd"
        exit -1
esac

echo "All done!"
exit 0
