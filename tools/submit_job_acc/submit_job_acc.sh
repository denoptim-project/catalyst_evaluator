#!/bin/bash
###############################################################################
#
# This is an example of a script for submitting AutoCompChem job that can 
# contain any number and type sub-jobs.
# Run it with no argument to get the documentation or with option "--help".
# The jobscript is ment for SLURM scheduler. 
# 
# Usage: run this script with no argument or with -h to get the documentation.
#
# Crated: Marco Foscato - June 2024
#
###############################################################################


# Configuration:

# You want to edit these parameters to adapt this script to your system!

# SLURM account paying for cpu time
account=_you_must_change_this_string_

# Pathname to the AutoCompChem JAR archive
ACCJar=_you_must_change_this_path_/autocompchem-2.0.1-jar-with-dependencies.jar
javaForACC=_you_must_change_this_path_/bin/java

# Pathname to file system location where to actually run the job (work space)
scratchdisk=_you_must_change_this_path_

###############################################################################

#
# Function printing the usage instructions
#

function printUsage() {
cat <<EOF 

 This is the command for submission of AutoCompChem jobs. These jobs may include
 one or more other jobs that can run any what-so-ever sorftware.

 Usage
 ===== 
 
 <command> -j <jd> [-f <file>[,<files>]] [...] 

 Where:

 <command>    is the command you have already used to get this message, and can
              be any link or alias to execute script '$0'.

 Options (if omitted, default values will be used):

 -h (--help)  prints this help message.

 -c NN        set the number of cores per node to NN (default: 40).
              As a rule, when submitting other computational tasks from within
              an ACC job, we leave one core to ACC and use the rest for the
              chiled task.

 -t NN        set the walltime to NN hours (default 10 min)

 -j FILE      specifies that the job details and parameters are to be taken 
              from file FILE.

 -f FILE[,FILE] specifies pathname of the input file (or many files - comma 
              separated list) to be made available to the job, i.e., 
              copied to the work space.

 --jobname NAME sets the name of this job. By default, we try to set this name
              using either the value of the -j option, or the -f option. If
              this is not possible (e.g., no -j and multiple files given to -f)
              the job will be gven a unique name starting with "accJob-".

 --accargs STRING use this to provide additional comman line argument when 
              running AutoCompChem. 
        
 --tcl FILE   requires the generation of a tack-completed label/flag file (i.e.,
              the TCL-file) named FILE. The TCL-file is used to signal that the
              job has completed, and to provide a list of files that you might
              want to recover or transfer somethere else for further processing.
              These files are specified, one by line, in the content of the 
              TCL-file. Use the --tclREGEX to control what files will be added 
              to the TCL-file.

 --tclREGEX REGEX[,REGEX]  provide a list of REGEX used to identify which files 
              to list in the TCL-file. REGEX will be use to match the wholepath
              in a "find . -regex" command.

EOF
}

###############################################################################
#
# Main
#
##############################################################################

# Then check input
if [ $# -lt 2 ]
then
    echo " "
    echo " ERROR! You must give me at least the pathname to the AutoCompChem job details/parameters (e.g., -j <pathname>)."
    printUsage
    exit 1;
fi

submitDir="$(pwd -P)"

numCores=40
wallTime=01
listInpFiles=""
inpFiles=()
jdFile=""
accArgs="" #NB: we threat them as a single string at this stage
populateTCLRegexString=""
jobName="accJob-$(date +%s)"

# Parse command line options
args=("$@")
for ((i=0; i<$#; i++))
do
    arg="${args[$i]}"
    case "$arg" in
        "-c") numCores=${args[$i+1]};;
        "-t") wallTime=${args[$i+1]};;
        "-j") jdFile=${args[$i+1]};;
        "-f") listInpFiles=${args[$i+1]};;
        "--tcl") tclFile=${args[$i+1]};;
        "--jobname") jobName=${args[$i+1]};;
        "--accargs") accArgs=${args[$i+1]};;
        "--tclREGEX") populateTCLRegexString=${args[$i+1]};;
        *);;
    esac
done

# Process list of input filenames
for inpFile in $(echo "$listInpFiles" | sed 's/,/ /g')
do
    inpFiles+=("$inpFile")
done
echo "Found ${#inpFiles[@]} links to input files."

echo "Job name set to '$jobName'"

# Process request to change the REGEX, first set the default
outputExclusionRules=()
outputExclusionRules+=("^xtb.err$")
outputExclusionRules+=("^.*property.txt$")
outputExclusionRules+=("^.*\.charges$")
outputExclusionRules+=("^.*\.prop$")
outputExclusionRules+=("^.*\.lastxtb$")
outputExclusionRules+=("^.*\.lastxtberr$")
outputExclusionRules+=("^.*\.engrad$")
outputExclusionRules+=("^.*\.gradient$")
outputExclusionRules+=("^.*\.proc[0-9]+\..*$")
outputExclusionRules+=("^core\..*$")
outputExclusionRules+=("^fort\.[0-9]$")
outputExclusionRules+=("^tmp$")
outputExclusionRules+=("\.tmp$")
outputExclusionRules+=("^[a-zA-Z0-9]+$")
outputExclusionRules+=("^Gau-.*$")

# Define regex to populate TCL file
rulesToListFilesInTCL=()
if [ ! -z "$populateTCLRegexString" ]
then
    rulesToListFilesInTCL=()
    for r in $(echo $populateTCLRegexString | sed 's/,/ /g')
    do
        rulesToListFilesInTCL+=("\"$r\"")
    done
fi

#Prepare job file to be submitted to SLURM
cat<<EOF>$jobName.job
#!/bin/bash -l
##################  AutoCompChem -- JOB-SCRIPT  #####################
#
# SLURM-section
#SBATCH --account=$account
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80000mb
#SBATCH --cpus-per-task=$numCores
#SBATCH --time=$wallTime:10:00
#SBATCH --job-name=$jobName
#SBATCH --output=$jobName.log
#SBATCH --error=$jobName.log
#
#####################################################################

jobName="$jobName"
jdFile="$jdFile"
listInpFiles="$listInpFiles"
accArgs="$accArgs"
submitDir=$submitDir
scratchdisk=$scratchdisk
tclFile=$tclFile

#####################################################################
echo "########### ACC - Job Start ###########"
machine=\`hostname\`
echo "Job started on \$machine at"
date

echo "The job was submitted from \$submitDir"

# Define and Make working space (scratch directory)
workDir=\$scratchdisk/\$LOGNAME/\$SLURM_JOB_ID/\$jobName
echo "The job will use scratch-directory \$workDir"
mkdir -p "\$workDir"

# Copy input files to scratch directory
if [ ! -z "\$jdFile" ]; then cp "\$jdFile" "\$workDir" ; fi
for f in \$(echo "\$listInpFiles" | sed 's/,/ /g')
do
  echo "Copying \$f to work space"
  cp "\$f" "\$workDir"
done

#Move to scratch directory
cd "\$workDir"
   
echo " "
# Call AutoCompChem
echo "********************** Calling ACC ********************"
java -jar $ACCJar --params \$jdFile \$accArgs
echo "********************* done with ACC *******************"
#
ls -Rlt

# We need rules to define what files to list in the TCL file
EOF

# We need to initialize an array of variables
echo "rulesToListFilesInTCL=()" >> $jobName.job
for r in "${rulesToListFilesInTCL[@]}"
do
    echo "rulesToListFilesInTCL+=($r)" >> $jobName.job
done

cat<<EOF>>$jobName.job

#
# function that says goodbye
#
function sayGoodbye() {
    echo "Job finished on \$machine at"
    date
    echo " "
    scontrol show job \$SLURM_JOB_ID
    echo " "
    if [ ! -z "\$tclFile" ] ; then
        touch "\$submitDir/\${tclFile}_pre"
        echo  "\$jobName.log" >> "\$submitDir/\${tclFile}_pre"
        for query in "\${rulesToListFilesInTCL[@]}"
        do
            find * -type f -regex "\$query" -print0 | while read -d \$'\0' f
            do
                echo "Adding \$f to TCL file"
                echo "\$f"  >> "\$submitDir/\${tclFile}_pre"
            done
        done
        mv "\$submitDir/\${tclFile}_pre"  "\$submitDir/\$tclFile"
    fi
    exit 0
}

# We cannot know what are the actual output files, but we need to bring them back home
# So, we copy all, excluding files that match the REGEX rules.
EOF

# We need to initialize an array of variables
echo "outputExclusionRules=()" >> $jobName.job
for r in "${outputExclusionRules[@]}"
do
    echo "outputExclusionRules+=($r)" >> $jobName.job
done

cat<<EOF>>$jobName.job
# If you need to avoid the copying of some crap, just add more exceptions to the list.
# or use the --cleanupREGEX option to reset this list at submission time.

find * -type d | while read d
do
    mkdir -p "\$submitDir"/"\$d"
done

echo " "
echo "  WARNING:"
echo "  We use these REGEX to identify files that are NOT copied back home:"
for regex in "\${outputExclusionRules[@]}"
do
    echo "  -> \$regex"
done
echo " "

find * -type f -print0 | while read -d \$'\0' f
do
    fileName="\$(basename "\$f")"
    extension="\${fineName##*.}"

    ignoreIt=false
    for regex in "\${outputExclusionRules[@]}"
    do
        if [[ "\$fileName" =~ \$regex ]]
        then
            ignoreIt=true
            break;
        fi
    done
    if \$ignoreIt
    then
        echo "File '\$f' is NOT copied back home!"
    else
        cp "\$f" "\$submitDir/\$f"
        if [ \$? -gt 0 ]
        then
            echo " "
            echo "   WARNING! File \$f could not be copied back home for"
            echo "   some problem in executing the cp command."
            echo "   Therefore, I'll not remove the work space directory."
            echo " "
            sayGoodbye
        fi
    fi
done

# Move back to the submitDir's home
cd \$submitDir

if [ -d \$workDir ]
then 
   echo "Removing the scratch-directory with rm -rf"
   rm -rf \$workDir
fi
if [ -d \$workDir ]
then 
   echo "For some strange reason the scratch-directory is still present"
   echo "Trying to remove it using rmdir"
   rmdir \$workDir
   rmdir \$scratchdisk/\$LOGNAME/\$SLURM_JOB_ID/
fi

    
# All done, leave.
sayGoodbye
EOF

chmod o-rwx "$jobName.job"
chmod g-w "$jobName.job"

sbatch $jobName.job

exit 0
