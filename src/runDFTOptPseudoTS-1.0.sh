#!/bin/bash
###############################################################################
#
# Script for running several DFT task possibly
#
# The script implements a lockfile mechanism to prevent overloading: 
# it starts running only if there are no more than N copies of this script 
# already running on this host machine. The value N is hard coded in variable
# 'maxParallelRuns'
#
# @param $1 comma-separated list of pathnames of the initial SDF file (i.e., source of chemical info)
# @param $2 comma-separated list of pathnames of the final output SDF file
# @param $3 path to work space
# @param $4 the numerical ID of this job 
#
# @author Marco Foscato
# @author Jonas Brattebø Ekeli
#
###############################################################################
# Tunable parameters
###############################################################################

# RootName for lockfiles that restrain parallel executions of this script
rootLock="/tmp/runDFTTask"
rootLockID="/tmp/activeDFTTask"
maxParallelRuns="100"
lckMaxWait=500000  # we will wait for a lockfile for up to this number of steps
lockStep=30       # length of waiting step
lockTimeUnit="m"  # time unit: s (seconds), m (minutes), h (hours)

# State dependent settings for states:
# A1: Hoveyda-Grubbs precursor - 1st stereoisomer
# A2: Hoveyda-Grubbs precursor - 2nd stereoisomer
# C1: Ru(IV) MCB from productive metathesis - 1st stereoisomer
# C2: Ru(IV) MCB from productive metathesis - 2nd stereoisomer
# E1: Hoveyda-Grubbs precursor with two Cl ligand - 1st stereoisomer
# E2: Hoveyda-Grubbs precursor with two Cl ligand - 2nd stereoisomer
# L0: free ligand
# X1: TS-guess for productive metathesis - 1st stereoisomer
# X2: TS-guess for productive metathesis - 2nd stereoisomer
# Z1: TS-guess for internal beta-H elimination - 1st stereoisomer
# Z2: TS-guess for internal beta-H elimination - 2st stereoisomer

labels=("X" "Z" "C")
charge=("0" "0" "0")
spinmult=("1" "1" "1")
jobDetailsFile="$WORKDIR/P4_DFTXZ.jd"
# WARNING: sp_DFT-DZ_all-6.jd include a pathname to the Gaussian job details sp_PBEPBE-GD3MBJ_DZ_singlet
gaussianJobDetailsFileX="$WORKDIR/P4_DFTX"
gaussianJobDetailsFileZ="$WORKDIR/P4_DFTZ"
gaussianJobDetailsFileC="$WORKDIR/P4_DFTC"

# Labels of states that MUST be succesful for the fitness to be completed
requiredStateLabels=("X" "Z" "C")

#Parameters for HPC (see documentation of interface to HPC)
hpcWorkersList="$REMOTEWORKERBRIDGEHOME/configuration"
lockPointerHPC="/tmp/lockDFT"
pointerHPC="/tmp/lastUsedDFT"
hpcUserNameLst=()
hpcMachineNameLst=()
projDirOnHPCLst=()
keyToHPCLst=()
#
stepWaitForHPC="30"
unitWaitForHPC="m"
maxWaitForHPC="2000"

# Cleanup 0:doit, 1:don't
cleanup=0

#Exit code for incomplete evaluation of fitness
E_OPTERROR=0 # '0' leads to rejection of this single molecule
E_FATAL=1    # '1' stops the master DEMOPTIM job

# Initialization of job management variables
tcl="tcl-"
lockFile="${rootLock}_notSet.lock"
lockIDFile="${rootLockID}_notSet.lock"
errMsg="Error not assigned."

###############################################################################
# Functions
###############################################################################

#
# Function that reads in the list of currently available HPC workers from file.
# The file is read by the RemoteWorkersBridge utils and must have the 
# syntax defined in RemoteWorkersBridge.
# @param $1 name of file to read in
#
function getListOfHPCWorkers(){
    if [ ! -f "$1" ]
    then
        echo " ERROR! File '$1' not found! "
        exit 1
    fi
    n=0
    while read -r line
    do
        if [[ $line == "#"* ]] || [[ $line == "" ]]
        then
            continue
        fi
        nw=$(echo "$line" | wc -w )
        if [ 4 -ne "$nw" ]
        then
            echo " ERROR! Unexpected syntax in file '$1'. Check line '$line'."
            exit 1
        fi
        hpcUserNameLst+=($(echo "$line" | awk '{print $1}'))
        hpcMachineNameLst+=($(echo "$line" | awk '{print $2}'))
        projDirOnHPCLst+=($(echo "$line" | awk '{print $3}'))
        keyToHPCLst+=($(echo "$line" | awk '{print $4}'))
        n=$((n+1))
    done < <("$REMOTEWORKERBRIDGEHOME/utils/parseConfiguration.sh" -f "$1" -k dft)
    if [ 0 -eq "$n" ]
    then
        echo " ERROR! No HPC worker found in file '$1' "
        exit 1
    fi
    echo "Imported HPC Workers:"
    for i in $(seq 0 $((n-1)))
    do
        echo "${hpcUserNameLst[$i]} ${hpcMachineNameLst[$i]} ${projDirOnHPCLst[$i]} ${keyToHPCLst[$i]}"
    done
}


#
# Function to append a property to an SDF file. Does not overwrite any existing
# property in the SDF-
# @param $1 the name of the property to append
# @param $2 the property value of the property
# @param $3 the SDF file to which append the given property
#
function addPropertyToSingleMolSDF() {
    propName="$1"
    propValue="$2"
    file="$3"
    # we want this to be sed-free...
    nn=$(grep -n '$$$$' "$file" | cut -f1 -d:)
    nn=$((nn-=1))
    uniquestring=$(date +%s)
    head -n "$nn" "$file" > "${file}_tmp$uniquestring"
    echo ">  <$propName>" >> "${file}_tmp$uniquestring"
    echo "$propValue" >> "${file}_tmp$uniquestring"
    echo "" >> "${file}_tmp$uniquestring"
    echo '$$$$' >> "${file}_tmp$uniquestring"
    mv "${file}_tmp$uniquestring" "$file"
}

#
# Function to append a property to a multi entry SDF file.
#  Does not overwrite any existing property in the SDF.
# @param $1 the name of the property to append
# @param $2 the property value of the property
# @param $3 the SDF file
#
function addPropertyToMultiMolSDF() {
    propNameM="$1"
    propValueM="$2"
    fileMOri="$3"
    uniquestring=$(date +%s)
    fileM="$fileMOri$uniquestring"
    cp "$fileMOri" "$fileM"
    awk -v f=$fileM '{print > out}; /\$\$\$\$/{out=f"_tmpAPrp"nm++".sdf"}' nm=1 out="${fileM}_tmpAPrp0.sdf" "$fileM"
    for tmpFile in $(ls "${fileM}_tmpAPrp"*".sdf" )
    do
        addPropertyToSingleMolSDF "$propNameM" "$propValueM" "$tmpFile"
        cat "$tmpFile" >> "${fileM}_allTmpAPrp"
        rm -f "$tmpFile"
    done
    rm -f "$fileM"
    mv "${fileM}_allTmpAPrp" "$fileMOri"
}

#
# Abandon this script due to some error, but first store the latest result
# in the output and append it with information on why we are abandoning
# @param $1 the pathname of the latest result (i.e., an SDF file)
# @param $2 exit status
#
function abandon {
    latestSDF="$1"
    es="$2"
    obabel -isdf "$latestSDF" -osdf -O "$outSDF" --property "MOL_ERROR" "$errMsg"
    obabel -isdf "$outSDF" -osdf -O "$outSDF" --property "EXIT_STATUS" "$es"
    res=$?
    if [ "$res" != 0 ]
    then 
        cp "$latestSDF" "$outSDF"
        addPropertyToSingleMolSDF "MOL_ERROR" "$errMsg" "$outSDF"
        addPropertyToSingleMolSDF "EXIT_SATUS" "$es" "$outSDF"
    fi
    echo " "
    echo "ERROR MESSAGE: "
    echo "$errMsg"
    echo " "
    #NB: the exit status is actually written and read in the outSDF
    # this one if for tests
    exit "$es" 
    #NB: we trap the EXIT signal so that tcl and lock are dealt with
}

#
# Writes the task completed label (TCL), which is an empty file telling to other
# processes that this script has completed its task.
# @param $1 the ID of this job
#
function writeTCL {
    echo "Preparing task completed label $tcl"
    touch "$tcl"
}

#
# Release the lock file and remove related info
#
function releaseLock {
    if [ -f "$lockFile" ]
    then
        echo "Removing lock $lockFile"
        rm -f "$lockFile"
    fi
    if [ -f "$lockIDFile" ]
    then
        echo "Removing lock details in $lockIDFile"
        rm -f "$lockIDFile"
    fi
}

#
# Perform all tasks to be done when exiting the script for whatever reason
#
function finish {
    releaseLock
    if [ ! -f "$tcl" ]
    then
        writeTCL
    fi
}

###############################################################################
# Main
###############################################################################

#
# Ensure that any exit writes a tcl file and remove the lockfile if already
# acquired
#
trap finish EXIT

# 
# Parse arguments
#
if [ "$#" -ne 4 ]
then
    errMsg="#DFTTask: wrong number of arguments"
    echo $errMsg
    tcl="${tcl}PID$$"
    exit "$E_FATAL"
fi
inpSDFs="$1"
inpFiles=()
for inpFile in $(echo "$inpSDFs" | sed 's/,/ /g')
do
    if [ -z "$inpFile" ] || [ ! -f "$inpFile" ]
    then
        echo "ERROR Missing input file"
        exit "$E_FATAL"
    fi
    inpFiles+=("$inpFile")
done

outSDFs="$2"
outFiles=()
for outFile in $(echo "$outSDFs" | sed 's/,/ /g')
do
    outFiles+=("$outFile")
done
# This is used only in case of having to report errors
outSDF="${outFiles[0]}"

wrkDir="$3"
jobId="$4"
tcl="$wrkDir/$tcl${jobId}"

molName=`basename "${inpFiles[0]}" .sdf`
# WARNING: the "_outSub" string is hard-coded in fitness_provider.sh
molNum=$(echo "$molName" | awk -F"_outSub" '{print $1}')

#
# Setup Log file
#
log="$jobId.log"
exec > $log
exec 2>&1
echo "Log for multiple DFT tasks"
echo "Input: ${inpFiles[@]}" 
echo "molName: $molName"
echo "molNum:  $molNum"
echo "WorkDir: $wrkDir"
echo "Output: ${outFiles[@]}"
echo "JobID: $jobId"
echo "TCL: $tcl"

#
# Wait for available token 
#
echo "Acquiring lock..."
for i in $(seq 1 $lckMaxWait)
do
    for j in $(seq 1 $maxParallelRuns)
    do
        lockFile="$rootLock$j.lock"
        (set -o noclobber ; echo > "$lockFile" ) 2> /dev/null
        res=$?
        if [ "$res" == 0 ]
        then
            lockIDFile="$rootLockID$j"
            dateTime=$(date)
            echo "$jobId (PID:$$) owns lock $lockFile from $dateTime" > "$lockIDFile"
            echo "$dateTime"
            echo "Lockfile $j acquired: execution allowed!"
            break 2
        fi
    done
    if [ $i == $lckMaxWait ]
    then
        errMsg="#LockDFT: too log wait for lockfile"
        abandon "$inpSDF" "$E_OPTERROR"
    else
        echo "All $maxParallelRuns lock files already taken:  waiting (i: $i)"
    fi
    sleep "$lockStep$lockTimeUnit"
done

#
# Prepare input with chemical info for DFT job
#
inpFilesToDFT=""
for i in $(seq 0 $((${#labels[@]}-1)))
do
    label="${labels[$i]}"
    # NB name syntax assumed in other places (i.e., cleanup, and scripts on remote)
    sdfToDFT="$wrkDir/${molNum}_${label}_DFT.sdf"
    locInpSDF="noname"
    for inpFile in "${inpFiles[@]}"
    do
        if [[ "$inpFile" == *-${label}* ]]
        then
            locInpSDF="$inpFile"
            break
        fi
    done
    cp "$locInpSDF" "$sdfToDFT"
    addPropertyToMultiMolSDF "CHARGE" "${charge[$i]}" "$sdfToDFT"
    addPropertyToMultiMolSDF "SPIN_MULTIPLICITY" "${spinmult[$i]}" "$sdfToDFT"

    inpFilesToDFT="$inpFilesToDFT $sdfToDFT "

    #
    # Keep a copy of the input to DFT
    #
    cp "$sdfToDFT" "$wrkDir/${molNum}_${label}_toDFT.sdf"
done

# Adding "A" to inpFiles for molName extraction in HPC jobscript
molA1="$wrkDir/${molNum}_outSub-A1.sdf"
if [ -z "$molA1" ] || [ ! -f "$molA1" ]
then
    echo "ERROR Missing sdf file for A1 state '$molA1'."
    exit "$E_FATAL"
fi
cp "$molA1" "$wrkDir/${molNum}_A_DFT.sdf"
molDFTA="$wrkDir/${molNum}_A_DFT.sdf"
inpFilesToDFT="$inpFilesToDFT $molDFTA"

#
# Choose HPC
#
getListOfHPCWorkers "$hpcWorkersList"
idHPC=-1
while true
do
    (set -o noclobber ; echo > "$lockPointerHPC" ) 2> /dev/null
    res=$?
    if [ "$res" == 0 ]
    then
        if [ -f "$pointerHPC" ]
        then
            idHPC=$(cat "$pointerHPC")
            echo "LAST USED: $idHPC"
        fi
        if [ "$idHPC" -ge $((${#hpcMachineNameLst[@]}-1)) ]
        then
            idHPC=-1
        fi
        idHPC=$((idHPC+1))
        echo "$idHPC" > "$pointerHPC"
        rm -f "$lockPointerHPC"
        break
    else
        sleep 1
    fi
done
#
# Un/Comment out this line to impose a specific HPC
#idHPC=2
#
hpcUserName=${hpcUserNameLst[$idHPC]}
hpcMachineName=${hpcMachineNameLst[$idHPC]}
projDirOnHPC=${projDirOnHPCLst[$idHPC]}
keyToHPC=${keyToHPCLst[$idHPC]}
gaussianJobDetailsFile="$gaussianJobDetailsFileX $gaussianJobDetailsFileZ $gaussianJobDetailsFileC"

#
# Send to HPC
#
if [ "$HIGHFREQUENCY" == 1 ]
then
    echo "Calling python interface (HighFreq) to $hpcMachineName for $molNum"
    python "$REMOTEWORKERBRIDGE" -u "$hpcUserName" -m "$hpcMachineName" -i "$inpFilesToDFT $jobDetailsFile $gaussianJobDetailsFile" -p "$projDirOnHPC" -d 5 -t s -x 17280 -k "$keyToHPC" -K md
else
    echo "Calling python interface to $hpcMachineName for $molNum"
    python "$REMOTEWORKERBRIDGE" -u "$hpcUserName" -m "$hpcMachineName" -i "$inpFilesToDFT $jobDetailsFile $gaussianJobDetailsFile" -p "$projDirOnHPC" -d "$stepWaitForHPC" -t "$unitWaitForHPC" -x "$maxWaitForHPC" -k "$keyToHPC" -K md
fi
echo "Done with python interface."


#
# Evaluate Output
#
outputDFT="$wrkDir/${molNum}_DFT.log" #NB: this name is defined in the script running on HPC (the only one permitted to execute from remote)

echo "Evaluating outcome of DFT jobs"
if [ ! -f "$outputDFT" ]; then
    errMsg="#EvaluationDFTOutput: $outputDFT not found."
    abandon "$sdfToDFT" "$E_OPTERROR"
fi
if ! grep -q "Final message: Normal Termination" "$outputDFT" ; then
    errMsg="#EvaluationDFTOutput: Abnormal termination of AutoCompChem."
    abandon "$sdfToDFT" "$E_OPTERROR"
fi

for i in $(seq 0 $((${#labels[@]}-1)))
do
    label="${labels[$i]}"
    echo "Evaluating results for ${label}..."
    outputDFTSdf="$wrkDir/${molNum}_${label}_last-DFT.sdf" #NB: this name is defined in the jobsetails and in the script running on HPC 
    if [ ! -f "$outputDFTSdf" ]; then
        errMsg="#EvaluationDFTOutput: $outputDFTSdf not found."
        abandon "$wrkDir/${molNum}_${label}_DFT.sdf" "$E_OPTERROR"
    fi

    lineInitLog=$(grep -n "File:${molNum}_${label}_DFT.out" "$outputDFT" | awk -F":" '{print $1}')
    linesToKeep=3
    if ! tail -n +"$lineInitLog" "$outputDFT" | head -n "$linesToKeep" | grep -q "Energy"
    then
        errMsg="#EvaluationDFTOutput: energy for $label not found in $outputDFT"
        abandon "$wrkDir/${molNum}_${label}_DFT.sdf" "$E_OPTERROR"
    fi 
    energy=$(tail -n +"$lineInitLog" "$outputDFT" | head -n "$linesToKeep" | grep "Energy" | awk '{print $3}')
    echo "Energy = $energy"

    #
    # Prepare final output
    #
    for ii in $(seq 0 $((${#inpFiles[@]}-1)))
    do
        inpFile="${inpFiles[$ii]}"
        if [[ "$inpFile" == *-${label}* ]]
        then
            cp "$outputDFTSdf" "${outFiles[$ii]}"
            addPropertyToSingleMolSDF "DFT-ENERGY" "$energy" "${outFiles[$ii]}"
            break
        fi
    done
done
echo "Succesful completion of all DFT jobs!"


#
# Cleanup
#
if [ "$cleanup" == 0 ]
then
    echo "Cleaning up"
    # we keep a copy of the SDF to DFT file (see "_toDFT.sdf")
    for label in "${labels[@]}"
    do 
        rm -f "${molNum}_${label}_DFT.sdf"
        rm -f "${molNum}_${label}_last-DFT.sdf"
        rm -f "${molNum}_DFT.log"
    done
    rm -f "$evalDFTPar"
    rm -f "$evalDFTLog"
    rm -f "$extractDFTPar"
    rm -f "$extractDFTLog"
    rm -f "$geomDescPar"
    rm -f "$geomDescLog"
    rm -f "$geomDescOutPar"
    rm -f "$geomDescOutLog"
    echo "Done!"
fi


#
# Create task completed label
#
writeTCL

if ! [ -f $tcl ]
then
    echo "tcl-filei ( $tcl ) was not written... Retrying"
    echo "" > $tcl
fi

exit 0
