#!/bin/bash
###############################################################################
#
# Script for running one DFT task possibly in parallel with other DFT tasks.
#
# The script implements a lockfile mechanism to prevent overloading: 
# it starts running only if there are no more than N copies of this script 
# already running on this host machine. The value N is hard coded in variable
# 'maxParallelRuns'
#
# @param $1 pathname of the initial SDF file (i.e., source of chemical info)
# @param $2 pathname of the final output SDF file
# @param $3 path to work space
# @param $4 the alphanumerical ID of this job; this code is parsed and used
# to assign proper job-dependent parameters, thus must adhere to the syntax
# <numeridalID>.<label1><label2>
# where:
# <numeridalID> is an integer
# <label1> is a one-character alphanumerical label
# <label2> is a one-character alphanumerical label
# Have a look at the section "State dependent settings" to see
# how these label are interpreted to change setting of this script. 
# Future adaptations of this script are likely to alter that very section
# according to needs.
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
maxParallelRuns="200"
lckMaxWait=50000  # we will wait for a lockfile for up to this number of steps
lockStep=30       # length of waiting step
lockTimeUnit="m"  # time unit: s (seconds), m (minutes), h (hours)

# State dependent settings for states:
# A1: Hoveyda-Grubbs precursor - 1st stereoisomer
# A2: Hoveyda-Grubbs precursor - 2nd stereoisomer
# B1: Ru(IV) alkylidene from productive metathesis - 1st stereoisomer
# B2: Ru(IV) alkylidene from productive metathesis - 2nd stereoisomer
# C1: Ru(IV) MCB from productive metathesis - 1st stereoisomer
# C2: Ru(IV) MCB from productive metathesis - 2nd stereoisomer
# D1: Ru(II) decomposition product - 1st stereoisomer
# D2: Ru(II) decomposition product - 2nd stereoisomer
# E1: Hoveyda-Grubbs precursor with two Cl ligand - 1st stereoisomer
# E2: Hoveyda-Grubbs precursor with two Cl ligand - 2nd stereoisomer
# L0: free ligand

# Labels corresponding to transition states
transitionStateLabels=()

# Labels corresponding to ZTS jobs
ztsLabels=()

# Labels of states that MUST be succesful for the fitness to be completed
requiredStateLabels=("A" "B" "C" "D" "E" "L")

# Geometrical descriptors used to confirm the nature of the species
# NB: SMARTS are used on structures with connectivity defined by DENOPTIM
geomDescLabel="DIST-TS"
geomDescDefB="$geomDescLabel [\$([#6](~[#1])(~[#1])~[Ru])] [\$([#6](~[#6])(~[#6])(~[#1])~[#6](~[#1])(~[#1])~[Ru])]"
geomDescMinB=(2.00)
geomDescMaxB=(2.70)
geomDescDefC="$geomDescLabel [\$([#1]~[Ru])] [\$([#6]~[#1]~[Ru])]
$geomDescLabel [\$([#1]~[Ru])] [Ru]"
geomDescMinC=(0.50 2.20)
geomDescMaxC=(1.50 3.00)
geomDescDefD="$geomDescLabel [\$([#1]~[Ru])] [\$([#6]~[#1]~[Ru])]
$geomDescLabel [\$([#1]~[Ru])] [Ru]"
geomDescMinD=(0.5 1.4)
geomDescMaxD=(1.7 2.2)
#NB: 'geomDescDef' is "" for all other labels

#Charge and spin multiplicity of the entire system in DFT model
chargeA="0"
spinmultA="1"
chargeB="0"
spinmultB="1"
chargeC="0"
spinmultC="1"
chargeD="0"
spinmultD="3"
chargeE="0"
spinmultE="1"
chargeL="0"
spinmultL="1"

#DFT job details (NB: machine depepndent pathnames are dealt with in the
# respectinve machines by the ~/remote/runJob.sh script) 
jdDFTmin="$WORKDIR/P3_DFT-G16_Min_spinRestricted_1.S.jd"
jdDFTminU="$WORKDIR/P3_DFT-G16_Min_spinUnRestricted_1.S.jd"
jdDFTminL="$WORKDIR/P3_DFT-G16_L-only_spinRestricted_1.S.jd"
jdDFTts="$WORKDIR/DFT-NWChem_TS_1.0.jd"
jdDFTzts="$WORKDIR/DFT-NWChem_ZTS_1.S.jd"

#Parameters for HPC (see documentation of interface to HPC)
hpcWorkersList="$REMOTEWORKERBRIDGEHOME/configuration"
lockPointerHPC="/tmp/lockHPC"
pointerHPC="/tmp/lastUsedHPC"
hpcUserNameLst=()
hpcMachineNameLst=()
projDirOnHPCLst=()
keyToHPCLst=()
#
stepWaitForHPC="30"
unitWaitForHPC="m"
maxWaitForHPC="5000"
#

# Minimum number of ZTS refinements of the minimum energy path; below this number the ZTS calculation is considered as failed.
minZtsPaths=10

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
# Portable 'sed -i': avoids dealing with the '-i' option that is not portable.
#
# @param: the expression to apply (e.g., "s/$old/$new/g")
# @param: the pathname to the file to edit
#
function portablesed() {
    sed -e "$1" "$2" > "$2.newFromSed"
    mv "$2.newFromSed" "$2"
}

#
# Function that reads in the list of currently available HPC workers from file.
# The file is read by the RemoteWorkersBridge utils and must have the
# syntax defined in RemoteWorkersBridge.
# @param $1 name of file to read in
#
function getListOfHPCWorkers(){
    if [ ! -f $1 ]
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
# Function calculating the free energy from a properly completed DFT freq. job
# @param $1 SDF file with 'DFT_DATA' property containing:
#           -> Total DFT energy [a.u.]
#           -> Thermal correction to Enthalpy [a.u.]
#           -> Temperature [K]
#           -> Total Entropy  [cal/(mol*K)]
#

function calculateFreeEnergy(){
    vector=$(grep -A1 "DFT_DATA" "$1" | tail -n 1)
    totE=$(echo "$vector" | awk '{print $1}')
    corrH=$(echo "$vector" | awk '{print $2}')
    temp=$(echo "$vector" | awk '{print $3}')
    totS=$(echo "$vector" | awk '{print $5}')
    stdCorr="0.0030188039" #[a.u.]
    cf="627.5095" #conversion factor a.u.-->kcal/(mol*K)
    freeEng=$(bc -l <<< "$totE + $corrH - $temp * ($totS / ($cf*1000)) + $stdCorr")
    addPropertyToSingleMolSDF "FREE_ENERGY" "$freeEng" "$1"
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
inpSDF="$1"
outSDF="$2"
wrkDir="$3"
jobId="$4"
tcl="$wrkDir/$tcl$jobId"

molName=`basename "$inpSDF" .sdf`
# WARNING: the "_outSub" string is hard-coded in fitness_provider.sh
molNum=$(echo "$molName" | awk -F"_outSub" '{print $1}')


#
# Setup Log file
#
log="$jobId.log"
exec > $log
exec 2>&1
echo "Log for single DFT task"
echo "Input: $inpSDF" 
echo "molName: $molName"
echo "molNum:  $molNum"
echo "WorkDir: $wrkDir"
echo "Output: $outSDF"
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
# Parse jobID
#
echo "Partsing jobID..."
jobNum=0
stateLab=0
stereo=0
if [[ "$jobId" == *.* ]] 
then
    pA="${jobId%.*}"
    pB="${jobId#*.}"
    if [ "$pA" == "" ] || [ "$pB" == "" ]
    then
        errMsg="#DFTTask: error parsing jobId"
        abandon "$inpSDF" "$E_FATAL"
    else
        jobNum="$pA"
        if [[ "${#pB}" == 1 ]] 
        then
            stateLab="$pB"
        else
            stateLab="${pB:0:1}"
            stereo="${pB:1:2}"
        fi
    fi
else
    errMsg="#DFTTask: error parsing jobId"
    abandon "$inpSDF" "$E_FATAL"
fi
echo "JobNumber:       $jobNum"
echo "StateLab:        $stateLab"
echo "Stereochemistry: $stereo"


#
# Pick job settings for the type of state (minimum/transition state) 
#
jobDetailsFile="$jdDFTmin"
stateType="MINIMUM"
for tsLab in "${transitionStateLabels[@]}"
do
    if [ "$stateLab" == "$tsLab" ]
    then
        jobDetailsFile="$jdDFTts"
        stateType="TRANSITION STATE"
        break;
    fi
done
for ztsLab in "${ztsLabels[@]}"
do
    if [ "$stateLab" == "$ztsLab" ]
    then
        jobDetailsFile="$jdDFTzts"
        stateType="ZTS"
        break;
    fi
done
isRequired=false
for rsLab in "${requiredStateLabels[@]}"
do
    if [ "$stateLab" == "$rsLab" ]
    then
        isRequired=true
        break;
    fi
done
geomDescDef="none"
geomDescMin=()
geomDescMax=()
charge="0"
spinmult="1"
case "$stateLab" in
    "A")
        charge="$chargeA"
        spinmult="$spinmultA"
        ;;
    "B")
#        geomDescDef="$geomDescDefB"
#        geomDescMin=("${geomDescMinB[@]}")
#        geomDescMax=("${geomDescMaxB[@]}")
        charge="$chargeB"
        spinmult="$spinmultB"
        ;;
    "C")
#        geomDescDef="$geomDescDefC"
#        geomDescMin=("${geomDescMinC[@]}")
#        geomDescMax=("${geomDescMaxC[@]}")
        charge="$chargeC"
        spinmult="$spinmultC"
        ;;
    "D")
#        geomDescDef="$geomDescDefD"
#        geomDescMin=("${geomDescMinD[@]}")
#        geomDescMax=("${geomDescMaxD[@]}")
        charge="$chargeD"
        spinmult="$spinmultD"
        jobDetailsFile="$jdDFTminU"
        ;;
    "E")
#        geomDescDef="$geomDescDefE"
#        geomDescMin=("${geomDescMinD[@]}")
#        geomDescMax=("${geomDescMaxD[@]}")
        charge="$chargeE"
        spinmult="$spinmultE"
        ;;
    "L")
        charge="$chargeL"
        spinmult="$spinmultL"
        jobDetailsFile="$jdDFTminL"
        ;;
esac


#
# WARNING! Multi entry SDF must be expected in case of ZTS jobs.
#
sdfToDFT="$wrkDir/${molNum}_toDFT_$stateLab$stereo.sdf"
cp "$inpSDF" "$sdfToDFT"
addPropertyToMultiMolSDF "CHARGE" "$charge" "$sdfToDFT"
addPropertyToMultiMolSDF "SPIN_MULTIPLICITY" "$spinmult" "$sdfToDFT"


#
# Keep a copy of the input to DFT
#
cp "$sdfToDFT" "$wrkDir/${molNum}_${stateLab}_toDFT.sdf"


# 
# Identify key atoms that will be used to check the geometrical descriptors
#
geomDescAtmIds=()
if [ "$geomDescDef" != "none" ]
then
    geomDescPar="$wrkDir/${molNum}_${stateLab}_accGeomDesc.par"
    geomDescLog="$wrkDir/${molNum}_${stateLab}_accGeomDesc.log"
    echo "Detecting key atoms for definition of geometrical descriptors..."
    echo "VERBOSITY: 1" > "$geomDescPar"
    echo "TASK: MeasureGeomDescriptors" >> "$geomDescPar"
    echo "INFILE: $sdfToDFT" >> "$geomDescPar"
    echo "\$STARTSMARTS: $geomDescDef" >> "$geomDescPar"
    echo "\$END" >> "$geomDescPar"
    echo "ONLYBONDED: true" >> "$geomDescPar"
    # Lauch 
    autocompchem -p "$geomDescPar" > "$geomDescLog"
    # Check outcome
    if [ ! -f "$geomDescLog" ]; then
        errMsg="#MeasureGeomDescriptors: $geomDescLog not found."
        abandon "$sdfToDFT" "$E_OPTERROR"
    fi
    if ! grep -q "Termination status: 0" "$geomDescLog" ; then
        errMsg="#MeasureGeomDescriptors: non-zero exit status from AutoCompChem."
        abandon "$sdfToDFT" "$E_OPTERROR"
    fi
    if ! grep -q "$geomDescLabel" "$geomDescLog" ; then
        errMsg="#MeasureGeomDescriptors: descriptors not found in $geomDescLog."
        abandon "$sdfToDFT" "$E_OPTERROR"
    fi
    # Remember atom indexes for later
    geomDescAtmIds=($(grep "$geomDescLabel" "$geomDescLog" | awk '{print $4}'))
fi


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


#
# Send to HPC
#
locHighFreq="$HIGHFREQUENCY"
for hfLabRoot in "${labelsHighFreqStates[@]}"
do
    if [[ "$stateLab" == "$hfLabRoot" ]]
    then
        locHighFreq=1
        break
    fi
done
if [ "$locHighFreq" == 1 ]
then
    echo "Calling python interface to $hpcMachineName for $stateType $stateLab"
    python "$REMOTEWORKERBRIDGE" -u "$hpcUserName" -m "$hpcMachineName" -i "$sdfToDFT $jobDetailsFile" -p "$projDirOnHPC" -d 5 -t s -x 17280 -k "$keyToHPC"
else
    echo "Calling python interface to $hpcMachineName for $stateType $stateLab"
    python "$REMOTEWORKERBRIDGE" -u "$hpcUserName" -m "$hpcMachineName" -i "$sdfToDFT $jobDetailsFile" -p "$projDirOnHPC" -d "$stepWaitForHPC" -t "$unitWaitForHPC" -x "$maxWaitForHPC" -k "$keyToHPC"
fi
echo "Done with python interface."

#
# Evaluate Output
#
outputDFT="$wrkDir/${molNum}_${stateLab}_DFT.out" #NB: this name is defined in the script running on HPC (the only one permitted to execute from remote)
outputZTSPath="$wrkDir/${molNum}_${stateLab}_DFT-ZTSpath.xyz" #NB: this name is defined in the script running on HPC (the only one permitted to execute from remote)
outputDFTSdf="$wrkDir/${molNum}_${stateLab}_DFT.sdf" #NB: the syntax of outputDFTSdf is recalled for calculating free energy (see calculateFreeEnergy)
mv "./${molNum}_${stateLab}_DFT.out" "$outputDFT"
if [ "$stateType" != "ZTS" ]
then
    echo "Evaluating outcome of DFT (MIN or TS)..."
    evalDFTPar="$wrkDir/${molNum}_${stateLab}_accEvalDFT.par"
    evalDFTLog="$wrkDir/${molNum}_${stateLab}_accEvalDFT.log"
    echo "VERBOSITY: 1" > "$evalDFTPar"
    echo "TASK: EvaluateGaussianOutput" >> "$evalDFTPar"
    echo "INFILE: $outputDFT" >> "$evalDFTPar"
    echo "ANALYSE: ENERGY IMGFREQ" >> "$evalDFTPar"
    echo "QUASIHARM: 100.0" >> "$evalDFTPar"
    echo "IGNOREIM: 5.0" >> "$evalDFTPar"
    echo "LOWESTFREQ: 0.01" >> "$evalDFTPar"
    echo "PRINTLASTGEOMETRY: no_arg" >> "$evalDFTPar"
    echo "OUTFILE: $outputDFTSdf" >> "$evalDFTPar"
    echo "OUTFORMAT: SDF" >> "$evalDFTPar"
    echo "TEMPLATECONNECTIVITY: $sdfToDFT" >> "$evalDFTPar"
    # Launch
    autocompchem -p "$evalDFTPar" > "$evalDFTLog"
    # Check outcome
    if [ ! -f "$evalDFTLog" ]; then
        errMsg="#EvaluationDFTOutput: $evalDFTLog not found."
        abandon "$sdfToDFT" "$E_OPTERROR"
    fi
    if ! grep -q "Termination status: 0" "$evalDFTLog" ; then
        errMsg="#EvaluationDFTOutput: non-zero exit status from AutoCompChem."
        abandon "$sdfToDFT" "$E_OPTERROR"
    fi

    #
    # Analyse geometry
    #
    geometryIsOK="0"
    if [ "$geomDescDef" != "none" ]
    then
        #
        # Use atom index from before
        #
        geomDescIDs=""
        nl='
'
        for atmsRef in "${geomDescAtmIds[@]}"
        do
            #
            # WARNING! For now, assumed to be pairs of indexes (i.e., DISTANCE)
            #
            ids=$(echo "$atmsRef" | sed -e 's/\:/ /g' -e 's/[a-zA-Z]//g')
            geomDescIDs="$geomDescIDs$geomDescLabel $ids${nl}"
        done
        geomDescOutPar="$wrkDir/${molNum}_${stateLab}_accGeomDescOut.par"
        geomDescOutLog="$wrkDir/${molNum}_${stateLab}_accGeomDescOut.log"
        echo "Analysis of refined geometry..."
        echo "VERBOSITY: 1" > "$geomDescOutPar"
        echo "TASK: MeasureGeomDescriptors" >> "$geomDescOutPar"
        echo "INFILE: $outputDFTSdf" >> "$geomDescOutPar"
        echo "\$STARTATOMINDEXES: ${geomDescIDs}\$END" >> "$geomDescOutPar"
        # Lauch
        autocompchem -p "$geomDescOutPar" > "$geomDescOutLog"
        # Check outcome
        if [ ! -f "$geomDescOutLog" ]; then
            errMsg="#MeasureGeomDescriptors: $geomDescOutLog not found."
            abandon "$sdfToDFT" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$geomDescOutLog" ; then
            errMsg="#MeasureGeomDescriptors: non-zero exit status from AutoCompChem."
            abandon "$sdfToDFT" "$E_OPTERROR"
        fi
        if ! grep -q "$geomDescLabel" "$geomDescOutLog" ; then
            errMsg="#MeasureGeomDescriptors: descriptors not found in $geomDescOutLog."
            abandon "$sdfToDFT" "$E_OPTERROR"
        fi
    
    
        #
        # Compare with threshold values
        #
        labels=($(grep "$geomDescLabel" "$geomDescOutLog" | awk '{print $4}'))
        values=($(grep "$geomDescLabel" "$geomDescOutLog" | awk '{print $6}'))
        for iv in $(seq 0 $((${#values[@]}-1))) 
        do
            if [ "${geomDescMin[iv]}_" == "_" ] || [ "${geomDescMax[iv]}_" == "_" ]
            then
                errMsg="WARNING! No min/max limits for ${labels[iv]}. Ignoring!"
            else
                echo "Checking geometrical descriptor '${labels[iv]}': ${values[iv]} (Min=${geomDescMin[iv]} <-> Max:${geomDescMax[iv]})"
                if (( $(bc -l <<< "${values[iv]} < ${geomDescMin[iv]}") ))
                then
                    geometryIsOK="1"
#                    errMsg="#GeomDescriptors: descriptor lower than min. value:'${labels[iv]}' < ${geomDescMin[iv]}"
#                    abandon "$sdfToDFT" "$E_OPTERROR"
                elif  (( $(bc -l <<< "${values[iv]} > ${geomDescMax[iv]}") ))
                then
                    geometryIsOK="1"
#                    errMsg="#GeomDescriptors: descriptor larger than max. value:'${labels[iv]}' > ${geomDescMax[iv]}"
#                    abandon "$sdfToDFT" "$E_OPTERROR"
                fi
            fi
        done
    fi


    #
    # Make final decision on this point
    #
    echo "geometryIsOK: $geometryIsOK"
    echo "energy:"
    grep "$stateType" "$evalDFTLog"
    if [ "$geometryIsOK" == "0" ] && grep -q "$stateType" "$evalDFTLog"
    then
        echo "Succesful completion of DFT job (label: $stateLab)"
        energyData=$(grep "$stateType" "$evalDFTLog" | awk '{print $3,$4,$5,$6,$7,$8}')
        addPropertyToSingleMolSDF "DFT_DATA" "$energyData" "$outputDFTSdf"
        calculateFreeEnergy "$outputDFTSdf"
    else
        echo "Unsuccestul job (label: $stateLab)"
        if $isRequired
        then
            errMsg="#DFTJobs: required job $stateLab was unsuccessful."
	    addPropertyToSingleMolSDF "MOL_ERROR" "$errMsg" "$outputDFTSdf"
        fi
    fi
else
    #
    # Evaluate ZTS 
    # NOTE: ZTS will normally fail, but we only want to run it for a while.
    #
    echo "Evaluating outcome of DFT (ZTS)..."
    numPaths=$(grep -c "^@zts   " "$outputDFT")
    if [ "$numPaths" -gt "$minZtsPaths" ]
    then
        if [ ! -f "$outputZTSPath" ]
        then
            errMsg="#DFTJobs: ZTS path not found."
            abandon "$sdfToDFT" "$E_OPTERROR"
        fi
        echo "Succesful completion of DFT job (label: $stateLab)"
        energyData=$(grep "^@zts   " "$outputDFT" | tail -n1 | awk '{print $8}')


        #
        # Identify best guess for TS and use it as only molecular structure
        #
        #numDigits=$(echo $energyData | wc -c)
        #numDigits=$((numDigits-3))
        #patternTSNRG=$(echo $energyData | cut -c 1-$numDigits)
        patternTSNRG=$(printf "%0.6f\n" "$energyData")
        patternTSNRG=${patternTSNRG/\-/\\-}
        patternTSNRG=${patternTSNRG/\./\\.}
        beadClosestToTS=$(grep "^string:" "$outputDFT" | grep "$patternTSNRG" | tail -n1 | awk '{print $4}')
        echo "Best guess for TS structure is bead #$beadClosestToTS"
        portablesed 's/\(^ [A-Z,a-z]*\)\([0-9]* \)/\1 /g' "$outputZTSPath"
        obabel -ixyz "$outputZTSPath" -f "$beadClosestToTS" -l "$beadClosestToTS" -osdf -O "$outputDFTSdf"
        if [ "$?" != 0 ]
        then
            errMsg="#ExtractDFTOutput: obabel couldn't extract ZTS results."
            abandon "$sdfToDFT" "$E_OPTERROR"
        fi
        addPropertyToSingleMolSDF "DFT_DATA" "$energyData" "$outputDFTSdf"
    else
        echo "Unsuccestul job (label: $stateLab)"
        if [[ $isRequired ]]
        then
            errMsg="#DFTJobs: required job $stateLab was unsuccessful."
            abandon "$sdfToDFT" "$E_OPTERROR"
        fi
    fi
fi


#
# Prepare final output
#
mv "$outputDFTSdf" "$outSDF"


#
# Free-up token: done by funtion 'finish'
#
#releaseLock "$lockIDFile" "$lockFile"


#
# Cleanup
#
if [ "$cleanup" == 0 ]
then
    echo "Cleaning up"
    # we keep a copy of the SDF to DFT file (see "_toDFT.sdf")
    rm -f "$sdfToDFT"
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
writeTCL "$jobId"

exit 0
