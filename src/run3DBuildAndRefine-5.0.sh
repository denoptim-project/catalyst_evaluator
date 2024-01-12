#!/bin/bash
###############################################################################
#
# Script for construction of the 3D models and relaxation of the ligand 
# conformation, but does not check for atom clashes. It also considers the 
# generation and refinement of modified geometries deriving from the application
# of a given Cartesian move.
#
# The script implements a lockfile mechanism to prevent overloading: 
# it starts running only if there are no more than N copies of this script 
# already running on this host machine.
#
# @param $1 pathname of the initial SDF file (i.e., source of DENOPTIMGraph)
# @param $2 pathname of the final output SDF file
# @param $3 path to work space
# @param $4 the alphanumerical ID of this job; this code is parsed and used
# for to assign proper job-dependent parameters, thus must adhere to the syntax
# <numeridalID>.<label1><label2>
# where:
# <numeridalID> is an integer
# <label1> is a one-character alphanumerical label
# <label2> is a one-character alphanumerical label
# Have a look at the section "Setup the job type-dependent parameters" to see
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

#Tinker (with extensions)
pathToTinkerBin="$TINKERBIN"

#Label and extension of input file
inputLabelExt="_core-ligGraph.sdf" #NB: this extension is defined also in other scripts

# RootName for lockfiles that restrain parallel executions of this script
rootLock="/tmp/runConfSearchOrRefine"
rootLockID="/tmp/activeConfSearchOrRefine"
lckMaxWait=10000   # we will wait for a lockfile for up to this number of steps
lockStep=30       # length of waiting step
lockTimeUnit="s"  # time unit: s (seconds), m (minutes), h (hours)

# GraphEditor
editGraph=1 # 0:doit, 1:don´t
# Graph editing tasks
grpEdtRemoveXLigands="$WORKDIR/removeXligands"
grpEdtReplaceLigands="$WORKDIR/replaceLigandsAndScaffoldEditingTask"

# Fragment space
compMatrix="$WORKDIR/CPMap_core-lig.par"
#Defined below based on the input
libFrags="notSet"
#Not needed, but defined below based on the input
libCaps="notSet"
#Defined below based on the input
libScaff="notSet"
libScaffA1="$WORKDIR/scaff_A1_v3.sdf"
libScaffA2="$WORKDIR/scaff_A2_v3.sdf"
libScaffB1="$WORKDIR/scaff_B1_v3.sdf"
libScaffC1="$WORKDIR/scaff_C1_v3.sdf"
libScaffD1="$WORKDIR/scaff_D1_v3.sdf"
libScaffE1="$WORKDIR/scaff_E1_v3.sdf"
libScaffB2="$WORKDIR/scaff_B2_v3.sdf"
libScaffC2="$WORKDIR/scaff_C2_v3.sdf"
libScaffC3="$WORKDIR/scaff_C3_v3.sdf"
libScaffD2="$WORKDIR/scaff_D2_v3.sdf"
libScaffE2="$WORKDIR/scaff_E2_v3.sdf"
libScaffX1="$WORKDIR/scaff_X1_v3.sdf"
libScaffY1="$WORKDIR/scaff_Y1_v3.sdf"
libScaffZ1="$WORKDIR/scaff_Z1_v3.sdf"
libScaffX2="$WORKDIR/scaff_X2_v3.sdf"
libScaffY2="$WORKDIR/scaff_Y2_v3.sdf"
libScaffZ2="$WORKDIR/scaff_Z2_v3.sdf"
#Define the Cartesian move 
crtMove="notSet"
crtMoveX1="$WORKDIR/crtMove_X1.dat"
crtMoveX2="$WORKDIR/crtMove_X2.dat"
crtMoveY1="$WORKDIR/crtMove_Y1.dat"
crtMoveY2="$WORKDIR/crtMove_Y2.dat"
crtMoveZ1="$WORKDIR/crtMove_Z1.dat"
crtMoveZ2="$WORKDIR/crtMove_Z2.dat"
crtMvRef="notSet"
crtMvRefX1="$WORKDIR/crtMvRef_X1.sdf"
crtMvRefX2="$WORKDIR/crtMvRef_X2.sdf"
crtMvRefY1="$WORKDIR/crtMvRef_Y1.sdf"
crtMvRefY2="$WORKDIR/crtMvRef_Y2.sdf"
crtMvRefZ1="$WORKDIR/crtMvRef_Z1.sdf"
crtMvRefZ2="$WORKDIR/crtMvRef_Z2.sdf"
#Define the scaling factors for the Cartesian move
scalingMove="notSet"
scalingMoveX1="-3 4" 
scalingMoveX2="-3 4" 
scalingMoveY1="-5 6" 
scalingMoveY2="-5 6"
scalingMoveZ1="-6 4"
scalingMoveZ2="-6 4"

# state (i.e., A,B,C...) specific molecular parameters (may need to add chargeA, chargeB,...)
charge=0
spinMult=1

# 3D builder
tinkerForceField="$WORKDIR/uff_vdw.prm"
rotSpaceDef="$WORKDIR/rotatableBonds-RefineOnly"
tinkerKeyFile="$WORKDIR/build_uff_Refine.key"
tinkerSubmitFile="$WORKDIR/submit_Refine"

# Option for checking for atom clashes
checkAtomClashes=1 # 0:check, 1:don't check

# Spartan conformation-constrained geometry optimization
#sprtGOFF="$WORKDIR/params.MMFF94_CS-2.2" we use a global params file ~/PARAMS.MMFF.DENOPTIM
sprtGOFrozen="
[*]~[Ru]
[\$([#6,#1]~[#6]~[Ru]);!\$([#6,#1]~[#6](~[Ru])~[#7])]
[\$([*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[#6](~[Ru])~[#7]);!\$([*]~[#7]~[#6](~[Ru])~[#6,#7])]
[\$([*]~[*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[*]~[#6]~[#7]);!\$([*]~[*]~[#7]~[#6](~[Ru])~[#6,#7])]
[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8](~[Ru])~[#6]
[#6]([#6])([#6])~[#8]~[Ru]"
sprtGOAtmSpecKey="
[Ru] FFHINT= ~~210
[\$([#6;X3]([#1])(~[Ru])~[#6;X3](@~[#6;X3])@~[#6;X3])] FFHINT= ~~206
[\$([#6;X4]([#1])([#1])~[Ru])] FFHINT= ~~207
[\$([#7;X3](!@~[c])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#6;X4]),\$([#7;X3](!@~[c])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#1])] FFHINT= ~~208
[\$([#7;X3](!@~[c])~[#6](~[Ru])~[#6;X4](~[#6])(~[#6])~[#6]),\$([#7;X3]1~[#6;X3](~[Ru])~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2~[#6]~1)] FFHINT= ~~208
[\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6;X3](~[Ru])~[#7;X3]!@~[c]),\$([#7;X3](~[#6;X4])([#1])~[#6;X3](~[Ru])~[#7;X3]!@~[c])] FFHINT= ~~209
[\$([#7;X3](!@~[c])~[#6;X3](~[Ru])~[#7;X3]!@~[c]),\$([#7;X3](!@~[c])~[#6;X3](~[Ru])~[#7;X3]([#1])[#1]),\$([#7;X3]([#1])([#1])~[#6;X3](~[Ru])~[#7;X3]!@~[c])] FFHINT= ~~211
[\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#6;X4]),\$([#7;X3](!@~[#6;X4])(@~[#6;X3])@~[#6;X3](~[Ru])@~[#7;X3](@~[#6;X3])!@~[#6;X4]),\$([#7;X3](~[#6;X4])([#1])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#6;X4]),\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6;X3](~[Ru])~[#7;X3]([#1])~[#6;X4]),\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6;X3](~[Ru])~[#7;X3]([#1])[#1]),\$([#7;X3]([#1])([#1])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#6;X4]),\$([#7;X3](~[#6;X4])([#1])~[#6;X3](~[Ru])~[#7;X3]([#1])~[#6;X4]),\$([#7;X3](~[#6;X4])([#1])~[#6;X3](~[Ru])~[#7;X3]([#1])[#1]),\$([#7;X3]([#1])([#1])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#1]),\$([#7;X3]([#1])([#1])~[#6;X3](~[Ru])~[#7;X3]([#1])[#1]),\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6](~[Ru])~[#6;X4](~[#6])(~[#6])~[#6]),\$([#7;X3](~[#6;X4])(~[#1])~[#6](~[Ru])~[#6;X4](~[#6])(~[#6])~[#6])] FFHINT= ~~212"
sprtConstrGO="
[Ru] [\$([#6](~[Ru])@~[#7])] [\$([#7]@~[#6]~[Ru])] [\$([#6]!@~[#7]@~[#6]~[Ru])]"

# Spartan semi-empirical refinement
sprtGOSEFrozen="
[Ru]
[#9,#17,#35,#53,#7,#8,#15,#16,#1]~[Ru]
[\$([#6]~[Ru]);!\$([#6](~[Ru])~[#7])]
[\$([#6,#1]~[#6]~[Ru]);!\$([#6,#1]~[#6](~[Ru])~[#7])]
[\$([*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[#6](~[Ru])~[#7]);!\$([*]~[#7]~[#6](~[Ru])~[#6,#7])]
[\$([*]~[*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[*]~[#6]~[#7]);!\$([*]~[*]~[#7]~[#6](~[Ru])~[#6,#7])]
[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8](~[Ru])~[#6]
[#6]([#6])([#6])~[#8]~[Ru]"

# Spartan generation of acyclic Carbene-N-R conformers
skipDynCon=1 # if 0 then the search for conformers Carbene-N-R is skipped
sprtDYNConstraints="[\$([#6](~[#1])~[Ru])] [Ru] [\$([#6](~[#7])~[Ru]),\$([#6]1(~[#44])~[#6;X3]~[#6;X3]~[#7;X3]~[#6;X3]~[#6;X3]~1)] [\$([#7]~[#6]~[Ru]),\$([#6;X3]1~[#6;X3]~[#7;X3]~[#6;X3]~[#6;X3]~[#6;X3](~[#44])~1)]  0.0 0.0 180.0 1
[\$([#6]~[#15])] [#15] [Ru] [\$([#6](~[#1])~[Ru])] 0.0 0.0 180 1"
# ^Inclusion of PCy3 and PPh3

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
    head -n "$nn" "$file" > "${file}_tmp"
    echo "> <$propName>" >> "${file}_tmp"
    echo "$propValue" >> "${file}_tmp"
    echo "" >> "${file}_tmp"
    echo '$$$$' >> "${file}_tmp"
    mv "${file}_tmp" "$file"
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
    fileM="$3"
    awk -v f=$fileM '{print > out}; /\$\$\$\$/{out=f"_tmpAPrp"nm++".sdf"}' nm=1 out="${fileM}_tmpAPrp0.sdf" "$fileM"
    for tmpFile in $(ls "${fileM}_tmpAPrp"*".sdf" )
    do
        addPropertyToSingleMolSDF "$propNameM" "$propValueM" "$tmpFile"
        cat "$tmpFile" >> "${fileM}_allTmpAPrp"
        rm -f "$tmpFile"
    done
    mv "${fileM}_allTmpAPrp" "$fileM"
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
    echo "Finishing job!"
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
    errMsg="#3DBuild: wrong number of arguments"
    tcl="${tcl}PID$$"
    exit "$E_FATAL"
fi
inpSDF="$1"
outSDF="$2"
wrkDir="$3"
jobId="$4"
tcl="$wrkDir/$tcl$jobId"

molName=`basename "$inpSDF" .sdf`
molNum=`basename "$inpSDF" "$inputLabelExt"`

#
# Setup Log file
#
log="$jobId.log"
exec > $log
exec 2>&1
echo "Log for 3D-model Build and Refine" 
echo "Input: $inpSDF" 
echo "WorkDir: $wrkDir"
echo "Output: $outSDF"
echo "JobID: $jobId"
echo "TCL: $tcl"
echo "malName: $molName"
echo "molNum: $molNum"

#
# Wait for available token 
#
for i in $(seq 1 $lckMaxWait)
do
    for j in $(seq 1 $MAXPARALLELPROCESSES)
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
        errMsg="#LockConfSearch: too log wait for lockfile"
        abandon "$inpSDF" "$E_OPTERROR"
    else
        echo "All $MAXPARALLELPROCESSES lock files already taken:  waiting (i: $i)"
    fi
    sleep "$lockStep$lockTimeUnit"
done


#
# Parse jobID
#
jobNum=0
jobTyp=0
stereo=0
if [[ "$jobId" == *.* ]] 
then
    pA="${jobId%.*}"
    pB="${jobId#*.}"
    if [ "$pA" == "" ] || [ "$pB" == "" ]
    then
        errMsg="#3DBuild: error parsing jobId"
        abandon "$inpSDF" "$E_FATAL"
    else
        jobNum="$pA"
        if [[ "${#pB}" == 1 ]] 
        then
            jobTyp="$pB"
        else
            jobTyp="${pB:0:1}"
            stereo="${pB:1:2}"
        fi
    fi
else
    errMsg="#3DBuild: error parsing jobId"
    abandon "$inpSDF" "$E_FATAL"
fi
echo "JobNumber:      $jobNum"
echo "JobType:        $jobTyp"
echo "Stereochemistry: $stereo"


#
# Setup the job type-dependent parameters
#
graphEditTasks=""
applyCartesianMove=1
libFrags="$wrkDir/${molNum}_libFragsAndLigand.sdf" #NB: defined also in other scripts!
libCaps="$wrkDir/${molNum}_libFragsAndLigand.sdf" #NB: defined also in other scripts!
echo "Using fragment library $libFrags"
case "$jobTyp$stereo" in
    "A1") 
        libScaff="$libScaffA1" 
        ;;
    "A2")
        libScaff="$libScaffA2"
        ;;
    "F1") 
        libScaff="$libScaffA1" 
        ;;
    "F2")
        libScaff="$libScaffA2"
        ;;
    "B1") 
        libScaff="$libScaffB1"
        ;;
    "B2")
        libScaff="$libScaffB2"
        ;;
    "C1") 
        libScaff="$libScaffC1"
        skipDynCon=0
        ;;
    "C2")
        libScaff="$libScaffC2"
        skipDynCon=0
        ;;
    "D1") 
        libScaff="$libScaffD1"
        skipDynCon=0
        ;;
    "D2")
        libScaff="$libScaffD2"
        skipDynCon=0
        ;;
    "E1")
        editGraph=0
        graphEditTasks="$grpEdtRemoveXLigands"
        libScaff="$libScaffE1"
        ;;
    "E2")
        editGraph=0
        graphEditTasks="$grpEdtRemoveXLigands"
        libScaff="$libScaffE2"
        ;;
    "X1") 
        libScaff="$libScaffX1"
        ;;
    "X2")
        libScaff="$libScaffX2"
        ;;
    "Y1") 
        libScaff="$libScaffY1"
        applyCartesianMove=0
        crtMove="$crtMoveY1"
        scalingMove="$scalingMoveY1"
        crtMvRef=$crtMvRefY1
        ;;
    "Y2")
        libScaff="$libScaffY2"
        applyCartesianMove=0
        crtMove="$crtMoveY2"
        scalingMove="$scalingMoveY2"
        crtMvRef=$crtMvRefY2
        ;;
    "Z1") 
        libScaff="$libScaffZ1"
        ;;
    "Z2")
        libScaff="$libScaffZ2"
        ;;
    *) 
        errMsg="#3DBuild: error interpreting jobId $jobTyp$stereo"
        abandon "$inpSDF" "$E_FATAL"
        ;;
esac
echo "Using scaffolds library $libScaff"
molName="${molNum}_$jobTyp${stereo}"
locInpSDF="$molName$inputLabelExt"
locInpSDFpre="${molName}_pre$inputLabelExt"
echo "Molecule name set to $molName"


#
# Edit graph upon request
#
if [[ "$editGraph" == 0 ]]
then
    GrEdPar="$wrkDir/${molName}_GrEd.par"
    GrEdLog="$wrkDir/${molName}_GrEd.log"
    echo "GRAPHEDIT-INPUTGRAPHS=$inpSDF" > "$GrEdPar"
    echo "GRAPHEDIT-GRAPHSEDITSFILE=$graphEditTasks" >> "$GrEdPar"
    echo "GRAPHEDIT-OUTPUTGRAPHS=$locInpSDFpre" >> "$GrEdPar"
    echo "GRAPHEDIT-OUTPUTGRAPHSFORMAT=SDF" >> "$GrEdPar"
    echo "GRAPHEDIT-ENFORCESYMMETRY=yes" >> "$GrEdPar"
    echo "GRAPHEDIT-VERBOSITY=0" >> "$GrEdPar"
    echo "FS-scaffoldLibFile=$libScaff" >> "$GrEdPar"
    echo "FS-fragmentLibFile=$libFrags" >> "$GrEdPar"
    echo "FS-cappingFragmentLibFile=$libCaps" >> "$GrEdPar"
    echo "FS-compMatrixFile=$compMatrix" >> "$GrEdPar"
    # Launch GraphEditor
    denoptim -r GE "$GrEdPar" > "$GrEdLog" 2>&1 
    # Check outcome
    if [ ! -f "$GrEdLog" ]; then
        echo "$GrEdLog not found."
        errMsg="#GraphEditor: log file $GrEdLog not found."
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    if [ ! -f "$locInpSDFpre" ]; then
        echo "$locInpSDF not found."
        errMsg="#GraphEditor: output file $locInpSDFpre not found."
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    if ! grep -q "Completed GraphEditor" "$GrEdLog" ; then
        errMsg="#GraphEditor: non-zero exit status"
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    if [ "$cleanup" == 0 ]
    then
        rm -f "$GrEdPar"
        rm -f "$GrEdLog"
    fi
else
    cp "$inpSDF" "$locInpSDFpre"
fi
 
# Creating finalized input sdf
GrEdPar2="$wrkDir/${molName}_GrEd2.par"
GrEdLog2="$wrkDir/${molName}_GrEd2.log"
graphEditTasks="$grpEdtReplaceLigands"
echo "GRAPHEDIT-INPUTGRAPHS=$locInpSDFpre" > "$GrEdPar2"
echo "GRAPHEDIT-GRAPHSEDITSFILE=$graphEditTasks" >> "$GrEdPar2"
echo "GRAPHEDIT-OUTPUTGRAPHS=$locInpSDF" >> "$GrEdPar2"
echo "GRAPHEDIT-OUTPUTGRAPHSFORMAT=SDF" >> "$GrEdPar2"
echo "GRAPHEDIT-ENFORCESYMMETRY=yes" >> "$GrEdPar2"
echo "GRAPHEDIT-VERBOSITY=0" >> "$GrEdPar2"
echo "FS-scaffoldLibFile=$libScaff" >> "$GrEdPar2"
echo "FS-fragmentLibFile=$libFrags" >> "$GrEdPar2"
# Launch Graph Editor
denoptim -r GE "$GrEdPar2" > "$GrEdLog2" 2>&1 

# Check outcome
if [ ! -f "$GrEdLog2" ]; then
    echo "$GrEdLog2 not found."
    errMsg="#GraphEditor: log file $GrEdLog2 not found."
    abandon "$inpSDF" "$E_OPTERROR"
fi
if [ ! -f "$locInpSDF" ]; then
    echo "$locInpSDF not found."
    errMsg="#GraphEditor: output file $locInpSDF not found."
    abandon "$inpSDF" "$E_OPTERROR"
fi
if ! grep -q "Completed GraphEditor" "$GrEdLog2" ; then
    errMsg="#GraphEditor: non-zero exit status"
    abandon "$inpSDF" "$E_OPTERROR"
fi
if [ "$cleanup" == 0 ]
then
    rm -f "$GrEdPar2"
    rm -f "$GrEdLog2"
    rm -f "$locInpSDFpre"
fi

#
# Build 3D model
#
echo "Starting DenoptimCG..."
# Setting new params
DnCG3Dout="$wrkDir/${molName}_3D-CG.sdf"
DnCGParFile="$wrkDir/${molName}_DnCG.par"
# prepare input
portablesed "1s/.*/$molName/" "$locInpSDF"
# prepare param file
echo "3DB-VERBOSITY=1" > "$DnCGParFile"
echo "3DB-inpSDF=$locInpSDF" >> "$DnCGParFile"
echo "3DB-outSDF=$DnCG3Dout" >> "$DnCGParFile"
echo "FS-scaffoldLibFile=$libScaff" >> "$DnCGParFile"
echo "FS-fragmentLibFile=$libFrags" >> "$DnCGParFile"
#echo "FS-cappingFragmentLibFile=$libCaps" >> "$DnCGParFile"
echo "FS-compMatrixFile=$compMatrix" >> "$DnCGParFile"
echo "FS-ROTBONDSDEFFILE=$rotSpaceDef" >> "$DnCGParFile"
echo "3DB-atomOrderingScheme=1" >> "$DnCGParFile"
echo "3DB-workDir=$wrkDir" >> "$DnCGParFile"
echo "3DB-TOOLPSSROT=$pathToTinkerBin/pssrot" >> "$DnCGParFile"
echo "3DB-TOOLXYZINT=$pathToTinkerBin/xyzint" >> "$DnCGParFile"
echo "3DB-TOOLINTXYZ=$pathToTinkerBin/intxyz" >> "$DnCGParFile"
echo "3DB-FORCEFIELDFILE=$tinkerForceField" >> "$DnCGParFile"
echo "3DB-KEYFILE=$tinkerKeyFile" >> "$DnCGParFile"
echo "3DB-PSSROTPARAMS=$tinkerSubmitFile" >> "$DnCGParFile"
# launch DenoptimCG
denoptim -r B3D "$DnCGParFile" 
# Cleanup files tmp files
echo "Removing $wrkDir/$molNum"_cs0."*"
rm -f "$wrkDir/${molNum}_cs0."*
# Adding molName to output sdf
portablesed "1s/.*/$molName/" "$DnCG3Dout"
# Check outcome
if [ -f "$DnCG3Dout" ]; then
    if grep -q "MOL_ERROR" "$DnCG3Dout" ; then
        echo "Exiting! Found ERROR in output of DenoptimCG"
        cat "$DnCG3Dout" > "$outSDF"
        writeTCL "$jobId" 
        releaseLock "$lockFile" "$lockIDFile"
        exit "$E_OPTERROR"
    fi
else
    echo "$DnCG3Dout not found."
    errMsg="#DenoptimCG: $DnCG3Dout not found."
    abandon "$inpSDF" "$E_OPTERROR"
fi


#
# Check for exceeding atom crowding
#
if [ "$checkAtomClashes" -eq 0 ]
then
    echo "Starting Atom Clash detector..."
    # Setting new params
    atmClshParFile="$wrkDir/${molName}_AtmClsh.par"
    atmClshLog="$wrkDir/${molName}_AtmClsh.log"
    # Prepare param file
    echo "VERBOSITY: 1" > "$atmClshParFile"
    echo "TASK: AnalyzeVDWClashes" >> "$atmClshParFile"
    echo "INFILE: $DnCG3Dout" >> "$atmClshParFile"
    echo "ALLOWANCE13: 1.0" >> "$atmClshParFile"
    echo "CUTOFF: 0.75" >> "$atmClshParFile"
    echo "ALLOWANCE: 0.40" >> "$atmClshParFile"
    # Launch ACC
    autocompchem -p "$atmClshParFile" > "$atmClshLog"
    # Check outcome
    if [ ! -f "$atmClshLog" ]; then
        errMsg="#AtomClash: $atmClshLog not found."
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    if ! grep -q "Termination status: 0" "$atmClshLog" ; then
        errMsg="#AtomClash: non-zero exit status from AutoCompChem"
        abandon "$DnCG3Dout" "$E_OPTERROR" 
    fi
    if grep -q 'Found 0 mols with one or more atom clashes' "$atmClshLog" ; then
        echo "No atom clashes"
    else
        errMsg="#AtomClash: Found Atom Clashes"
        abandon "$DnCG3Dout" "$E_OPTERROR"
    fi
    if [ "$cleanup" == 0 ]
    then
        rm -f "$atmClshParFile"
        rm -f "$atmClshLog"
    fi
fi

#
# Prepare input for Ru-NHC conformation-constrained optimization
#
echo "Prepare input for Ru-NHC conformation-contrained geometry optimization..."
# Setting new params
accSprGOParFile="$wrkDir/${molName}_SprtGOIn.par"
accSprGOLog="$wrkDir/${molName}_SprtGOIn.log"
accSprGOInp="$wrkDir/${molName}_SprtGOIn.sdf"
sprtGOInp="$wrkDir/${molName}_GO.spartan"
cp "$DnCG3Dout" "$accSprGOInp"
portablesed "s/${molName}.*/${molName}_GO/g" "$accSprGOInp"
# Prepare param file
echo "VERBOSITY: 1" > "$accSprGOParFile"
echo "TASK: PrepareInputSpartan" >> "$accSprGOParFile"
echo "INFILE: $accSprGOInp" >> "$accSprGOParFile"
echo "KEYWORDS: OPT MMFF PRINTLEV=3 OPTCYCLE=26000 USEFFPARMS=~/PARAMS.MMFF.DENOPTIM" >> "$accSprGOParFile"
echo "CHARGE: $charge" >> "$accSprGOParFile"
echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprGOParFile"
echo "OUTFILE: $sprtGOInp" >> "$accSprGOParFile"
echo "\$STARTATOMSPECIFICKEYWORDS: $sprtGOAtmSpecKey" >> "$accSprGOParFile"
echo "\$END" >> "$accSprGOParFile"
echo "\$STARTFREEZEATOMS: $sprtGOFrozen" >> "$accSprGOParFile"
echo "\$END" >> "$accSprGOParFile"
# Launch ACC
autocompchem -p "$accSprGOParFile" > "$accSprGOLog"
# Check outcome
if [ ! -f "$accSprGOLog" ]; then
    errMsg="#MakeSpartanGOInput: $accSprGOLog not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -d "$sprtGOInp" ]; then
    errMsg="#MakeSpartanGOInput: $sprtGOInp directory not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Termination status: 0" "$accSprGOLog" ; then
    errMsg="#MakeSpartanGOInput: non-zero exit status from AutoCompChem."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi


#
# Submit Ru-NHC conformation-constrained geometry optimization
#
echo "Submitting Spartan Ru-NHC conformation-constrained optimization..."
# Setting new parameters
sprtGOOutLog="$sprtGOInp/${molName}_GO/output"
sprtGOLog="$wrkDir/${molName}_SprtGORun.log"
# Lauch Spartan
"$SPARTANEXE" --foreground-submit "$sprtGOInp" > "$sprtGOLog"
# Check outcome
if [ ! "$?" == 0 ]; then
    errMsg="#RunSpartanGO: non-zero exit status from Spartan call."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -f "$sprtGOOutLog" ]; then
    errMsg="#RunSpartanGO: $sprtGOOutLog file not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Successful completion" "$sprtGOOutLog" ; then
    errMsg="#RunSpartanGO: unsuccessful completion."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi

#
# Extract optimized geometry
#
echo "Extracting results from Spartan constrained geometry optimization..."
# Setting new params
accSprGOOutParFile="$wrkDir/${molName}_SprtGOOut.par"
accSprGOOutLog="$wrkDir/${molName}_SprtGOOut.log"
sprtGOOut="$wrkDir/${molName}_SprtGOOut.sdf"
# Prepare param file
echo "VERBOSITY: 1" > "$accSprGOOutParFile"
echo "TASK: ExtractLastGeometryFromSpartanTree" >> "$accSprGOOutParFile"
echo "INFILE: $sprtGOInp" >> "$accSprGOOutParFile"
echo "OUTFORMAT: SDF" >> "$accSprGOOutParFile"
echo "OUTFILE: $sprtGOOut" >> "$accSprGOOutParFile"
# Launch ACC
autocompchem -p "$accSprGOOutParFile" > "$accSprGOOutLog"
# Check outcome
if [ ! -f "$accSprGOOutLog" ]; then
    errMsg="#ExtractSpartanGOOutput: $accSprGOOutLog not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -f "$sprtGOOut" ]; then
    errMsg="#ExtractSpartanGOOutput: $sprtGOOut file not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Termination status: 0" "$accSprGOOutLog" ; then
    errMsg="#ExtractSpartanGOOutput: non-zero exit status from AutoCompChem."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
portablesed "1 s/^.*$/$molName/g" "$sprtGOOut"

    
#
# Prepare input semi-empirical refinement
#
echo "Prepare input for semi-empirical geometry optimization ..."
# Setting new params
accSprGOSEParFile="$wrkDir/${molName}_SprtGOSE0In.par"
accSprGOSELog="$wrkDir/${molName}_SprtGOSE0In.log"
accSprGOSEInp="$wrkDir/${molName}_SprtGOSE0In.sdf"
sprtGOSEInp="$wrkDir/${molName}_GOSE0.spartan"
cp "$sprtGOOut" "$accSprGOSEInp"
portablesed "s/${molName}.*/${molName}_GOSE0/g" "$accSprGOSEInp"
# Prepare param file
echo "VERBOSITY: 1" > "$accSprGOSEParFile"
echo "TASK: PrepareInputSpartan" >> "$accSprGOSEParFile"
echo "INFILE: $accSprGOSEInp" >> "$accSprGOSEParFile"
echo "KEYWORDS: OPT PM6 PRINTLEV=1 OPTCYCLE=26000" >> "$accSprGOSEParFile"
echo "CHARGE: $charge" >> "$accSprGOSEParFile"
echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprGOSEParFile"
echo "OUTFILE: $sprtGOSEInp" >> "$accSprGOSEParFile"
echo "\$STARTFREEZEATOMS: $sprtGOSEFrozen" >> "$accSprGOSEParFile"
echo "\$END" >> "$accSprGOSEParFile"
# Launch ACC
autocompchem -p "$accSprGOSEParFile" > "$accSprGOSELog"
# Check outcome
if [ ! -f "$accSprGOSELog" ]; then
    errMsg="#MakeSpartanGOSEInput: $accSprGOSELog not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -d "$sprtGOSEInp" ]; then
    errMsg="#MakeSpartanGOSEInput: $sprtGOSEInp directory not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Termination status: 0" "$accSprGOSELog" ; then
    errMsg="#MakeSpartanGOSEInput: non-zero exit status from AutoCompChem."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi


#
# Submit semi-empirical refinement
#
echo "Submitting Spartan semi-empirical optimization..."
# Setting new parameters
sprtGOSEOutLog="$sprtGOSEInp/${molName}_GOSE0/output"
sprtGOSELog="$wrkDir/${molName}_SprtGOSE0Run.log"
# Lauch Spartan
"$SPARTANEXE" --foreground-submit "$sprtGOSEInp" > "$sprtGOSELog"
# Check outcome
if [ ! "$?" == 0 ]; then
    errMsg="#RunSpartanGOSE: non-zero exit status from Spartan call."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -f "$sprtGOSEOutLog" ]; then
    errMsg="#RunSpartanGOSE: $sprtGOSEOutLog file not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Successful completion" "$sprtGOSEOutLog" ; then
    errMsg="#RunSpartanGOSE: unsuccessful completion."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi


#
# Extract optimized geometry
#
echo "Extracting results from Spartan semi-empirical geometry optimization..."
# Setting new params
accSprGOSEOutParFile="$wrkDir/${molName}_SprtGOSE0Out.par"
accSprGOSEOutLog="$wrkDir/${molName}_SprtGOSE0Out.log"
sprtGOSEOut="$wrkDir/${molName}_SprtGOSE0Out.sdf"
# Prepare param file
echo "VERBOSITY: 1" > "$accSprGOSEOutParFile"
echo "TASK: ExtractLastGeometryFromSpartanTree" >> "$accSprGOSEOutParFile"
echo "INFILE: $sprtGOSEInp" >> "$accSprGOSEOutParFile"
echo "OUTFORMAT: SDF" >> "$accSprGOSEOutParFile"
echo "OUTFILE: $sprtGOSEOut" >> "$accSprGOSEOutParFile"
# Launch ACC
autocompchem -p "$accSprGOSEOutParFile" > "$accSprGOSEOutLog"
# Check outcome
if [ ! -f "$accSprGOSEOutLog" ]; then
    errMsg="#ExtractSpartanGOSEOutput: $accSprGOSEOutLog not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -f "$sprtGOSEOut" ]; then
    errMsg="#ExtractSpartanGOSEOutput: $sprtGOSEOut file not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Termination status: 0" "$accSprGOSEOutLog" ; then
    errMsg="#ExtractSpartanGOSEOutput: non-zero exit status from AutoCompChem."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
portablesed "1 s/^.*$/$molName/g" "$sprtGOSEOut"


#
# Extract energy
#
echo "Extracting energy of PM6-optimized geometry..."
if grep -q " SCF did not converge - abort optimization" "$sprtGOSEOutLog"
then
    pm6GOEnergy=1000000000
    echo "$n-PM6 Energy not found (FAILED SCF)"
else
    if grep -q " Heat of Formation: " "$sprtGOSEOutLog"
    then
        pm6GOEnergy=$(grep " Heat of Formation: " "$sprtGOSEOutLog" | tail -n 1 | awk '{print $4}')
        pm6GOEnergy=$(echo ${pm6GOEnergy} | sed -e 's/[eE]+*/\*10\^/')
        pm6GOEnergy=$(bc -l <<< "($pm6GOEnergy)")
        echo "PM6 Energy of refined conformation: $pm6GOEnergy (from GOSE0)"
    fi
fi
lowestEnergy="$pm6GOEnergy"
lowestEnergyConfSDF="$sprtGOSEOut"


#
# Evaluate alternative conformations along specific bonds
#
sprtDYNOutDir="$wrkDir/${molName}_DYN.Prof.${molName}_DYN.spartan"
if [ 1 == "$skipDynCon" ]
then
    #
    # Prepare input for a possible generation of conformers
    #
    echo "Checking possibility for generating specific conformers..."
    # Setting new params
    accSprDYNParFile="$wrkDir/${molName}_SprtDYNIn.par"
    accSprDYNLog="$wrkDir/${molName}_SprtDYNIn.log"
    sprtDYNInp="$wrkDir/${molName}_DYN.spartan"
    accSprDYNInp="$wrkDir/${molName}_SprtDYNIn.sdf"
    cp "$lowestEnergyConfSDF" "$accSprDYNInp"
    portablesed "s/${molName}.*/${molName}_DYN/g" "$accSprDYNInp"
    # Prepare param file
    echo "VERBOSITY: 1" > "$accSprDYNParFile"
    echo "TASK: PrepareInputSpartan" >> "$accSprDYNParFile"
    echo "INFILE: $accSprDYNInp" >> "$accSprDYNParFile"
    echo "KEYWORDS: DYNCON MMFF PRINTLEV=2 DYNCONMETHOD=GRID OPTCYCLE=2600 USEFFPARMS=~/PARAMS.MMFF.DENOPTIM" >> "$accSprDYNParFile"
    #
    echo "CHARGE: $charge" >> "$accSprDYNParFile"
    echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprDYNParFile"
    echo "OUTFILE: $sprtDYNInp" >> "$accSprDYNParFile"
    echo "\$STARTATOMSPECIFICKEYWORDS: $sprtGOAtmSpecKey" >> "$accSprDYNParFile"
    echo "\$END" >> "$accSprDYNParFile"
    echo "\$STARTFREEZEATOMS: $sprtGOSEFrozen" >> "$accSprDYNParFile"
    echo "\$END" >> "$accSprDYNParFile"
    echo "\$STARTDYNAMICCONSTRAINT: $sprtDYNConstraints" >> "$accSprDYNParFile"
    echo "\$END" >> "$accSprDYNParFile"
    # Launch ACC
    autocompchem -p "$accSprDYNParFile" > "$accSprDYNLog"
    # Check outcome
    if [ ! -f "$accSprDYNLog" ]; then
        errMsg="#MakeSpartanInput: $accSprDYNLog not found."
        abandon "$DnCG3Dout" "$E_OPTERROR"
    fi
    if [ ! -d "$sprtDYNInp" ]; then
        errMsg="#MakeSpartanInput: $sprtDYNInp directory not found."
        abandon "$DnCG3Dout" "$E_OPTERROR"
    fi
    if ! grep -q "Termination status: 0" "$accSprDYNLog" ; then
        errMsg="#MakeSpartanInput: non-zero exit status from AutoCompChem."
        abandon "$DnCG3Dout" "$E_OPTERROR"
    fi
    # Need to append the preference or Spartan will crash
    echo "BEGINPREFERENCES" >> "$sprtDYNInp/${molName}_DYN/input"
    echo " MM:CONF_SELECTION_RULE=2" >> "$sprtDYNInp/${molName}_DYN/input"
    echo "ENDPREFERENCES" >> "$sprtDYNInp/${molName}_DYN/input"
    
    
    #
    # Submit DYNamic constraint
    #
    echo "Submitting DYNCON generation of alternative conformation..."
    # Setting new parameters
    sprtDYNOutLog="$sprtDYNInp/${molName}_DYN/output"
    # thid one is defined outside if block
    #sprtDYNOutDir="$wrkDir/${molName}_DYN.Prof.${molName}_DYN.spartan"
    sprtDYNLog="$wrkDir/${molName}_SprtDYNRun.log"
    # Lauch Spartan
    "$SPARTANEXE" --foreground-submit "$sprtDYNInp" > "$sprtDYNLog"
    # Check outcome
    if [ ! "$?" == 0 ]; then
        errMsg="#RunSpartan: non-zero exit status from Spartan call."
        abandon "$DnCG3Dout" "$E_OPTERROR"
    fi
    if [ ! -f "$sprtDYNOutLog" ]; then
        errMsg="#RunSpartan: $sprtDYNOutLog file not found."
        abandon "$DnCG3Dout" "$E_OPTERROR"
    fi
    #NB: Not mandatory to complete this step successfully
    #if [ ! -d "$sprtDYNOutDir" ]; then
    #    errMsg="#RunSpartan: $sprtDYNOutDir folder not found."
    #    abandon "$DnCG3Dout" "$E_OPTERROR"
    #fi
    if ! grep -q "Successful completion" "$sprtDYNOutLog" ; then
        errMsg="#RunSpartan: unsuccessful completion."
        abandon "$DnCG3Dout" "$E_OPTERROR"
    fi
fi
    
    
#
# Refine the alternative conformers
#
n=0
if [ -d "$sprtDYNOutDir" ]
then
    lowestEnergyConfSDF="$lowestEnergyConfSDF"
    lowestEnergy="$lowestEnergy"
    sprtDYNConfs=($(ls -d "$sprtDYNOutDir/Profile."*))
    for conformer in ${sprtDYNConfs[@]}
    do
        n=$((n+1))
        #
        # Extracting structure of conformer
        #
        echo "$n-Extracting input geometry of conformer..."
        # Setting new params
        accSprDYNOutParFile="$wrkDir/${molName}_SprtDYNOut${n}.par"
        accSprDYNOutLog="$wrkDir/${molName}_SprtDYNOut${n}.log"
        sprtDYNOut="$wrkDir/${molName}_SprtDYNOut${n}.sdf"
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprDYNOutParFile"
        echo "TASK: ExtractLastGeometryFromSpartanTree" >> "$accSprDYNOutParFile"
        echo "INFILE: $sprtDYNOutDir" >> "$accSprDYNOutParFile"
        echo "MOLNAME: Profile.${n}" >> "$accSprDYNOutParFile"
        echo "OUTFORMAT: SDF" >> "$accSprDYNOutParFile"
        echo "OUTFILE: $sprtDYNOut" >> "$accSprDYNOutParFile"
        # Launch ACC
        autocompchem -p "$accSprDYNOutParFile" > "$accSprDYNOutLog"
        # Check outcome
        if [ ! -f "$accSprDYNOutLog" ]; then
            errMsg="#ExtractSpartanOutput: $accSprDYNOutLog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtDYNOut" ]; then
            errMsg="#ExtractSpartanOutput: $sprtDYNOut file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprDYNOutLog" ; then
            errMsg="#ExtractSpartanOutput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        portablesed "s/Profile.*/$molName/g" "$sprtDYNOut"


        #
        # Prepare input molecular mechanics refinement
        #
        echo "$n-Prepare input for molecular mechanics geometry optimization..."
        # Setting new params
        accSprGOMMParFile="$wrkDir/${molName}_SprtGOMM${n}In.par"
        accSprGOMMLog="$wrkDir/${molName}_SprtGOMM${n}In.log"
        accSprGOMMInp="$wrkDir/${molName}_SprtGOMM${n}In.sdf"
        sprtGOMMInp="$wrkDir/${molName}_GOMM${n}.spartan"
        cp "$sprtDYNOut" "$accSprGOMMInp"
        portablesed "s/${molName}.*/${molName}_GOMM${n}/g" "$accSprGOMMInp"
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprGOMMParFile"
        echo "TASK: PrepareInputSpartan" >> "$accSprGOMMParFile"
        echo "INFILE: $accSprGOMMInp" >> "$accSprGOMMParFile"
        echo "KEYWORDS: OPT MMFF PRINTLEV=2  OPTCYCLE=26000 USEFFPARMS=~/PARAMS.MMFF.DENOPTIM" >> "$accSprGOMMParFile"
        echo "CHARGE: $charge" >> "$accSprGOMMParFile"
        echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprGOMMParFile"
        echo "OUTFILE: $sprtGOMMInp" >> "$accSprGOMMParFile"
        echo "\$STARTATOMSPECIFICKEYWORDS: $sprtGOAtmSpecKey" >> "$accSprGOMMParFile"
        echo "\$END" >> "$accSprGOMMParFile"
        echo "\$STARTFREEZEATOMS: $sprtGOFrozen" >> "$accSprGOMMParFile"
        echo "\$END" >> "$accSprGOMMParFile"
        # Launch ACC
        autocompchem -p "$accSprGOMMParFile" > "$accSprGOMMLog"
        # Check outcome
        if [ ! -f "$accSprGOMMLog" ]; then
            errMsg="#MakeSpartanGOMMInput: $accSprGOMMLog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -d "$sprtGOMMInp" ]; then
            errMsg="#MakeSpartanGOMMInput: $sprtGOMMInp directory not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprGOMMLog" ; then
            errMsg="#MakeSpartanGOMMInput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi


        #
        # Submit molecular mechanics refinement
        #
        echo "$n-Submitting Spartan molecular mechanics optimization..."
        # Setting new parameters
        sprtGOMMOutLog="$sprtGOMMInp/${molName}_GOMM${n}/output"
        sprtGOMMLog="$wrkDir/${molName}_SprtGOMM${n}Run.log"
        # Lauch Spartan
        "$SPARTANEXE" --foreground-submit "$sprtGOMMInp" > "$sprtGOMMLog"
        # Check outcome
        if [ ! "$?" == 0 ]; then
            errMsg="#RunSpartanGOMM: non-zero exit status from Spartan call."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtGOMMOutLog" ]; then
            errMsg="#RunSpartanGOMM: $sprtGOMMOutLog file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Successful completion" "$sprtGOMMOutLog" ; then
            errMsg="#RunSpartanGOMM: unsuccessful completion."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi

        #
        # Extract optimized geometry
        #
        echo "$n-Extracting results from Spartan molecular mechanics geometry optimization..."
        # Setting new params
        accSprGOMMOutParFile="$wrkDir/${molName}_SprtGOMM${n}Out.par"
        accSprGOMMOutLog="$wrkDir/${molName}_SprtGOMM${n}Out.log"
        sprtGOMMOut="$wrkDir/${molName}_SprtGOMM${n}Out.sdf"
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprGOMMOutParFile"
        echo "TASK: ExtractLastGeometryFromSpartanTree" >> "$accSprGOMMOutParFile"
        echo "INFILE: $sprtGOMMInp" >> "$accSprGOMMOutParFile"
        echo "OUTFORMAT: SDF" >> "$accSprGOMMOutParFile"
        echo "OUTFILE: $sprtGOMMOut" >> "$accSprGOMMOutParFile"
        # Launch ACC
        autocompchem -p "$accSprGOMMOutParFile" > "$accSprGOMMOutLog"
        # Check outcome
        if [ ! -f "$accSprGOMMOutLog" ]; then
            errMsg="#ExtractSpartanGOMMOutput: $accSprGOMMOutLog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtGOMMOut" ]; then
            errMsg="#ExtractSpartanGOMMOutput: $sprtGOMMOut file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprGOMMOutLog" ; then
            errMsg="#ExtractSpartanGOMMOutput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        portablesed "1 s/^.*$/$molName/g" "$sprtGOMMOut"

 
        #
        # Prepare input semi-empirical refinement
        #
        echo "$n-Prepare input for semi-empirical geometry optimization..."
        # Setting new params
        accSprGOSEParFile="$wrkDir/${molName}_SprtGOSE${n}In.par"
        accSprGOSELog="$wrkDir/${molName}_SprtGOSE${n}In.log"
        accSprGOSEInp="$wrkDir/${molName}_SprtGOSE${n}In.sdf"
        sprtGOSEInp="$wrkDir/${molName}_GOSE${n}.spartan"
        cp "$sprtGOMMOut" "$accSprGOSEInp"
        portablesed "s/${molName}.*/${molName}_GOSE${n}/g" "$accSprGOSEInp"
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprGOSEParFile"
        echo "TASK: PrepareInputSpartan" >> "$accSprGOSEParFile"
        echo "INFILE: $accSprGOSEInp" >> "$accSprGOSEParFile"
        echo "KEYWORDS: OPT PM6 PRINTLEV=1 OPTCYCLE=2600" >> "$accSprGOSEParFile"
        echo "CHARGE: $charge" >> "$accSprGOSEParFile"
        echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprGOSEParFile"
        echo "OUTFILE: $sprtGOSEInp" >> "$accSprGOSEParFile"
        echo "\$STARTFREEZEATOMS: $sprtGOSEFrozen" >> "$accSprGOSEParFile"
        echo "\$END" >> "$accSprGOSEParFile"
        # Launch ACC
        autocompchem -p "$accSprGOSEParFile" > "$accSprGOSELog"
        # Check outcome
        if [ ! -f "$accSprGOSELog" ]; then
            errMsg="#MakeSpartanGOSEInput: $accSprGOSELog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -d "$sprtGOSEInp" ]; then
            errMsg="#MakeSpartanGOSEInput: $sprtGOSEInp directory not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprGOSELog" ; then
            errMsg="#MakeSpartanGOSEInput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi


        #
        # Submit semi-empirical refinement
        #
        echo "$n-Submitting Spartan semi-empirical optimization..."
        # Setting new parameters
        sprtGOSEOutLog="$sprtGOSEInp/${molName}_GOSE${n}/output"
        sprtGOSELog="$wrkDir/${molName}_SprtGOSE${n}Run.log"
        # Lauch Spartan
        "$SPARTANEXE" --foreground-submit "$sprtGOSEInp" > "$sprtGOSELog"
        # Check outcome
        if [ ! "$?" == 0 ]; then
            errMsg="#RunSpartanGOSE: non-zero exit status from Spartan call."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtGOSEOutLog" ]; then
            errMsg="#RunSpartanGOSE: $sprtGOSEOutLog file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Successful completion" "$sprtGOSEOutLog" ; then
            errMsg="#RunSpartanGOSE: unsuccessful completion."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi


        #
        # Extract optimized geometry
        #
        echo "$n-Extracting results from Spartan semi-empirical geometry optimization..."
        # Setting new params
        accSprGOSEOutParFile="$wrkDir/${molName}_SprtGOSE${n}Out.par"
        accSprGOSEOutLog="$wrkDir/${molName}_SprtGOSE${n}Out.log"
        sprtGOSEOut="$wrkDir/${molName}_SprtGOSE${n}Out.sdf"
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprGOSEOutParFile"
        echo "TASK: ExtractLastGeometryFromSpartanTree" >> "$accSprGOSEOutParFile"
        echo "INFILE: $sprtGOSEInp" >> "$accSprGOSEOutParFile"
        echo "OUTFORMAT: SDF" >> "$accSprGOSEOutParFile"
        echo "OUTFILE: $sprtGOSEOut" >> "$accSprGOSEOutParFile"
        # Launch ACC
        autocompchem -p "$accSprGOSEOutParFile" > "$accSprGOSEOutLog"
        # Check outcome
        if [ ! -f "$accSprGOSEOutLog" ]; then
            errMsg="#ExtractSpartanGOSEOutput: $accSprGOSEOutLog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtGOSEOut" ]; then
            errMsg="#ExtractSpartanGOSEOutput: $sprtGOSEOut file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprGOSEOutLog" ; then
            errMsg="#ExtractSpartanGOSEOutput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        portablesed "1 s/^.*$/$molName/g" "$sprtGOSEOut"


        #
        # Extract energy
        #
        if grep -q " SCF did not converge - abort optimization" "$sprtGOSEOutLog"
        then
            pm6GOEnergy=1000000000
            echo "$n-PM6 Energy not found (FAILED SCF)"
        else 
            if grep -q " Heat of Formation: " "$sprtGOSEOutLog"
            then
            echo "$n-Extracting PM6 energy of optimized geometry..."
                pm6GOEnergy=$(grep " Heat of Formation: " "$sprtGOSEOutLog" | tail -n 1 | awk '{print $4}')
                pm6GOEnergy=$(echo ${pm6GOEnergy} | sed -e 's/[eE]+*/\*10\^/')
                pm6GOEnergy=$(bc -l <<< "($pm6GOEnergy)")
                echo "$n-PM6 Energy of refined conformation: $pm6GOEnergy (from GOSE$n)"
    
                #
                # Possibly update lowest value
                #
                newLowest=$( bc -l <<< "$pm6GOEnergy < $lowestEnergy")
                if [ "$newLowest" == 1 ]
                then
                    echo "   ==> New lowest PM6 energy conformation found! <=="
                    lowestEnergy="$pm6GOEnergy"
                    lowestEnergyConfSDF="$sprtGOSEOut"
                fi
            fi
        fi
    done
fi


#
# Edit geometry according to Cartesian move
#
if [ "$applyCartesianMove" == 0 ]
then
    #
    # Generate modified geometry
    #
    echo "Editing geometry according to cartesina move..."
    # Setting new parameters
    accGeomEdtParFile="$wrkDir/${molName}_geomEdt.par"
    accGeomEdtLogFile="$wrkDir/${molName}_geomEdt.log"
    accGeomEdtOutFile="$wrkDir/${molName}_geomEdt.sdf"
    echo "VERBOSITY: 1" > "$accGeomEdtParFile"
    echo "TASK: ModifyGeometry" >> "$accGeomEdtParFile"
    echo "INFILE: $lowestEnergyConfSDF" >> "$accGeomEdtParFile"
    echo "CARTESIANMOVE: $crtMove" >> "$accGeomEdtParFile"
    echo "CARTESIANSCALINGFACTORS: $scalingMove" >> "$accGeomEdtParFile"
    echo "REFERENCESUBSTRUCTUREFILE: $crtMvRef" >> "$accGeomEdtParFile"
    echo "OUTFILE: $accGeomEdtOutFile" >> "$accGeomEdtParFile"
    # Launch ACC
    autocompchem -p "$accGeomEdtParFile" > "$accGeomEdtLogFile"
    # Check outcome
    if [ ! -f "$accGeomEdtLogFile" ]; then
        errMsg="#ACCGeomEdit: $accGeomEdtLogFile not found."
        abandon "$lowestEnergyConfSDF" "$E_OPTERROR"
    fi
    if [ ! -f "$accGeomEdtOutFile" ]; then
        errMsg="#ACCGeomEdit: $accGeomEdtOutFile not found."
        abandon "$lowestEnergyConfSDF" "$E_OPTERROR"
    fi
    if ! grep -q "Termination status: 0" "$accGeomEdtLogFile" ; then
        errMsg="#ACCGeomEdit: non-zero exit status from AutoCompChem."
        abandon "$lowestEnergyConfSDF" "$E_OPTERROR"
    fi

   
    #
    # Split the modified geometries
    #
    modGeomFileRoot="$wrkDir/${molName}_stp"
    obabel -isdf "$accGeomEdtOutFile" -osdf -O "$modGeomFileRoot.sdf" -m
    # Check outcome
    if [ "$?" -ne 0 ]; then
        errMsg="#OBabelSplit: non-zero exit status from OpenBabel."
        abandon "$lowestEnergyConfSDF" "$E_OPTERROR"
    fi


    #
    # Refine all modified geometries
    #
    sumNrgModGeoms="0.0"
    allModGeomRefinedSDF="$wrkDir/${molName}_allRefStp.sdf"
    modGeomFiles=($(ls -d "$modGeomFileRoot"*))
    for modGeomInpSDF in ${modGeomFiles[@]}
    do
        n=$((n+1))
        #
        # Prepare input semi-empirical refinement
        #
        echo "$n-Prepare input for semi-empirical geometry optimization..."
        # Setting new params
        accSprGOSEParFile="$wrkDir/${molName}_SprtGOSE${n}In.par"
        accSprGOSELog="$wrkDir/${molName}_SprtGOSE${n}In.log"
        accSprGOSEInp="$wrkDir/${molName}_SprtGOSE${n}In.sdf"
        sprtGOSEInp="$wrkDir/${molName}_GOSE${n}.spartan"
        cp "$modGeomInpSDF" "$accSprGOSEInp"
        portablesed "s/${molName}.*/${molName}_GOSE${n}/g" "$accSprGOSEInp"
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprGOSEParFile"
        echo "TASK: PrepareInputSpartan" >> "$accSprGOSEParFile"
        echo "INFILE: $accSprGOSEInp" >> "$accSprGOSEParFile"
        echo "KEYWORDS: OPT PM6 PRINTLEV=1 OPTCYCLE=2600" >> "$accSprGOSEParFile"
        echo "CHARGE: $charge" >> "$accSprGOSEParFile"
        echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprGOSEParFile"
        echo "OUTFILE: $sprtGOSEInp" >> "$accSprGOSEParFile"
        echo "\$STARTFREEZEATOMS: $sprtGOSEFrozen" >> "$accSprGOSEParFile"
        echo "\$END" >> "$accSprGOSEParFile"
        # Launch ACC
        autocompchem -p "$accSprGOSEParFile" > "$accSprGOSELog"
        # Check outcome
        if [ ! -f "$accSprGOSELog" ]; then
            errMsg="#MakeSpartanGOSEInput: $accSprGOSELog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -d "$sprtGOSEInp" ]; then
            errMsg="#MakeSpartanGOSEInput: $sprtGOSEInp directory not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprGOSELog" ; then
            errMsg="#MakeSpartanGOSEInput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi


        #
        # Submit semi-empirical refinement
        #
        echo "$n-Submitting Spartan semi-empirical optimization..."
        # Setting new parameters
        sprtGOSEOutLog="$sprtGOSEInp/${molName}_GOSE${n}/output"
        sprtGOSELog="$wrkDir/${molName}_SprtGOSE${n}Run.log"
        # Lauch Spartan
        "$SPARTANEXE" --foreground-submit "$sprtGOSEInp" > "$sprtGOSELog"
        # Check outcome
        if [ ! "$?" == 0 ]; then
            errMsg="#RunSpartanGOSE: non-zero exit status from Spartan call."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtGOSEOutLog" ]; then
            errMsg="#RunSpartanGOSE: $sprtGOSEOutLog file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Successful completion" "$sprtGOSEOutLog" ; then
            errMsg="#RunSpartanGOSE: unsuccessful completion."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi


        #
        # Extract optimized geometry
        #
        echo "$n-Extracting results from Spartan semi-empirical geometry optimization..."
        # Setting new params
        accSprGOSEOutParFile="$wrkDir/${molName}_SprtGOSE${n}Out.par"
        accSprGOSEOutLog="$wrkDir/${molName}_SprtGOSE${n}Out.log"
        sprtGOSEOut="$wrkDir/${molName}_SprtGOSE${n}Out.sdf"
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprGOSEOutParFile"
        echo "TASK: ExtractLastGeometryFromSpartanTree" >> "$accSprGOSEOutParFile"
        echo "INFILE: $sprtGOSEInp" >> "$accSprGOSEOutParFile"
        echo "OUTFORMAT: SDF" >> "$accSprGOSEOutParFile"
        echo "OUTFILE: $sprtGOSEOut" >> "$accSprGOSEOutParFile"
        # Launch ACC
        autocompchem -p "$accSprGOSEOutParFile" > "$accSprGOSEOutLog"
        # Check outcome
        if [ ! -f "$accSprGOSEOutLog" ]; then
            errMsg="#ExtractSpartanGOSEOutput: $accSprGOSEOutLog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtGOSEOut" ]; then
            errMsg="#ExtractSpartanGOSEOutput: $sprtGOSEOut file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprGOSEOutLog" ; then
            errMsg="#ExtractSpartanGOSEOutput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        portablesed "1 s/^${molName}.*$/$molName/g" "$sprtGOSEOut"
        portablesed "1 s/^.*$/$molName/g" "$sprtGOSEOut"


        #
        # Extract energy
        #
        if grep -q " SCF did not converge - abort optimization" "$sprtGOSEOutLog"
        then
            pm6GOEnergy=1000000000
            echo "$n-PM6 Energy not found (FAILED SCF)"
        else
            if grep -q " Heat of Formation: " "$sprtGOSEOutLog"
            then
            echo "$n-Extracting PM6 energy of optimized geometry..."
                pm6GOEnergy=$(grep " Heat of Formation: " "$sprtGOSEOutLog" | tail -n 1 | awk '{print $4}')
                pm6GOEnergy=$(echo ${pm6GOEnergy} | sed -e 's/[eE]+*/\*10\^/')
                pm6GOEnergy=$(bc -l <<< "($pm6GOEnergy)")
                echo "$n-PM6 Energy of refined conformation: $pm6GOEnergy (from GOSE$n)"
                addPropertyToSingleMolSDF "RefinedStep-PM6_ENERGY" "$pm6GOEnergy" "$sprtGOSEOut"
                sumNrgModGeoms=$(bc -l <<< "($pm6GOEnergy) + ($sumNrgModGeoms)")
            else
                errMsg="#ExtractSpartanEnergy: cannot find heat of formation."
                abandon "$DnCG3Dout" "$E_OPTERROR"
            fi
        fi

        #
        # Put all in one file
        #
        cat "$sprtGOSEOut" >> "$allModGeomRefinedSDF"
    done


    #
    # Prepare output
    #
    echo "Storing sum of PM6 energies is $sumNrgModGeoms as FrozenCore-PM6_ENERGY."
    cp "$allModGeomRefinedSDF" "$outSDF"
    addPropertyToMultiMolSDF "FrozenCore-PM6_ENERGY" "$sumNrgModGeoms" "$outSDF"

    #
    # Clean up
    #
    if [ "$cleanup" == 0 ]
    then
        rm -f "$accGeomEdtParFile"
        rm -f "$accGeomEdtLogFile"
        rm -f "$accGeomEdtOutFile"
        rm -f "allModGeomRefinedSDF"
    fi
else
    #
    # Prepare output
    #
    echo "Preparing output from $lowestEnergyConfSDF"
    cp "$lowestEnergyConfSDF" "$outSDF"
    addPropertyToSingleMolSDF "FrozenCore-PM6_ENERGY" "$lowestEnergy" "$outSDF"
fi


#
# Free-up token is done by function 'finish'
#
#releaseLock "$lockIDFile" "$lockFile"


#
# Cleanup
#
if [ "$cleanup" == 0 ]
then
    echo "Cleaning up"
    rm -f "$locInpSDF"
    rm -f  "$DnCG3Dout"
    rm -f  "$DnCGParFile"
    rm -f  "$accSprGOParFile"
    rm -f  "$accSprGOLog"
    rm -f  "$accSprGOInp"
    rm -f  "$sprtGOLog"
    rm -f  "$accSprGOOutParFile"
    rm -f  "$accSprGOOutLog"
    rm -f  "$sprtGOOut"
    rm -fr  "$sprtGOInp"
    rm -fr  "$sprtGOSEInp"
    rm -f  "$accSprDYNParFile"
    rm -fr  "$sprtDYNInp"
    rm -f  "$accSprDYNLog"
    rm -f  "$accSprDYNInp"
    rm -fr  "$sprtDYNOutDir"
    rm -f  "$sprtDYNLog"
    for i in $(seq 0 $((n+1)))
    do
        rm -f  "$wrkDir/${molName}_SprtDYNOut${i}.par"
        rm -f  "$wrkDir/${molName}_SprtDYNOut${i}.log"
        rm -f  "$wrkDir/${molName}_SprtDYNOut${i}.sdf"
        rm -f  "$wrkDir/${molName}_SprtGOSE${i}In.par"
        rm -f  "$wrkDir/${molName}_SprtGOSE${i}In.log"
        rm -f  "$wrkDir/${molName}_SprtGOSE${i}In.sdf"
        rm -f  "$wrkDir/${molName}_SprtGOSE${i}Run.log"
        rm -f  "$wrkDir/${molName}_SprtGOSE${i}Out.par"
        rm -f  "$wrkDir/${molName}_SprtGOSE${i}Out.log"
        rm -f  "$wrkDir/${molName}_SprtGOSE${i}Out.sdf"
        rm -fr  "$wrkDir/${molName}_GOSE${i}.spartan"
        rm -fr  "$wrkDir/${molName}_GOMM${i}.spartan"
    done
    echo "Done!"
fi

#
# Create task completed label
#
writeTCL "$jobId"

exit 0
