#!/bin/bash
###############################################################################
#
# Script for construction of the 3D models and conformational search.
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

#Fragmentation
coreLigandCutRules="$WORKDIR/core-lig_CutRule.rul"

#Tinker (with extensions)
pathToTinkerBin="$TINKERBIN"

#Label and extension of input file
inputLabelExt="_inp"

# RootName for lockfiles that restrain parallel executions of this script
rootLock="/tmp/runConfSearchOrRefine"
rootLockID="/tmp/activeConfSearchOrRefine"
lckMaxWait=10000  # we will wait for a lockfile for up to this number of steps
lockStep=30       # length of waiting step
lockTimeUnit="s"  # time unit: s (seconds), m (minutes), h (hours)

# Fragment space
compMatrix="$WORKDIR/CPMap_FS-C5.0.par"
libFrags="$WORKDIR/lib_frags_FS_C5.0_v3.sdf"
libCaps="$WORKDIR/lib_caps_FS_C5.0_v3.sdf"
libScaff="notSet"
libScaffA1="$WORKDIR/scaff_A1_v3.sdf"
libScaffA2="$WORKDIR/scaff_A2_v3.sdf"
libScaffB1="$WORKDIR/scaff_B1_v3.sdf"
libScaffC1="$WORKDIR/scaff_C1_v3.sdf"
libScaffD1="$WORKDIR/scaff_D1_v3.sdf"
libScaffB2="$WORKDIR/scaff_B2_v3.sdf"
libScaffC2="$WORKDIR/scaff_C2_v3.sdf"
libScaffC3="$WORKDIR/scaff_C3_v3.sdf"
libScaffD2="$WORKDIR/scaff_D2_v3.sdf"
libScaffE1="$WORKDIR/scaff_E1_v3.sdf"
libScaffF1="$WORKDIR/scaff_F1_v3.sdf"
libScaffH1="$WORKDIR/scaff_H1_v3.sdf"
libScaffE2="$WORKDIR/scaff_E2_v3.sdf"
libScaffF2="$WORKDIR/scaff_F2_v3.sdf"
libScaffH2="$WORKDIR/scaff_H2_v3.sdf"
libScaffL0="$WORKDIR/scaff_L0_v3.sdf"
libScaffP1="$WORKDIR/scaff_P1_v3.sdf"
libScaffP2="$WORKDIR/scaff_P2_v3.sdf"

# state (i.e., A,B,C...) specific molecular parameters (may need to add chargeA, chargeB,...)
charge=0
spinMult=1

# 3D builder
tinkerForceField="$WORKDIR/uff_vdw.prm"
rotSpaceDef="$WORKDIR/rotatableBonds-1.2"
tinkerKeyFile="$WORKDIR/build_uff_ConfSearch.key"
tinkerSubmitFile="$WORKDIR/submit_ConfSearch"

# Spartan single point energy calculations 
#sprtSPFF="$WORKDIR/params.MMFF94_CS-2.2" we use a global params file ~/PARAMS.MMFF.DENOPTIM

sprtSPAtmSpecKey="
[Ru] FFHINT= ~~210
[\$([#6;X3]([#1])(~[Ru])~[#6;X3](@~[#6;X3])@~[#6;X3])] FFHINT= ~~206
[\$([#6;X4]([#1])([#1])~[Ru])] FFHINT= ~~207
[\$([#7;X3](!@~[c])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#6;X4]),\$([#7;X3](!@~[c])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#1]),\$([#7;X3](!@~[c])~[#6;X3](~[Ru])~[#7;X3](@~[#6;X3])~[#6;X4])] FFHINT= ~~208
[\$([#7;X3](!@~[c])~[#6](~[Ru])~[#6;X4](~[#6])(~[#6])~[#6]),\$([#7;X3]1~[#6;X3](~[Ru])~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2~[#6]~1)] FFHINT= ~~208
[\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6;X3](~[Ru])~[#7;X3]!@~[c]),\$([#7;X3](~[#6;X4])([#1])~[#6;X3](~[Ru])~[#7;X3]!@~[c]),\$([#7;X3](~[#6;X4])(@~[#6;X3])~[#6;X3](~[Ru])~[#7;X3]!@~[c])] FFHINT= ~~209
[\$([#7;X3](!@~[c])~[#6;X3](~[Ru])~[#7;X3]!@~[c]),\$([#7;X3](!@~[c])~[#6;X3](~[Ru])~[#7;X3]([#1])[#1]),\$([#7;X3]([#1])([#1])~[#6;X3](~[Ru])~[#7;X3]!@~[c])] FFHINT= ~~211
[\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#6;X4]),\$([#7;X3](!@~[#6;X4])(@~[#6;X3])@~[#6;X3](~[Ru])@~[#7;X3](@~[#6;X3])!@~[#6;X4]),\$([#7;X3](~[#6;X4])([#1])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#6;X4]),\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6;X3](~[Ru])~[#7;X3]([#1])~[#6;X4]),\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6;X3](~[Ru])~[#7;X3]([#1])[#1]),\$([#7;X3]([#1])([#1])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#6;X4]),\$([#7;X3](~[#6;X4])([#1])~[#6;X3](~[Ru])~[#7;X3]([#1])~[#6;X4]),\$([#7;X3](~[#6;X4])([#1])~[#6;X3](~[Ru])~[#7;X3]([#1])[#1]),\$([#7;X3]([#1])([#1])~[#6;X3](~[Ru])~[#7;X3](~[#6;X4])~[#1]),\$([#7;X3]([#1])([#1])~[#6;X3](~[Ru])~[#7;X3]([#1])[#1]),\$([#7;X3](~[#6;X4])(~[#6;X4])~[#6](~[Ru])~[#6;X4](~[#6])(~[#6])~[#6]),\$([#7;X3](~[#6;X4])(~[#1])~[#6](~[Ru])~[#6;X4](~[#6])(~[#6])~[#6])] FFHINT= ~~212
[\$([#7;X3]~[#6;X2]~[#7;X3]),\$([#7;X3]~[#6;X2]~[#6;X4](~[#6])(~[#6])~[#6]),\$([#7;X3]1~[#6;X2]~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2~[#6]~1)] FFHINT= ~~212"

# Spartan conformation-constrained geometry optimization
#sprtGOFF="$WORKDIR/params.MMFF94_CS-2.2" we use a global params file ~/PARAMS.MMFF.DENOPTIM
sprtGOFrozen="
[#7;X3]~[#6;X2]~[#7;X3]
[*]~[Ru]
[\$([#6,#1]~[#6]~[Ru]);!\$([#6,#1]~[#6](~[Ru])~[#7])]
[\$([*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[#6](~[Ru])~[#7]);!\$([*]~[#7]~[#6](~[Ru])~[#6,#7])]
[\$([*]~[*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[*]~[#6]~[#7]);!\$([*]~[*]~[#7]~[#6](~[Ru])~[#6,#7])]
[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8](~[Ru])~[#6]
[#6]([#6])([#6])~[#8]~[Ru]"
sprtGOAtmSpecKey="$sprtSPAtmSpecKey"
sprtConstrGO="
[Ru] [\$([#6](~[Ru])@~[#7])] [\$([#7]@~[#6]~[Ru])] [\$([#6]!@~[#7]@~[#6]~[Ru])]"

# Spartan conformational search 
numConformersSparse=500     # max num. of conformers examined
numConformersMC=500     # max num. of conformers examined
maxSparseConformers=5 # max num. of conformers kept and submitted to next step
#sprtCSFF="$WORKDIR/params.MMFF94_CS-2.2" we use a global params file ~/PARAMS.MMFF.DENOPTIM
sprtRotBnd="[#6]!@-[#6] 3
[#7]!@~[#6] 3
[#8,X2]!@-[#6] 3
[\$([#6]1-[#6]-[*]-[*]-[*]-[*]-[#6]1)] 2
[\$([#6]1-[#6]-[*]-[*]-[*]-[#6]1)] 4
[\$([#6]1-[#6]-[*]-[*]-[#6]1)] 2"
sprtCSAtmSpecKey="$sprtSPAtmSpecKey"
sprtConstrCS="[#7;X3] [#6;X2]
[#7;X3] [#6;X2] [#7;X3]
[#7;X3] [#6;X2] [#7;X3] [*]

[Ru] [*]
[\$([*]~[Ru])] [\$([*]~[*]~[Ru])]
[\$([*]~[*]~[Ru]);!\$([*]~[#6]~[#7])] [\$([*]~[*]~[*]~[Ru]);!\$([*]~[*]~[#6]~[#7])]
[\$([*]~[*]~[*]~[Ru]);!\$([*]~[*]~[#6]~[#7])] [\$([*]~[*]~[*]~[*]~[Ru]);!\$([*]~[*]~[*]~[#6]~[#7])]
[\$([*]~[*]~[*]~[*]~[Ru]);!\$([*]~[*]~[*]~[#6]~[#7])] [\$([*]~[*]~[*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[*]~[*]~[#6]~[#7])]
[\$([*]~[*]~[*]~[*]~[*]~[Ru]);!\$([*]~[*]~[*]~[*]~[#6]~[#7])] [\$([*]~[*]~[*]~[*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[*]~[*]~[*]~[#6]~[#7])]

[*] [Ru] [*]
[*] [\$([#6]~[Ru]);!\$([#6]~[#7])] [*]
[*] [\$([#8]~[Ru])] [*]
[*] [\$([#1]~[Ru])] [*]
[*] [\$([#6]~[#6]~[Ru]);!\$([#6]~[#6]~[#7])] [*]
[*] [\$([#6]~[#6]~[#6]~[Ru]);!\$([#6]~[#6]~[#6]~[#7])] [*]
[*] [\$([#6]~[#6]~[#6]~[#6]~[Ru]);!\$([#6]~[#6]~[#6]~[#6]~[#7])] [*]

[*] [\$([#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8](~[Ru])~[#6]),\$([#6]1~[#6]~[#6]~[#6]~[#6](~[#8](~[Ru])~[#6])~[#6]~1),\$([#6]1~[#6]~[#6](~[#8](~[Ru])~[#6])~[#6]~[#6]~[#6]~1)] [\$([#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8](~[Ru])~[#6]),\$([#6]1~[#6]~[#6]~[#6]~[#6](~[#8](~[Ru])~[#6])~[#6]~1),\$([#6]1~[#6]~[#6](~[#8](~[Ru])~[#6])~[#6]~[#6]~[#6]~1)] [*]
[\$([*]~[Ru]);!\$([#6](~[#7])~[Ru])] [Ru] [\$([#6](~[#7])~[Ru])] [\$([#6,#7]~[#6](~[#7])~[Ru])]"
#Possibly constrained torsions in CS
# [\$([#17]~[Ru])] [Ru] [\$([#6]@~[Ru])] [\$([*]~[#6]@~[Ru])]
# [\$([#1]~[Ru])] [Ru] [\$([#6]@~[Ru])] [\$([*]~[#6]@~[Ru])]
# [\$([#6]~[Ru])] [Ru] [\$([#6]@~[Ru])] [\$([*]~[#6]@~[Ru])]
# [Ru] [\$([#6]1~[Ru]~[#6]~[#6]~1)] [\$([#6]1~[#6]~[Ru]~[#6,#1]~1)] [*]
# [Ru] [\$([#6](~[Ru])~[#7])] [\$([#7]~[#6]~[Ru])] [\$([#6]!@~[#7]@~[#6]~[Ru])]

# Spartan generation of acyclic Carbene-N-R conformers
sprtDYNConstraints="
[\$([#6](~[#1])~[Ru])] [Ru] [\$([#6](~[#7])~[Ru])] [\$([#7]~[#6]~[Ru])]  0.0 0.0 180.0 1
[Ru] [\$([#6](!@~[#7])~[Ru])] [\$([#7]!@~[#6]~[Ru])] [\$([#6]~[#7]!@~[#6]~[Ru])]  0.0 0.0 180.0 1
[\$([#6]~[#15])] [#15] [Ru] [\$([#6](~[#1])~[Ru])] 0.0 0.0 180 1"
# ^Inclusion of PCy3 and PPh3

# Spartan semi-empirical refinement
sprtGOSEFrozen="
[#17,#8,#1]~[Ru]
[\$([#6]~[Ru]);!\$([#6](~[Ru])~[#7])]
[\$([#6,#1]~[#6]~[Ru]);!\$([#6,#1]~[#6](~[Ru])~[#7])]
[\$([*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[#6](~[Ru])~[#7]);!\$([*]~[#7]~[#6](~[Ru])~[#6,#7])]
[\$([*]~[*]~[*]~[#6]~[Ru]);!\$([*]~[*]~[*]~[#6]~[#7]);!\$([*]~[*]~[#7]~[#6](~[Ru])~[#6,#7])]
[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8](~[Ru])~[#6]
[#6]~[#8]~[Ru]"

# Cleanup 0:do it, 1:don't
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
# Function to change title of a single molecule SDF file
# @param $1 new title
# param $2 the SDF filename
#

function changeTitleToSingleMolSDF(){
    newName="$1"
    file="$2"
    n=$(wc -l "$file" | awk '{print $1}')
    n=$((n-=1))
    tail -n "$n" "$file" > "${file}_tmp"
    echo "$newName" > "$file"
    cat "${file}_tmp" >> "$file"
    rm "${file}_tmp"
    removePropertyFromSDF "cdk:Title" "$file"
    addPropertyToSingleMolSDF "cdk:Title" "$newName" "$file"
}

#
# Function to append a property to an SDF file. Does not overwrite any existing
# property in the SDF.
# @param $1 the name of the property to be removed
# @param $3 the SDF file to alter
#
function removePropertyFromSDF() {
    propName="$1"
    file="$2"
    awk -v nlines=2 -v ptrn=$propName '$0~ptrn {for (i=0; i<nlines; i++) {getline}; next} 1' "$file" > "${file}_tmp"
    mv "${file}_tmp" "$file"
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
    n=$(grep -n '$$$$' "$file" | cut -f1 -d:)
    n=$((n-=1))
    head -n "$n" "$file" > "${file}_tmp"
    echo "> <$propName>" >> "${file}_tmp"
    echo "$propValue" >> "${file}_tmp"
    echo "" >> "${file}_tmp"
    echo '$$$$' >> "${file}_tmp"
    mv "${file}_tmp" "$file"
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
    errMsg="#3DBuild: wrong number of arguments"
    tcl="${tcl}PID$$"
    exit "$E_FATAL"
fi
inpSDF="$1"
outSDF="$2"
wrkDir="$3"
jobId="$4"
locDir="$(pwd)"
tcl="$wrkDir/$tcl$jobId"

molName=`basename "$inpSDF" .sdf`
molNum=`basename "$molName" "$inputLabelExt"`

#
# Setup Log file
#
log="$jobId.log"
exec > $log
exec 2>&1
echo "Log for 3D-model Build and Conformational Search" 
echo "Input: $inpSDF" 
echo "WorkDir: $wrkDir"
echo "Output: $outSDF"
echo "JobID: $jobId"
echo "TCL: $tcl"

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
case "$jobTyp$stereo" in
    "A1") 
        libScaff="$libScaffA1" 
        ;;
    "A2")
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
        ;;
    "C2")
        libScaff="$libScaffC2"
        ;;
    "D1") 
        libScaff="$libScaffD1"
        ;;
    "D2")
        libScaff="$libScaffD2"
        ;;
    "E1") 
        libScaff="$libScaffE1"
        ;;
    "E2")
        libScaff="$libScaffE2"
        ;;
    "F1") 
        libScaff="$libScaffF1"
        ;;
    "F2")
        libScaff="$libScaffF2"
        ;;
    "H1") 
        libScaff="$libScaffH1"
        ;;
    "H2")
        libScaff="$libScaffH2"
        ;;
    "L0")
        libScaff="$libScaffL0"
        ;;
    "P1")
        libScaff="$libScaffP1"
        ;;
    "P2")
        libScaff="$libScaffP2"
        ;;
    *) 
        errMsg="#3DBuild: error interpreting jobId $jobTyp$stereo"
        abandon "$inpSDF" "$E_FATAL"
        ;;
esac
echo "Using scaffolds library $libScaff"
molName="${molNum}_$jobTyp${stereo}"
echo "Molecule name set to $molName"


#
# Build 3D model
#
echo "Starting DenoptimCG..."
# Setting new params
DnCG3Dinp="$wrkDir/${molNum}$jobTyp${stereo}_CGinp.sdf"
DnCG3Dout="$wrkDir/${molName}_3D-CG.sdf"
DnCGParFile="$wrkDir/${molName}_DnCG.par"
DnGEParFile="$wrkDir/replaceScaff_$molName.params"
# Changing scaffold with denoptim GE
# Setting GE params
echo "GRAPHEDIT-INPUTGRAPHS=$wrkDir/$inpSDF" >> "$DnGEParFile"
echo "GRAPHEDIT-GRAPHSEDITSFILE=$WORKDIR/replaceScaffoldEditingTask" >> "$DnGEParFile"
echo "GRAPHEDIT-OUTPUTGRAPHS=$DnCG3Dinp" >> "$DnGEParFile"
echo "FS-scaffoldLibFile=$libScaff" >> "$DnGEParFile"
# Launch denoptim GE
echo Calling denoptim -r GE
denoptim -r GE "$DnGEParFile" > "${molName}_GE.log" 2>&1
# Check completed status from GE log
if [ -f "$DnCG3Dinp" ]; then
    if ! grep -q 'Completed GraphEditor' "${molName}_GE.log"
    then
        echo "Exiting! Missing job completed status in output log from Denoptim GE"
        cat "$DnCG3Dinp" > "$outSDF"
        writeTCL "$jobId"
        releaseLock "$lockFile" "$lockIDFile"
        exit "$E_OPTERROR"
    fi
else
    echo "$DnCG3Dinp not found."
    errMsg="#DenoptimGE: $DnCG3Dinp not found."
    abandon "$inpSDF" "$E_OPTERROR"
fi
# Clean up GE
rm -f "${molName}_GE.log" "$DnGEParFile"  
changeTitleToSingleMolSDF "${molNum}-$jobTyp${stereo}" "$DnCG3Dinp"
# Prepare param file
echo "3DB-VERBOSITY=1" > "$DnCGParFile"
echo "3DB-inpSDF=$DnCG3Dinp" >> "$DnCGParFile"
echo "3DB-outSDF=$DnCG3Dout" >> "$DnCGParFile"
echo "FS-scaffoldLibFile=$libScaff" >> "$DnCGParFile"
echo "FS-fragmentLibFile=$libFrags" >> "$DnCGParFile"
echo "FS-cappingFragmentLibFile=$libCaps" >> "$DnCGParFile"
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

# Launch DenoptimB3D
denoptim -r B3D "$DnCGParFile"
# Cleanup files tmp files
echo "Removing $wrkDir/$molNum$jobTyp${stereo}_cs0.*"
rm -f "$wrkDir/$molNum$jobTyp${stereo}_cs0."*
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
# Further handling of only for free ligand 
#
if [ "$jobTyp" == "L" ]
then
    #
    # Chop system
    #
    echo "Extracting only desired portion..."
    absMinSDF="$DnCG3Dout"
    gm3dfInp="$wrkDir/${molNum}-$jobTyp${stereo}_extPart.sdf"
    gm3dfPar="$wrkDir/${molNum}-$jobTyp${stereo}_extPart.par"
    gm3dfLog="$wrkDir/${molNum}-$jobTyp${stereo}_extPart.log"
    ligFrag="$wrkDir/${molNum}-$jobTyp${stereo}_extPartLigFrag.sdf"
    cp "$absMinSDF" "$gm3dfInp"
    echo "FRG-STRUCTURESFILE=$gm3dfInp" > "$gm3dfPar"
    echo "FRG-VERBOSITY=1" >> "$gm3dfPar"
    echo "FRG-CUTTINGRULESFILE=$coreLigandCutRules" >> "$gm3dfPar"
    echo "FRG-MAXFRAGSIZE=5000" >> "$gm3dfPar"
    echo "FRG-REJECTELEMENT=Ru" >> "$gm3dfPar"
    # Submit GM3DFragmenter
    cd "$wrkDir"
    denoptim -r FRG "$gm3dfPar" > "$gm3dfLog" 2>&1 
    cd "$locDir"
    rm -f "$wrkDir/MolFrag-ratio_${molNum}-$jobTyp${stereo}_extPart"*
    rm -f "$wrkDir/CPMap_${molNum}-$jobTyp${stereo}_extPart"*
    # Check outcome
    if [ ! -f "$gm3dfLog" ]; then
        errMsg="#GM3DFragmenter: $gm3dfLog not found."
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    if ! grep -q " Completed Fragmenter" "$gm3dfLog" ; then
        errMsg="#GM3DFragmenter: non-zero exit status"
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    fragmenterDir="$( grep "Output" "$gm3dfLog" | awk '{print $11}' )"
    numStoredFrags="$( grep "\$\$\$\$" "$fragmenterDir/Fragments.sdf" | wc -l )"
    if [[ "$numStoredFrags" -eq "1" ]] ; then
        echo "Ligand fragment correctly extracted"
    else
        errMsg="#GM3DFragmenter: ligand fragment not extracted properly"
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    

    #
    # Convert fragment into molecule for further processing
    #
    mv "$fragmenterDir/Fragments.sdf" "$DnCG3Dout"
    removePropertyFromSDF "CLASS" "$DnCG3Dout"
    removePropertyFromSDF "ATTACHMENT_POINT" "$DnCG3Dout"
    removePropertyFromSDF "cdk:Remark" "$DnCG3Dout"
fi
changeTitleToSingleMolSDF "$molName" "$DnCG3Dout"

#
# Check for exceeding atom crowding
#
if [ ! "$jobTyp" == "F" ] ; then
    echo "Starting Atom Clash detector..."
    # Setting new params
    atmClshParFile="$wrkDir/${molName}_AtmClsh.par"
    atmClshLog="$wrkDir/${molName}_AtmClsh.log"
    # Prepare param file
    echo "VERBOSITY: 1" > "$atmClshParFile"
    echo "TASK: AnalyzeVDWClashes" >> "$atmClshParFile"
    echo "INFILE: $DnCG3Dout" >> "$atmClshParFile"
    echo "ALLOWANCE13: 1.0" >> "$atmClshParFile"
    echo "CUTOFF: 0.85" >> "$atmClshParFile"
    echo "ALLOWANCE: 0.50" >> "$atmClshParFile"
    # Launch ACC
    autocompchem -p "$atmClshParFile" > "$atmClshLog"
    # Check outcome
    if [ ! -f "$atmClshLog" ]; then
        errMsg="#AtomClash: $atmClshLog not found."
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    #atomClashF="0"
    if ! grep -q "Termination status: 0" "$atmClshLog" ; then
        #if [ "$jobTyp" == "F" ]; then
        #atomClashF="1"
        #echo "Warning: Atom clash detected in F."
        #else
        errMsg="#AtomClash: non-zero exit status from AutoCompChem"
        abandon "$DnCG3Dout" "$E_OPTERROR"
        #fi
    fi
    if grep -q 'Found 0 mols with one or more atom clashes' "$atmClshLog" ; then
        echo "No atom clashes"
    #elif [ "$atomClashF" == "1" ] ; then
    #    echo "Atom clash found in F."
    else
        errMsg="#AtomClash: Found Atom Clashes"
        abandon "$DnCG3Dout" "$E_OPTERROR"
    fi
fi

#
# Prepare geometry optimization
#
echo "Prepare first MMFF geometry optimization with frozen core..."
# Setting new params
accSprGOParFile="$wrkDir/${molName}_SprtGOIn.par"
accSprGOLog="$wrkDir/${molName}_SprtGOIn.log"
accSprGOInp="$wrkDir/${molName}_SprtGOIn.sdf"
sprtGOInp="$wrkDir/${molName}_GO.spartan"
cp "$DnCG3Dout" "$accSprGOInp"
if [ "$SEDSYNTAX" == "GNU" ]
then
    sed -i "s/${molName}.*/${molName}_GO/g" "$accSprGOInp"
elif [ "$SEDSYNTAX" == "BSD" ]
then
    sed -i '' "s/${molNum}.*/${molName}_GO/g" "$accSprGOInp"
fi
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
# Submit geometry optimization
#
echo "Submitting Spartan MMFF geometry optimization..."
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
echo "Extracting results from Spartan first MMFF geometry optimization..."
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
if [ "$SEDSYNTAX" == "GNU" ]
then
    sed -i "1 s/^.*$/$molName/g" "$sprtGOOut"
elif [ "$SEDSYNTAX" == "BSD" ]
then
    sed -i '' "1 s/^.*$/$molName/g" "$sprtGOOut"
fi


#
# Extract energy
#
echo "Extracting MMFF energy of geometry from pre. geom. optimization..."
mmffGOEnergy=$(grep " Totals " "$sprtGOOutLog" | tail -n 1 | awk '{print $2}')
mmffGOEnergy=$(echo ${mmffGOEnergy} | sed -e 's/[eE]+*/\*10\^/')
mmffGOEnergy=$(bc -l <<< "($mmffGOEnergy)")


#
# Prepare input for single point energy calculation
#
echo "Prepare input for single point MMFF energy calculation..."
# Setting new params
accSprSP0ParFile="$wrkDir/${molName}_SprtSP0In.par"
accSprSP0Log="$wrkDir/${molName}_SprtSP0In.log"
accSprSP0Inp="$wrkDir/${molName}_SprtSP0In.sdf"
sprtSP0Inp="$wrkDir/${molName}_SP0.spartan"
cp "$sprtGOOut" "$accSprSP0Inp"
if [ "$SEDSYNTAX" == "GNU" ]
then
    sed -i "s/${molName}.*/${molName}_SP0/g" "$accSprSP0Inp"
elif [ "$SEDSYNTAX" == "BSD" ]
then
    sed -i '' "s/${molNum}.*/${molName}_SP0/g" "$accSprSP0Inp"
fi
# Prepare param file
echo "VERBOSITY: 1" > "$accSprSP0ParFile"
echo "TASK: PrepareInputSpartan" >> "$accSprSP0ParFile"
echo "INFILE: $accSprSP0Inp" >> "$accSprSP0ParFile"
echo "KEYWORDS: SOPT MMFF PRINTLEV=3 USEFFPARMS=~/PARAMS.MMFF.DENOPTIM" >> "$accSprSP0ParFile"
echo "CHARGE: $charge" >> "$accSprSP0ParFile"
echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprSP0ParFile"
echo "OUTFILE: $sprtSP0Inp" >> "$accSprSP0ParFile"
echo "\$STARTATOMSPECIFICKEYWORDS: $sprtSPAtmSpecKey" >> "$accSprSP0ParFile"
echo "\$END" >> "$accSprSP0ParFile"
# Launch ACC
autocompchem -p "$accSprSP0ParFile" > "$accSprSP0Log"
# Check outcome
if [ ! -f "$accSprSP0Log" ]; then
    errMsg="#MakeSpartanSP0Input: $accSprSP0Log not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -d "$sprtSP0Inp" ]; then
    errMsg="#MakeSpartanSP0Input: $sprtSP0Inp directory not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Termination status: 0" "$accSprSP0Log" ; then
    errMsg="#MakeSpartanSP0Input: non-zero exit status from AutoCompChem."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi


#
# Submit Single point energy calculation
#
echo "Submitting Spartan single point MMFF energy calculation ..."
# Setting new parameters
sprtSP0OutLog="$sprtSP0Inp/${molName}_SP0/output"
sprtSP0Log="$wrkDir/${molName}_SprtSP0Run.log"
# Lauch Spartan
"$SPARTANEXE" --foreground-submit "$sprtSP0Inp" > "$sprtSP0Log"
# Check outcome
if [ ! "$?" == 0 ]; then
    errMsg="#RunSpartanSP0: non-zero exit status from Spartan call."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -f "$sprtSP0OutLog" ]; then
    errMsg="#RunSpartanSP0: $sprtSP0OutLog file not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Successful completion" "$sprtSP0OutLog" ; then
    errMsg="#RunSpartanSP0: unsuccessful completion."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi


#
# Extract energy
#
echo "Extracting MMFF energy of geometry from pre. geom. optimization..."
mmffGOEnergySP=$(grep " Totals " "$sprtSP0OutLog" | tail -n 1 | awk '{print $2}')
mmffGOEnergySP=$(echo ${mmffGOEnergySP} | sed -e 's/[eE]+*/\*10\^/')
mmffGOEnergySP=$(bc -l <<< "($mmffGOEnergySP)")
echo "SP Energy from GO: $mmffGOEnergySP (from GO $mmffGOEnergy)"


#
# Prepare input for SPARSE conformational search
#
echo "Preparing input for SPARSE conformational search with Spartan..."
# Setting new params
accSprCSSPRParFile="$wrkDir/${molName}_SprtCSSPRIn.par"
accSprCSSPRLog="$wrkDir/${molName}_SprtCSSPRIn.log"
sprtCSSPRInp="$wrkDir/${molName}_CSSPR.spartan"
accSprCSSPRInp="$wrkDir/${molName}_SprtCSSPRIn.sdf"
rnd=$(echo $RANDOM)
cp "$sprtGOOut" "$accSprCSSPRInp"
if [ "$SEDSYNTAX" == "GNU" ]
then
    sed -i "s/${molName}.*/${molName}_CSSPR/g" "$accSprCSSPRInp"
elif [ "$SEDSYNTAX" == "BSD" ]
then
    sed -i '' "s/${molNum}.*/${molName}_CSSPR/g" "$accSprCSSPRInp"
fi
# Prepare param file
echo "VERBOSITY: 1" > "$accSprCSSPRParFile"
echo "TASK: PrepareInputSpartan" >> "$accSprCSSPRParFile"
echo "INFILE: $accSprCSSPRInp" >> "$accSprCSSPRParFile"
echo "KEYWORDS: CONFANAL MMFF PRINTLEV=3 SEARCHMETHOD=SPARSE SPARSE=$numConformersSparse CONFSEED=$rnd CONFSKEPT=$maxSparseConformers CONF_SELECTION_RULE=4 USEFFPARMS=~/PARAMS.MMFF.DENOPTIM" >> "$accSprCSSPRParFile"
echo "CHARGE: $charge" >> "$accSprCSSPRParFile"
echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprCSSPRParFile"
echo "OUTFILE: $sprtCSSPRInp" >> "$accSprCSSPRParFile"
echo "\$STARTATOMSPECIFICKEYWORDS: $sprtCSAtmSpecKey" >> "$accSprCSSPRParFile"
echo "\$END" >> "$accSprCSSPRParFile"
#WARNING! causes  errors in linux version
#echo "\$STARTFREEZEATOMS: $sprtGOFrozen" >> "$accSprCSSPRParFile"
#echo "\$END" >> "$accSprCSSPRParFile"
echo "\$STARTCONSTRAINTS: $sprtConstrCS" >> "$accSprCSSPRParFile"
echo "\$END" >> "$accSprCSSPRParFile"
#
# Launch ACC
autocompchem -p "$accSprCSSPRParFile" > "$accSprCSSPRLog"
# Check outcome
if [ ! -f "$accSprCSSPRLog" ]; then
    errMsg="#MakeSpartanInput: $accSprCSSPRLog not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -d "$sprtCSSPRInp" ]; then
    errMsg="#MakeSpartanInput: $sprtCSSPRInp directory not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Termination status: 0" "$accSprCSSPRLog" ; then
    errMsg="#MakeSpartanInput: non-zero exit status from AutoCompChem."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi


#
# Submit SPARSE conformational search
#
echo "Submitting SPARSE Conformational Search..."
# Setting new parameters
sprtCSSPROutLog="$sprtCSSPRInp/${molName}_CSSPR/output"
sprtCSSPROutDir="$wrkDir/${molName}_CSSPR.Conf.${molName}_CSSPR.spartan"
sprtCSSPRLog="$wrkDir/${molName}_SprtCSSPRRun.log"
# Lauch Spartan
"$SPARTANEXE" --foreground-submit "$sprtCSSPRInp" > "$sprtCSSPRLog"
# Check outcome
if [ ! "$?" == 0 ]; then
    errMsg="#RunSpartan: non-zero exit status from Spartan call."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -f "$sprtCSSPROutLog" ]; then
    errMsg="#RunSpartan: $sprtCSSPROutLog file not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -d "$sprtCSSPROutDir" ]; then
    echo "WARNING! Conformational search with SPARSE did not produce anything, but we go on anyway."
    #errMsg="#RunSpartan: $sprtCSSPROutDir folder not found."
    #abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "Successful completion" "$sprtCSSPROutLog" ; then
    echo "WARNING! Conformational search with SPARSE did not complete succesfully, but we go on anyway."
    #errMsg="#RunSpartan: unsuccessful completion."
    #abandon "$DnCG3Dout" "$E_OPTERROR"
fi


#
# Run independent Monte Carlo searches
#
n=0
lowestEnergyConfSDF="$sprtGOOut"
lowestEnergy="$mmffGOEnergy"
if [ -d "$sprtCSSPROutDir" ]; then
    sprtCSSPRConfs=($(ls -d "$sprtCSSPROutDir/Conformer."*))
    echo "Running ${#sprtCSSPRConfs[@]} independent Monte Carlo conformational searches"
    for conformer in ${sprtCSSPRConfs[@]}
    do
        n=$((n+1))
        #
        # Extracting structure of conformer
        #
        echo "$n-Extracting geometry of conformer..."
        # Setting new params
        accSprCSSPROutParFile="$wrkDir/${molName}_SprtCSSPROut${n}.par"
        accSprCSSPROutLog="$wrkDir/${molName}_SprtCSSPROut${n}.log"
        sprtCSSPROut="$wrkDir/${molName}_SprtCSSPROut${n}.sdf"
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprCSSPROutParFile"
        echo "TASK: ExtractLastGeometryFromSpartanTree" >> "$accSprCSSPROutParFile"
        echo "INFILE: $sprtCSSPROutDir" >> "$accSprCSSPROutParFile"
        echo "MOLNAME: Conformer.${n}" >> "$accSprCSSPROutParFile"
        echo "OUTFORMAT: SDF" >> "$accSprCSSPROutParFile"
        echo "OUTFILE: $sprtCSSPROut" >> "$accSprCSSPROutParFile"
        # Launch ACC
        autocompchem -p "$accSprCSSPROutParFile" > "$accSprCSSPROutLog"
        # Check outcome
        if [ ! -f "$accSprCSSPROutLog" ]; then
            errMsg="#ExtractSpartanOutput: $accSprCSSPROutLog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtCSSPROut" ]; then
            errMsg="#ExtractSpartanOutput: $sprtCSSPROut file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprCSSPROutLog" ; then
            errMsg="#ExtractSpartanOutput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ "$SEDSYNTAX" == "GNU" ]
        then
            sed -i "s/Conformer.*/$molName/g" "$sprtCSSPROut"
        elif [ "$SEDSYNTAX" == "BSD" ]
        then
            sed -i '' "s/Conformer.*/$molName/g" "$sprtCSSPROut"
        fi
    
    
        #
        # Prepare input for MC conformational search
        #
        echo "$n-Preparing input for Monte Carlo search from conformer..."
        # Setting new params
        accSprCSMCNParFile="$wrkDir/${molName}_SprtCSMC${n}In.par"
        accSprCSMCNLog="$wrkDir/${molName}_SprtCSMC${n}In.log"
        sprtCSMCNInp="$wrkDir/${molName}_CSMC${n}.spartan"
        accSprCSMCNInp="$wrkDir/${molName}_SprtCSMC${n}In.sdf"
        rnd=$(echo $RANDOM)
        cp "$sprtCSSPROut" "$accSprCSMCNInp"
        if [ "$SEDSYNTAX" == "GNU" ]
        then
            sed -i "s/${molName}.*/${molName}_CSMC${n}/g" "$accSprCSMCNInp"
        elif [ "$SEDSYNTAX" == "BSD" ]
        then
            sed -i '' "s/${molNum}.*/${molName}_CSMC${n}/g" "$accSprCSMCNInp"
        fi
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprCSMCNParFile"
        echo "TASK: PrepareInputSpartan" >> "$accSprCSMCNParFile"
        echo "INFILE: $accSprCSMCNInp" >> "$accSprCSMCNParFile"
        #echo "KEYWORDS: SCONFANAL MMFF PRINTLEV=2 SEARCHMETHOD=MC CONFSEXAMINED=$numConformersMC CONFSEED=$rnd CONF_SELECTION_RULE=4 USEFFPARMS=~/PARAMS.MMFF.DENOPTIM" >> "$accSprCSMCNParFile"
        echo "KEYWORDS: SCONFANAL MMFF PRINTLEV=2 SEARCHMETHOD=MC CONFSEXAMINED=$numConformersMC CONFSEED=$rnd USEFFPARMS=~/PARAMS.MMFF.DENOPTIM" >> "$accSprCSMCNParFile"
        echo "CHARGE: $charge" >> "$accSprCSMCNParFile"
        echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprCSMCNParFile"
        echo "OUTFILE: $sprtCSMCNInp" >> "$accSprCSMCNParFile"
        echo "\$STARTATOMSPECIFICKEYWORDS: $sprtCSAtmSpecKey" >> "$accSprCSMCNParFile"
        echo "\$END" >> "$accSprCSMCNParFile"
        #WARNING! causes problems in Linux
        #echo "\$STARTFREEZEATOMS: $sprtGOFrozen" >> "$accSprCSMCNParFile"
        #echo "\$END" >> "$accSprCSMCNParFile"
        echo "\$STARTCONSTRAINTS: $sprtConstrCS" >> "$accSprCSMCNParFile"
        echo "\$END" >> "$accSprCSMCNParFile"
        echo "\$STARTROTATABLEBONDS: $sprtRotBnd" >> "$accSprCSMCNParFile"
        echo "\$END" >> "$accSprCSMCNParFile"
        # Launch ACC
        autocompchem -p "$accSprCSMCNParFile" > "$accSprCSMCNLog"
        # Check outcome
        if [ ! -f "$accSprCSMCNLog" ]; then
            errMsg="#MakeSpartanInput: $accSprCSMCNLog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -d "$sprtCSMCNInp" ]; then
            errMsg="#MakeSpartanInput: $sprtCSMCNInp directory not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprCSMCNLog" ; then
            errMsg="#MakeSpartanInput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        
        
        #
        # Submit conformational search
        #
        echo "$n-Submitting Spartan Conformational Search..."
        # Setting new parameters
        sprtCSMCNOutLog="$sprtCSMCNInp/${molName}_CSMC${n}/output"
        sprtCSMCNLog="$wrkDir/${molName}_SprtCSMC${n}Run.log"
        # Lauch Spartan
        "$SPARTANEXE" --foreground-submit "$sprtCSMCNInp" > "$sprtCSMCNLog"
        # Check outcome
        if [ ! "$?" == 0 ]; then
            errMsg="#RunSpartan: non-zero exit status from Spartan call."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtCSMCNOutLog" ]; then
            errMsg="#RunSpartan: $sprtCSMCNOutLog file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Successful completion" "$sprtCSMCNOutLog" ; then
            errMsg="#RunSpartan: unsuccessful completion."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        
        
        #
        # Extract lowest energy conformer
        #
        echo "$n-Extracting results from Spartan..."
        # Setting new params
        accSprCSMCNOutParFile="$wrkDir/${molName}_SprtCSMC${n}Out.par"
        accSprCSMCNOutLog="$wrkDir/${molName}_SprtCSMC${n}Out.log"
        sprtCSMCNOut="$wrkDir/${molName}_SprtCSMC${n}Out.sdf"
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprCSMCNOutParFile"
        echo "TASK: ExtractLastGeometryFromSpartanTree" >> "$accSprCSMCNOutParFile"
        echo "INFILE: $sprtCSMCNInp" >> "$accSprCSMCNOutParFile"
        echo "OUTFORMAT: SDF" >> "$accSprCSMCNOutParFile"
        echo "OUTFILE: $sprtCSMCNOut" >> "$accSprCSMCNOutParFile"
        # Launch ACC
        autocompchem -p "$accSprCSMCNOutParFile" > "$accSprCSMCNOutLog"
        # Check outcome
        if [ ! -f "$accSprCSMCNOutLog" ]; then
            errMsg="#ExtractSpartanOutput: $accSprCSMCNOutLog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtCSMCNOut" ]; then
            errMsg="#ExtractSpartanOutput: $sprtCSMCNOut file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprCSMCNOutLog" ; then
            errMsg="#ExtractSpartanOutput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ "$SEDSYNTAX" == "GNU" ]
        then
            sed -i "1 s/^.*$/$molName/g" "$sprtCSMCNOut"
        elif [ "$SEDSYNTAX" == "BSD" ]
        then
            sed -i '' "1 s/^.*$/$molName/g" "$sprtCSMCNOut"
        fi
        
        #
        # Extract energy
        #
        echo "$n-Extracting the MMFF energy of the selected conformation..."
        mmffCSMCNEnergy=$(grep "Lowest energy conformation" "$sprtCSMCNOutLog" | awk '{print $4}')
        mmffCSMCNEnergy=$(echo ${mmffCSMCNEnergy} | sed -e 's/[eE]+*/\\*10\\^/')
        mmffCSMCNEnergy=$(bc -l <<< "($mmffCSMCNEnergy)")
        
        
        #
        # Prepare input for single point energy calculation
        #
        echo "$n-Prepare input for single point MMFF energy calculation..."
        # Setting new params
        accSprSPNParFile="$wrkDir/${molName}_SprtSP${n}In.par"
        accSprSPNLog="$wrkDir/${molName}_SprtSP${n}In.log"
        accSprSPNInp="$wrkDir/${molName}_SprtSP${n}In.sdf"
        sprtSPNInp="$wrkDir/${molName}_SP${n}.spartan"
        cp "$sprtCSMCNOut" "$accSprSPNInp"
        if [ "$SEDSYNTAX" == "GNU" ]
        then
            sed -i "s/${molName}.*/${molName}_SP${n}/g" "$accSprSPNInp"
        elif [ "$SEDSYNTAX" == "BSD" ]
        then
            sed -i '' "s/${molNum}.*/${molName}_SP${n}/g" "$accSprSPNInp"
        fi
        # Prepare param file
        echo "VERBOSITY: 1" > "$accSprSPNParFile"
        echo "TASK: PrepareInputSpartan" >> "$accSprSPNParFile"
        echo "INFILE: $accSprSPNInp" >> "$accSprSPNParFile"
        echo "KEYWORDS: SOPT MMFF PRINTLEV=3 USEFFPARMS=~/PARAMS.MMFF.DENOPTIM" >> "$accSprSPNParFile"
        echo "CHARGE: $charge" >> "$accSprSPNParFile"
        echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprSPNParFile"
        echo "OUTFILE: $sprtSPNInp" >> "$accSprSPNParFile"
        echo "\$STARTATOMSPECIFICKEYWORDS: $sprtSPAtmSpecKey" >> "$accSprSPNParFile"
        echo "\$END" >> "$accSprSPNParFile"
        # Launch ACC
        autocompchem -p "$accSprSPNParFile" > "$accSprSPNLog"
        # Check outcome
        if [ ! -f "$accSprSPNLog" ]; then
            errMsg="#MakeSpartanSPNInput: $accSprSPNLog not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -d "$sprtSPNInp" ]; then
            errMsg="#MakeSpartanSPNInput: $sprtSPNInp directory not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Termination status: 0" "$accSprSPNLog" ; then
            errMsg="#MakeSpartanSPNInput: non-zero exit status from AutoCompChem."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        
        
        #
        # Submit Single point energy calculation
        #
        echo "$n-Submitting Spartan single point MMFF energy calculation..."
        # Setting new parameters
        sprtSPNOutLog="$sprtSPNInp/${molName}_SP${n}/output"
        sprtSPNLog="$wrkDir/${molName}_SprtSP${n}Run.log"
        # Lauch Spartan
        "$SPARTANEXE" --foreground-submit "$sprtSPNInp" > "$sprtSPNLog"
        # Check outcome
        if [ ! "$?" == 0 ]; then
            errMsg="#RunSpartanSPN: non-zero exit status from Spartan call."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if [ ! -f "$sprtSPNOutLog" ]; then
            errMsg="#RunSpartanSPN: $sprtSPNOutLog file not found."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        if ! grep -q "Successful completion" "$sprtSPNOutLog" ; then
            errMsg="#RunSpartanSPN: unsuccessful completion."
            abandon "$DnCG3Dout" "$E_OPTERROR"
        fi
        
        
        #
        # Extract energy
        #
        echo "$n-Extracting the SP-MMFF energy of the selected conformation..."
        mmffCSMCNEnergySP=$(grep " Totals " "$sprtSPNOutLog" | tail -n 1 | awk '{print $2}')
        mmffCSMCNEnergySP=$(echo ${mmffCSMCNEnergySP} | sed -e 's/[eE]+*/\*10\^/')
        mmffCSMCNEnergySP=$(bc -l <<< "($mmffCSMCNEnergySP)")
        echo "SP-MMFF Energy from CSMCN: $mmffCSMCNEnergySP (from CSMC$n $mmffCSMCNEnergy)"
        
        
        #
        # Possibly update lowest value
        #
        newLowest=$( bc -l <<< "$mmffCSMCNEnergySP < $lowestEnergy")
        if [ "$newLowest" == 1 ]
        then
            echo " New lowest MMFF energy conformation found!"
            lowestEnergy="$mmffCSMCNEnergySP"
            lowestEnergyConfSDF="$sprtCSMCNOut"
        fi
    done
fi
nmm="$n"


#
# Semiempirical refinement of the lowest MM energy so far 
#

#
# Prepare input for semi-empirical refinement
#
echo "Prepare input for first semi-empirical geometry optimization (GOSE0)..."
# Setting new params
accSprGOSEParFile="$wrkDir/${molName}_SprtGOSE0In.par"
accSprGOSELog="$wrkDir/${molName}_SprtGOSE0In.log"
accSprGOSEInp="$wrkDir/${molName}_SprtGOSE0In.sdf"
sprtGOSEInp="$wrkDir/${molName}_GOSE0.spartan"
cp "$lowestEnergyConfSDF" "$accSprGOSEInp"
if [ "$SEDSYNTAX" == "GNU" ]
then
    sed -i "s/${molName}.*/${molName}_GOSE0/g" "$accSprGOSEInp"
elif [ "$SEDSYNTAX" == "BSD" ]
then
    sed -i '' "s/${molNum}.*/${molName}_GOSE0/g" "$accSprGOSEInp"
fi
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
if [ "$SEDSYNTAX" == "GNU" ]
then
    sed -i "1 s/^.*$/$molName/g" "$sprtGOSEOut"
elif [ "$SEDSYNTAX" == "BSD" ]
then
    sed -i '' "1 s/^.*$/$molName/g" "$sprtGOSEOut"
fi


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
if [ "$SEDSYNTAX" == "GNU" ]
then
    sed -i "s/${molName}.*/${molName}_DYN/g" "$accSprDYNInp"
elif [ "$SEDSYNTAX" == "BSD" ]
then
    sed -i '' "s/${molNum}.*/${molName}_DYN/g" "$accSprDYNInp"
fi
# Prepare param file
echo "VERBOSITY: 1" > "$accSprDYNParFile"
echo "TASK: PrepareInputSpartan" >> "$accSprDYNParFile"
echo "INFILE: $accSprDYNInp" >> "$accSprDYNParFile"
echo "KEYWORDS: DYNCON MMFF PRINTLEV=2 DYNCONMETHOD=GRID OPTCYCLE=2600 USEFFPARMS=~/PARAMS.MMFF.DENOPTIM" >> "$accSprDYNParFile"
#
echo "CHARGE: $charge" >> "$accSprDYNParFile"
echo "SPIN_MULTIPLICITY: $spinMult" >> "$accSprDYNParFile"
echo "OUTFILE: $sprtDYNInp" >> "$accSprDYNParFile"
echo "\$STARTATOMSPECIFICKEYWORDS: $sprtCSAtmSpecKey" >> "$accSprDYNParFile"
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
# Need to appen the preference or Spartan will crash
echo "BEGINPREFERENCES" >> "$sprtDYNInp/${molName}_DYN/input"
echo " MM:CONF_SELECTION_RULE=2" >> "$sprtDYNInp/${molName}_DYN/input"
echo "ENDPREFERENCES" >> "$sprtDYNInp/${molName}_DYN/input"


#
# Submit DYNamic constraint
#
echo "Submitting DYNCON generation of alternative conformation..."
# Setting new parameters
sprtDYNOutLog="$sprtDYNInp/${molName}_DYN/output"
sprtDYNOutDir="$wrkDir/${molName}_DYN.Prof.${molName}_DYN.spartan"
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
    if ! grep -q "No Dynamic Constraints Found" "$sprtDYNOutLog" ; then
        errMsg="#RunSpartan: unsuccessful completion."
        abandon "$DnCG3Dout" "$E_OPTERROR"
    fi
fi

#
# Refine the alternative conformers
#
if [ -d "$sprtDYNOutDir"* ]
then
    n=0
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
        if [ "$SEDSYNTAX" == "GNU" ]
        then
            sed -i "s/Profile.*/$molName/g" "$sprtDYNOut"
        elif [ "$SEDSYNTAX" == "BSD" ]
        then
            sed -i '' "s/Profile.*/$molName/g" "$sprtDYNOut"
        fi

        
        #
        # Prepare input semi-empirical refinement
        #
        echo "$n-Prepare input for semi-empirical geometry optimization..."
        # Setting new params
        accSprGOSEParFile="$wrkDir/${molName}_SprtGOSE${n}In.par"
        accSprGOSELog="$wrkDir/${molName}_SprtGOSE${n}In.log"
        accSprGOSEInp="$wrkDir/${molName}_SprtGOSE${n}In.sdf"
        sprtGOSEInp="$wrkDir/${molName}_GOSE${n}.spartan"
        cp "$sprtDYNOut" "$accSprGOSEInp"
        if [ "$SEDSYNTAX" == "GNU" ]
        then
            sed -i "s/${molName}.*/${molName}_GOSE${n}/g" "$accSprGOSEInp"
        elif [ "$SEDSYNTAX" == "BSD" ]
        then
            sed -i '' "s/${molNum}.*/${molName}_GOSE${n}/g" "$accSprGOSEInp"
        fi
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
        if [ "$SEDSYNTAX" == "GNU" ]
        then
            sed -i "1 s/^.*$/$molName/g" "$sprtGOSEOut"
        elif [ "$SEDSYNTAX" == "BSD" ]
        then
            sed -i '' "1 s/^.*$/$molName/g" "$sprtGOSEOut"
        fi
       
     
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
else
    echo "No alternative rotamer for the specified bonds."
fi


#
# Prepare output
#
cp "$lowestEnergyConfSDF" "$outSDF"
addPropertyToSingleMolSDF "FrozenCore-PM6_ENERGY" "$lowestEnergy" "$outSDF"


#
# Calculate RMSD between conformations
#
echo "Calculating RMSD due to conformational search and refinement..."
# Setting new parameters
cgConfTnk="$wrkDir/${molName}_cg.xyz"
cgConfKey="$wrkDir/${molName}_cg.key"
sprtCSTnk="$wrkDir/${molName}_sprt.xyz"
sprtCSKey="$wrkDir/${molName}_sprt.key"
rmdsLogCGCS="$wrkDir/${molName}_rmsdConf.log"
# prepare input for Tinker's superpose
obabel -isdf "$DnCG3Dout" -otxyz -O "$cgConfTnk"
echo "parameters $tinkerForceField" > "$cgConfKey"
obabel -isdf "$lowestEnergyConfSDF" -otxyz -O "$sprtCSTnk"
echo "parameters $tinkerForceField" > "$sprtCSKey"
# Launch Tinker
"$pathToTinkerBin/superpose" "$sprtCSTnk" "$cgConfTnk" 1 Y U Y 0.0 > "$rmdsLogCGCS"
# Check outcome
if [ ! "$?" == 0 ]; then
    errMsg="#RMSD: non-zero exit status from RMSD calculation."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -f "$rmdsLogCGCS" ]; then
    errMsg="#RMSD: $rmdsLogCGCS file not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "" "$rmdsLogCGCS" ; then
    errMsg="#RMSD: unsuccessful completion."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
rmsdConfCGCS=$(grep "Root Mean Square Distance" "$rmdsLogCGCS" | tail -n 1 | awk '{print $6}')
addPropertyToSingleMolSDF "RMSD_CONF_CG-CS" "$rmsdConfCGCS" "$outSDF"

# Setting new parameters
sprtGOTnk="$wrkDir/${molName}_go.xyz"
sprtGOKey="$wrkDir/${molName}_go.key"
rmdsLogGOCS="$wrkDir/${molName}_rmsdConf.log"
# prepare input for Tinker's superpose
obabel -isdf "$lowestEnergyConfSDF" -otxyz -O "$sprtGOTnk"
echo "parameters $tinkerForceField" > "$sprtGOKey"
# Launch Tinker
"$pathToTinkerBin/superpose" "$sprtCSTnk" "$sprtGOTnk" 1 Y U Y 0.0 > "$rmdsLogGOCS"
# Check outcome
if [ ! "$?" == 0 ]; then
    errMsg="#RMSD: non-zero exit status from RMSD calculation."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if [ ! -f "$rmdsLogGOCS" ]; then
    errMsg="#RMSD: $rmdsLogGOCS file not found."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
if ! grep -q "" "$rmdsLogGOCS" ; then
    errMsg="#RMSD: unsuccessful completion."
    abandon "$DnCG3Dout" "$E_OPTERROR"
fi
rmsdConfGOCS=$(grep "Root Mean Square Distance" "$rmdsLogGOCS" | tail -n 1 | awk '{print $6}')
addPropertyToSingleMolSDF "RMSD_CONF_GO-CS" "$rmsdConfGOCS" "$outSDF"


#
# Free-up token: done by function 'finish'
#
#releaseLock "$lockIDFile" "$lockFile"


#
# Cleanup
#
if [ "$cleanup" == 0 ]
then
    echo "Cleaning up..."
    rm -f "$gm3dfInp"
    rm -f "$gm3dfPar"
    rm -f "$gm3dfLog"
    rm -f "$DnCG3Dinp"
    rm -f "$DnCG3Dout"
    rm -f "$DnCGParFile"
    rm -f "$atmClshParFile"
    rm -f "$atmClshLog"
    rm -f "$accSprGOParFile"
    rm -f "$accSprGOLog"
    rm -f "$accSprGOInp"
    rm -f "$sprtGOLog"
    rm -f "$accSprGOOutParFile"
    rm -f "$accSprGOOutLog"
    rm -f "$sprtGOOut"
    rm -f "$accSprSP0ParFile"
    rm -f "$accSprSP0Log"
    rm -f "$accSprSP0Inp"
    rm -f "$sprtSP0Log"
    rm -r -f "$sprtSP0Inp"
    rm -f "$accSprCSSPRParFile"
    rm -f "$accSprCSSPRLog"
    rm -f "$accSprCSSPRInp"
    rm -f "$sprtCSSPRLog"
#    rm -f "$lowestEnergyConfSDF"
    rm -f "$accSprCSSPROutParFile"
    rm -f "$accSprCSSPROutLog"
    rm -f "$sprtCSSPROut"
    rm -f "$accSprCSMCNParFile"
    rm -f "$accSprCSMCNLog"
    rm -f "$accSprCSMCNInp"
    rm -f "$sprtCSMCNLog"
    rm -f "$accSprCSMCNOutParFile"
    rm -f "$accSprCSMCNOutLog"
    rm -f "$sprtCSMCNOut"
    rm -f "$accSprSPNParFile"
    rm -f "$accSprSPNLog"
    rm -f "$accSprSPNInp"
    rm -f "$sprtSPNLog"
    rm -f "$accSprDYNParFile"
    rm -f "$accSprDYNLog"
    rm -r -f "$sprtDYNInp"
    rm -f "$accSprDYNInp"
#    rm -r -f "$sprtGOSEInp"
    rm -f "$sprtDYNOutLog"
    rm -r -f "$sprtDYNOutDir"
    rm -f "$sprtDYNLog"
    rm -f "$cgConfTnk"
    rm -f "$cgConfKey"
    rm -f "${cgConfTnk}_2"
    rm -f "$rmdsLogCGCS"
    rm -f "$sprtGOTnk"
    rm -f "${sprtGOTnk}_2"
    rm -f "$sprtGOKey"
    rm -f "$rmdsLogGOCS"
    rm -f "$sprtCSTnk"
    rm -f "$sprtCSKey"
    rm -rf "$sprtGOInp"
    rm -rf "$sprtCSSPRInp"
    #NEEDED! rm -rf "$sprtCSSPROutDir" 
    rm -rf "$sprtCSMCNInp"
    for i in $(seq 0 $((nmm+1)))
    do
        rm -f "$wrkDir/${molName}_SprtCSSPROut${i}.par"
        rm -f "$wrkDir/${molName}_SprtCSSPROut${i}.log"
        rm -f "$wrkDir/${molName}_SprtCSSPROut${i}.sdf"
        rm -f "$wrkDir/${molName}_SprtCSMC${i}In.par"
        rm -f "$wrkDir/${molName}_SprtCSMC${i}In.log"
        rm -f "$wrkDir/${molName}_SprtCSMC${i}In.sdf"
        rm -f "$wrkDir/${molName}_SprtCSMC${i}Run.log"
        rm -f "$wrkDir/${molName}_SprtCSMC${i}Out.par"
        rm -f "$wrkDir/${molName}_SprtCSMC${i}Out.log"
        rm -f "$wrkDir/${molName}_SprtCSMC${i}Out.sdf"
        rm -rf "$wrkDir/${molName}_CSMC${i}.spartan"
        rm -f "$wrkDir/${molName}_SprtSP${i}In.par"
        rm -f "$wrkDir/${molName}_SprtSP${i}In.log"
        rm -f "$wrkDir/${molName}_SprtSP${i}In.sdf"
        rm -rf "$wrkDir/${molName}_SP${i}.spartan"
        rm -f "$wrkDir/${molName}_SprtSP${i}Run.log"
    done
    for i in $(seq 0 $((n+1)))
    do
        rm -f "$wrkDir/${molName}_SprtDYNOut${i}.par"
        rm -f "$wrkDir/${molName}_SprtDYNOut${i}.log"
        rm -f "$wrkDir/${molName}_SprtDYNOut${i}.sdf"
        rm -f "$wrkDir/${molName}_SprtGOSE${i}In.par"
        rm -f "$wrkDir/${molName}_SprtGOSE${i}In.log"
        rm -f "$wrkDir/${molName}_SprtGOSE${i}In.sdf"
        rm -f "$wrkDir/${molName}_SprtGOSE${i}Run.log"
        rm -f "$wrkDir/${molName}_SprtGOSE${i}Out.par"
        rm -f "$wrkDir/${molName}_SprtGOSE${i}Out.log"
        rm -fr "$wrkDir/${molName}_GOSE${i}.spartan"
        rm -f "$wrkDir/${molName}_SprtGOSE${i}Out.sdf"
    done
    echo "Done!"
fi

#
# Create task completed label
#
writeTCL "$jobId"

exit 0
