#!/bin/bash
###############################################################################
#
# This is the fitness evaluation script for Ru-catalysts for olefin metathesis.
# It computed the fitness by elaborating free energy values for
# the precursor, free L-ligand, TS-guesses for catalysis and decomposition.
#
# @param $1 pathname of input SDF file: it must contain the graph 
#        representation of the candidate system of which we are calculating 
#        the fitness.
# @param $2 pathname of the output SDF file where to report the fitness 
#        (possibly with some additional information attached) or the error/s
#        justifying rejecting of this candidate.
# @param $3 pathname of the working space, i.e., a directory.
# @param $4 numerical ID of the task this script is asked to perform.
# @param $5 pathname to the file collecting the candidates unique identifiers 
#        collected up to the moment when this script is asked to evaluate a new
#        candidate
#
# @author Marco Foscato
# @author Jonas Brattebø Ekeli
#
###############################################################################
# Tunable Parameters
###############################################################################

#FakeFitness: assign a fake fitness value doing the least work possible
# This is only used for testing and debugging purposes
useFakeFitness=1 # 0:doit, 1:don't
#Fast test
testRun=1 #0:this is a test, 1:this is NOT a test

#Path to state-specific scripts
build3DConfScript="$WORKDIR/run3DBuildAndConfSearch-4.0.sh"
build3DAndRefine="$WORKDIR/run3DBuildAndRefine-5.0.sh"
runXTBScript="$WORKDIR/runXTBTask-1.0.sh"
runManyXTBScript="$WORKDIR/runXTBMultiTask-2.0.sh"
runDFTScript="$WORKDIR/runDFTTask-2.0.sh"
runManyDFTScript="$WORKDIR/runDFTManyTask-1.0.sh"
runPseudoTSDFTPreOpt="$WORKDIR/runDFTOptPseudoTS-1.0.sh"

#Searh for (MOLUID), and import candidates that have already been evaluated.
importOld="1" # "0" = do it, "1" = don't
#Paths to old RUN#### folders to search for previously completed candidates.
pathToOld=""

#Label and extension of input file
inputLabelExt="_inp"

# Define labels of alternative sub-jobs with extensive conformational search
# that are used to select the ligand conformation. Only the results of the
# best of these sub-jobs will be kept and used to define the conformation of
# the ligand used for other intermediates/transition states.
# LABEL SYNTAX: 
# -> the first character is the label root, and identifies the species
# -> the rest is used for alternative models of the same species
# For example, labels A1, A2, and A3 they all pertain to the same
# species and only one of the alternatives models will be eventually selected
# for further modeling.
#
firstSubJobsLabelList=("F1" "F2" "A1" "A2" "L0")
#
labelForLigandExtraxtion="A"
#Define which label in the list above is screened to find the lowest energy from
# result from which the ligand conformatio is extracted
#
# Defines the labels of alternative sub-jobs for the states to be modeled 
# using the ligand's conformation defined from first batch.
secondSubJobsLabelList=("C1" "C2" "E1" "E2" "X1" "X2" "Z1" "Z2")
#
# A1: Hoveyda-Grubbs precursor any X-ligand - 1st stereoisomer
# A2: Hoveyda-Grubbs precursor any X-ligand - 2nd stereoisomer
# F1: Hoveyda-Grubbs precursor cis X-ligand - 1nd stereoisomer
# F2: Hoveyda-Grubbs precursor cis X-ligand - 2nd stereoisomer
# E1: Hoveyda-Grubbs precursor with Cl - 1st stereoisomer
# E2: Hoveyda-Grubbs precursor with Cl - 2nd stereoisomer
# C1: Ru(IV) MCB from productive metathesis - 1st stereoisomer
# C2: Ru(IV) MCB from productive metathesis - 2nd stereoisomer
# X1: TS-guess for productive metathesis - 1st stereoisomer
# X2: TS-guess for productive metathesis - 2nd stereoisomer
# Z1: TS-guess for internal beta-H elimination - 1st stereoisomer
# Z2: TS-guess for internal beta-H elimination - 2st stereoisomer
# L0: free ligand
#
# WARNING! The libraries of scaffolds and other state-dependent parameters are
# set by the run3DBuildAnd******-*.*.sh script based on the label given to the
# sub-job.
#

# Labels for pre DFT geom opt.
preDFTOptLabels=("X" "Z" "C")
#
preDFTtoSPLabels=("X" "Z" "D")
# Labels corresponding to transition states
transitionStateLabels=("X1" "X2" "Z1" "Z2")
#
# Sort the label roots (i.e., the first letter of the labels) in order of 
# challenging calculation (the worst first, the easiest last):
sortedLabelRoots=("A" "F" "E" "C" "X" "Z" "L")
#
# Labels of states that are submitted to XTB Opt
toXTBStateLabels=("A" "E" "C" "L" "F")
#
# Labels of states that are submitted to DFT (D is the DFT optimized version of C
# and comes from the same conformational search as C)
toDFTStateLabels=("A" "E" "C" "L" "F" "X" "Z" "D")
# Labels of states that MUST be successful for the fitness to be completed
#
requiredStateLabels=("A" "E" "C" "X" "Z" "L" "F")

#
# WARNING! Check the part of the script that calculated the delta_G: loss of
# generality there, and assumption based on the current relation between
# molecules and labels.
#

#Wait for conformational search
csMinTime=1       # delay of first checking iteration
csMinTimeUnit="m" # time unit: s (seconds), m (minutes), h (hours)
csStep=2          # delay of each checking iteration, but the first
csTimeUnit="m"    # time unit: s (seconds), m (minutes), h (hours)
csMaxWait=5000     # maximum number of checking iterations

# Fragment space
compMatrix="$WORKDIR/CPMap_FS-C5.0.par" 
libFrags="$WORKDIR/lib_frags_FS_C5.0_v3.sdf"
libCaps="$WORKDIR/lib_caps_FS_C5.0_v3.sdf"
libScaff="$WORKDIR/scaff_du_v3.sdf"

# Graph editing tasks 
graphEditingTasks="$WORKDIR/collapsLigand"

#Fragmentation
coreLigandCutRules="$WORKDIR/core-lig_CutRule.rul"
coreXLigandCutRules="$WORKDIR/core-Xlig_CutRule.rul"

#Graph for building of reaction steps
coreLigGraphTmpl="$WORKDIR/core-lig_Graph.sdf"

#Back door to hold new executions of scripts: create this file to make any 
#new submission of this script hang in "hold status" until the file is removed
flagHoldFitnessProvider="/tmp/hold_fitness_provider"
listHeldPIDs="/tmp/activeConfSearchOrRefine"
holdMaxWait=10000  # we will wait for up to this number of steps
holdStep=2         # length of waiting step
holdTimeUnit="m"   # time unit: s (seconds), m (minutes), h (hours)

#Wait for tasks run on HPC
stepFirstWaitForHPC=300       # delay of first checking iteration
stepFirstWaitForHPCXTB=30
unitFirstWaitForHPC="s" # time unit: s (seconds), m (minutes), h (hours)
stepWaitForHPC=30          # delay of each checking iteration, but the first
stepWaitForHPCXTB=2
unitWaitForHPC="m"    # time unit: s (seconds), m (minutes), h (hours)
maxWaitForHPC=4500    # maximum number of waiting cycles
if [ "$HIGHFREQUENCY" == 1 ]; then
  echo "(Checking for completion of external tasks with high frequncy)"
  csMinTime=2       # delay of first checking iteration
  csMinTimeUnit="s" # time unit: s (seconds), m (minutes), h (hours)
  csStep=2          # delay of each checking iteration, but the first
  csTimeUnit="s"    # time unit: s (seconds), m (minutes), h (hours)
  stepFirstWaitForHPC=2
  stepFirstWaitForHPCXTB=2
  unitFirstWaitForHPC="s" # time unit: s (seconds), m (minutes), h (hours)
  stepWaitForHPC=2          # delay of each checking iteration, but the first
  stepWaitForHPCXTB=2
  unitWaitForHPC="s"    # time unit: s (seconds), m (minutes), h (hours)
fi

# Back door to recover partial results from previous runs: use the optional
# parameter that gives a pathname to a file from which the information is taken
# and stored in these vectors:
isRestartRun=1    # 0: yes 1: no  => this is set by readPrevDataSettings
labelsCompletedStates=()
pathCompletedStatesDir=()
labelsTakeLowLevStates=()
labelsHighFreqStates=()

#Cleanup 0:doit, 1:don't
cleanup=0

#Exit code for incomplete evaluation of fitness
E_OPTERROR=0 # '0' leads to rejection of this single molecule
E_FATAL=1    # '1' stops the master DEMOPTIM job

#Initialization of job management variables
beginTime=$(date +%s)
errMsg="Error not assigned."


# Alter settings for test and debug mode
if [ "$testRun" == 0 ]
then
    csMinTime=10       # delay of first checking iteration
    csMinTimeUnit="s" # time unit: s (seconds), m (minutes), h (hours)
    csStep=10          # delay of each checking iteration, but the first
    csTimeUnit="s"    # time unit: s (seconds), m (minutes), h (hours)
    csMaxWait=5000     # maximum number of checking iterations 
    stepFirstWaitForHPC=15       # delay of first checking iteration
    unitFirstWaitForHPC="s" # time unit: s (seconds), m (minutes), h (hours)
    stepWaitForHPC=30          # delay of each checking iteration, but the first
    unitWaitForHPC="s"    # time unit: s (seconds), m (minutes), h (hours)
    stepWaitForHPCXTB="120" 
    maxWaitForHPC=1000    # manimum number of waiting cycles
    cleanup=1             # if 0 then remove tmp files
fi
if [ "$useFakeFitness" == 0 ]
then
    csMinTime=30        # delay of first checking iteration
    csMinTimeUnit="s"   # time unit: s (seconds), m (minutes), h (hours)
    csStep=1            # delay of each checking iteration, but the first
    csTimeUnit="m"      # time unit: s (seconds), m (minutes), h (hours)
    csMaxWait=100       # maximum number of checking iterations
    cleanup=1           # 1 means we don't remove tmp files 
fi

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
# Function that reads the state labels and pathnames needed to recover partial
# results from previous runs
# @param $1 the pathname of the text file to read
#

function readPrevDataSettings(){
    if [ ! -f $1 ]
    then
        echo " ERROR! File '$1' not found! "
        exit 1
    fi 
    while read -r line
    do 
        if [[ $line == "#"* ]] || [[ $line == "" ]]
        then
            continue
        fi
        oldIFS="$IFS"
        IFS=' ' read -ra words <<< "$line"
        IFS="$oldIFS"
        case ${words[0]} in
            "STATECOMPLETED" )
            if [ ${#words[@]} -ne 3 ]
            then
                echo "ERROR! Wrong syntax in line '$line'"
                exit 1
            fi
            isRestartRun=0
            labelsCompletedStates+=(${words[1]})
            pathCompletedStatesDir+=(${words[2]})
            ;;
            "STATEDFTHIGHFREQ" )
            if [ ${#words[@]} -ne 2 ]
            then
                echo "ERROR! Wrong syntax in line '$line'"
                exit 1
            fi
            labelsHighFreqStates+=(${words[1]})
            ;;
            *)
            echo " Keyword '${words[0]}' not recognized! Check file '$1'"
            exit 1
            ;;
        esac
    done < $1

    #check consistency
    firstPath="${pathCompletedStatesDir[0]}"
    for locPath in "${pathCompletedStatesDir[@]}"
    do
        if [[ "$locPath" != "$firstPath" ]]
        then
            echo "ERROR! Inconsistent pathnames of previous run folder."
            echo "$locPath != $firstPath"
            exit 1
        fi
    done
}


#
# Function to submit a batch of parallel MM+SemiEmp sub-jobs and wait for 
# completion.
# @param $1 batch ID
# @param $2 vector of labels, one per each sub-jobs. Must be given as vector[@]
#           with no double quotes
# @param $3 pathname of the BASH script to be called for this task
#

function submitParallelBatchMMSE() {
    batchID="$1"
    locInpSDF="$2"
    locSubJobLabelList=("${!3}")
    subJobScript="$4"
        
    # Submit all parallel sub-jobs
    echo "Submitting batch ${batchID}: ${locSubJobLabelList[@]}"
    locSubJobIDList=()
    locSubJobOutList=()
    for sjLab in "${locSubJobLabelList[@]}"
    do
        sjId="$taskId.$sjLab"
        # WARNING the "_outSub" string is assumed in downstread scripts
        sjOut="$wrkDir/${molNum}_outSub-$sjLab.sdf"
        submitSubJob "$locInpSDF" "$sjOut" "$wrkDir" "$subJobScript" "$sjId"
        locSubJobIDList+=("$sjId")
        locSubJobOutList+=("$sjOut")
        allSubJobsLabelList+=("$sjLab")
        allSubJobsIDList+=("$sjId")
        allSubJobsOutList+=("$sjOut")
    done
    
    # Wait for completion of parallel conformational searches
    echo "Start waiting for completion of submitted jobs (batch ${batchID})..."
    allDone=0
    for i in $(seq 1 $csMaxWait)
    do
        if [ $i == 1 ]
        then
            sleep $csMinTime$csMinTimeUnit
        else
            sleep $csStep$csTimeUnit
        fi
        dateTime=$(date)
        running=0
        for j in $(seq 0 $((${#locSubJobLabelList[@]}-1)) )
        do
            tcl="$wrkDir/tcl-${locSubJobIDList[$j]}"
            if [ ! -f "$tcl" ]
            then
                running=$((running+1))
                echo "$dateTime Still waiting for $tcl"
                if [ $i == $csMaxWait ]
                then
                    errMsg="#WaitingParallelBatch: time limit reached (task abandoned)"
                    abandon "$locInpSDF" "$E_OPTERROR"
                fi
            else
                if [ $j == $((${#locSubJobLabelList[@]}-1)) ] && [ $running == 0 ]
                then
                    allDone=1
                    break;
                fi
            fi
        done
        if [ $allDone == 1 ]
        then
            echo "$dateTime All done: stop waiting."
            break;
        fi 
    done
    
    # Check outcome, and collect energies from submitted tasks
    echo "Checking outcome of sub-jobs (batch ${batchID})..."
    energiesSubjobs=()
    for j in $(seq 0 $((${#locSubJobLabelList[@]}-1)) )
    do
        # Check output is found and not containing error
        out="${locSubJobOutList[$j]}"
        if [ ! -f "$out" ]
        then
            errMsg="#RetrieveConf: output file $out not found."
            abandon "$locInpSDF" "$E_OPTERROR"
        else
            # NB: here one could want to keep running if some, specific subjob
            # does not terminate properly.
            if grep -q "MOL_ERROR" "$out" > /dev/null
            then
                abandonDueToChild "$out"
            fi
        fi
       
        energy=$(grep -A1 FrozenCore-PM6_ENERGY "$out" | tail -n 1)
        energiesSubjobs+=("$energy")
     
        mv "$wrkDir/${locSubJobIDList[$j]}.log" "$wrkDir/${molNum}_MMSE${batchID}_${locSubJobIDList[$j]}.log"
        if [ "$cleanup" == 0 ]
        then
            rm -f "$wrkDir/tcl-${locSubJobIDList[$j]}"
            rm -f "$wrkDir/${locSubJobIDList[$j]}.nho"
            #rm -f "$wrkDir/${locSubJobIDList[$j]}.log"
        fi
    done
    
    # Selection among alternatives
    echo "Identification of lowest energy points (batch ${batchID})..."
    absMinValue=10000000000.0
    absMinLabel=""
    doneLabels=()
    for i in $(seq 0 $((${#locSubJobLabelList[@]}-1)) )
    do
        labI="${locSubJobLabelList[$i]}"
        rootI="${labI:0:1}"
        skip=1
        for doneLab in "${doneLabels[@]}"
        do
            if [ "$doneLab" == "$labI" ] 
            then
                skip=0
                break;
            fi
        done
        if [ "$skip" == 0 ]
        then
            continue
        fi
        groupValues=()
        groupLabels=()
        doneLabels+=("$labI")
        groupLabels+=("$labI")
        groupValues+=("${energiesSubjobs[$i]}")
        for j in $(seq $(($i+1)) $((${#locSubJobLabelList[@]}-1)) )
        do
            labJ="${locSubJobLabelList[$j]}"
            for doneLab in "${doneLabels[@]}"
            do
                if [ "$doneLab" == "$labJ" ]
                then
                    skip=0
                    break;
                fi
            done
            if [ "$skip" == 0 ]
            then
                continue
            fi
            rootJ="${labJ:0:1}"
            if [ "$rootI" == "$rootJ" ]
            then
                doneLabels+=("$labJ")
                groupLabels+=("$labJ")
                groupValues+=("${energiesSubjobs[$j]}")
            fi
        done
    
        echo "For species $rootI:"
        minValue=10000000000000.0
        minIdx=-1
        for j in $(seq 0 $((${#groupValues[@]}-1)) )
        do
            candidate="${groupValues[$j]}"
            echo " ${groupLabels[$j]} ($candidate)"
            # Evaluation within the candidate group
            newMin=$( bc -l <<< "$candidate < $minValue")
            if [ "$newMin" == 1 ]
            then
                minValue="$candidate"
                minIdx=$j
            fi
	    if [ "$rootI" == "$labelForLigandExtraxtion" ]
	    then
                # Evaluation in the entire block of sub-jobs
                newAbsMin=$( bc -l <<< "$candidate < $absMinValue")
                if [ "$newAbsMin" == 1 ]
                then
                    absMinValue="$candidate"
                    absMinLabel="${groupLabels[$j]}"
	        fi
            fi
        done
        selectedValues+=("$minValue")
        selectedLabels+=("${groupLabels[$minIdx]}")
        echo " Lowest energy among ${rootI}'s is ${groupLabels[$minIdx]} ($minValue)"
    done
    echo "====================================================================="
    echo "Lowest energy (batch $batchID) from $absMinLabel ($absMinValue)"
    echo "====================================================================="
}


#
# Function to submit a batch of parallel XTB jobs and wait for completion. No
# argument required: this function uses global variables as input. Namely:
# 'allSdfToXTB' and 'toDFTStateLabels'
#
# Many variable names refer to DFT because, still...
#

function submitParallelBatchXTB(){

    # Submit all parallel sub-jobs
    echo "Submitting batch of XTB jobs: ${toDFTStateLabels[@]}"
    locSubJobXTBIDList=()
    locSubJobXTBOutList=()
    for idft in $(seq 0 $((${#toDFTStateLabels[@]}-1)) )
    do
        sjLab="${toDFTStateLabels[idft]}"
        sjId="$taskId.${sjLab}#"
        # WARNING the "_outSub" string is assumed in downstread scripts
        sjOut="$wrkDir/${molNum}_outSubXTB-$sjLab.sdf"
        locInpSDF="${allSdfToXTB[idft]}"
        submitSubJob "$locInpSDF" "$sjOut" "$wrkDir" "$runXTBScript" "$sjId"
        locSubJobXTBIDList+=("$sjId")
        locSubJobXTBOutList+=("$sjOut")
        allXTBSubJobsLabelList+=("$sjLab")
        allXTBSubJobsIDList+=("$sjId")
        allXTBSubJobsOutList+=("$sjOut")
    done

    # Wait for completion of parallel jobs
    echo "Start waiting for completion of XTB jobs..."
    allXTBDone=0
    for i in $(seq 1 "$maxWaitForHPC")
    do
        if [ "$i" == $maxWaitForHPC ]
        then
            sleep "$stepFirstWaitForHPCXTB$unitFirstWaitForHPC"
        else
            sleep "$stepWaitForHPCXTB$unitWaitForHPC"
        fi
        dateTime=$(date)
        running=0
        for j in $(seq 0 $((${#toDFTStateLabels[@]}-1)) )
        do
            tcl="$wrkDir/tcl-${locSubJobXTBIDList[$j]}"
            if [ ! -f "$tcl" ]
            then
                running=$((running+1))
                echo "$dateTime Still waiting for $tcl"
                if [ "$i" == $maxWaitForHPC ]
                then
                    errMsg="#WaitingXTB: time limit reached (task abandoned)"
                    abandon "$locInpSDF" "$E_OPTERROR"
                fi
            else
                if [ "$j" == $((${#toDFTStateLabels[@]}-1)) ] && [ "$running" == 0 ]
                then
                    allXTBDone=1
                    break;
                fi
            fi
        done
        if [ "$allXTBDone" == 1 ]
        then
            echo "$dateTime All done: stop waiting."
            break;
        fi
    done

    # Check outcome, and collect results
    echo "Checking output files of XTB jobs..."
    for j in $(seq 0 $((${#toDFTStateLabels[@]}-1)) )
    do
        # Check output is found and not containing error
        out="${locSubJobXTBOutList[$j]}"
        if [ ! -f "$out" ]
        then
            errMsg="#RetrieveXTBResults: output file $out not found."
            abandon "$locInpSDF" "$E_OPTERROR"
        else
            # NB: here one could want to keep running if some, specific subjob
            # does not terminate properly.
            if grep -q "MOL_ERROR" "$out" > /dev/null
            then
                abandonDueToChild "$out"
            fi
        fi

        mv "$wrkDir/${locSubJobXTBIDList[$j]}.log" "$wrkDir/${molNum}_XTB_${locSubJobXTBIDList[$j]}.log"
        if [ "$cleanup" == 0 ]
        then
            rm -f "$wrkDir/tcl-${locSubJobXTBIDList[$j]}"
            rm -f "$wrkDir/${locSubJobXTBIDList[$j]}.nho"
        fi
    done
}

#
# Function to submit a batch of XTB jobs to be executed sequentially by a single
# child job, and wait for completion. No
# argument required: this function uses global variables as input. Namely:
# 'allSdfToXTB' and 'toDFTStateLabels'
#
# Many variable names refer to DFT because the functionality was originally meant
# for DFT, but your should ignore this labeling.
#

function submitXTBAllInOnce(){
    # Submit all parallel sub-jobs
    echo "Submitting batch of XTB jobs: ${toXTBStateLabels[@]}"
    locSubJobXTBOutList=()
    inpFilesList=""
    outFilesList=""
    for idft in $(seq 0 $((${#toXTBStateLabels[@]}-1)) )
    do
        sjLab="${toXTBStateLabels[idft]}"
        # WARNING the "_outSub" string is assumed in downstream scripts
        sjOut="$wrkDir/${molNum}_outSubXTB-$sjLab.sdf"
        inpFilesList="$inpFilesList,${allSdfToXTB[idft]}"
        outFilesList="$outFilesList,$sjOut"
        allXTBSubJobsOutList+=("$sjOut")
        locSubJobXTBOutList+=("$sjOut")
    done

    tcl="$wrkDir/tcl-${taskId}_runXTB"
    submitSubJob "$inpFilesList" "$outFilesList" "$wrkDir" "$runManyXTBScript" "${taskId}_runXTB"

    # Wait for completion of the job
    echo "Start waiting for completion of XTB master job..."
    for i in $(seq 1 "$maxWaitForHPC")
    do
        if [ "$i" == $maxWaitForHPC ]
        then
            sleep "$stepFirstWaitForHPCXTB$unitFirstWaitForHPC"
        else
            sleep "$stepWaitForHPCXTB$unitWaitForHPC"
        fi
        dateTime=$(date)
        if [ ! -f "$tcl" ]
        then
            echo "$dateTime Waiting for $tcl"
            if [ "$i" == $maxWaitForHPC ]
            then
                errMsg="#WaitingXTB: time limit reached (task abandoned)"
                abandon "$locInpSDF" "$E_OPTERROR"
            fi
        else
            echo "$dateTime All done: stop waiting."
            break;
        fi
    done

    # Check outcome, and collect results
    echo "Checking output files of XTB jobs..."
    for j in $(seq 0 $((${#toXTBStateLabels[@]}-1)) )
    do
        # Check output is found and not containing error
        out="${locSubJobXTBOutList[$j]}"
        if [ ! -f "$out" ]
        then
            errMsg="#RetrieveXTBResults: output file $out not found."
            abandon "$locInpSDF" "$E_OPTERROR"
        else
            # NB: here one could want to keep running if some, specific subjob
            # does not terminate properly.
            if grep -q "MOL_ERROR" "$out" > /dev/null
            then
                abandonDueToChild "$out"
            fi
        fi

        if [ "$cleanup" == 0 ]
        then
            rm -f "$wrkDir/tcl-${taskId}_runXTB"
            rm -f "$wrkDir/${taskId}_runXTB.nho"
        fi
    done

    mv "$wrkDir/${taskId}_runXTB.log" "$wrkDir/${molNum}_run_XTB.log"
}


#
# Function to submit a batch of parallel DFT jobs and wait for completion. No
# argument required: this function uses global variables as input. Namely:
# 'allSdfToDFT' and 'toDFTStateLabels'
#

function submitParallelBatchDFT(){
    # Submit all parallel sub-jobs
    echo "Submitting batch of DFT jobs: ${toDFTStateLabels[@]}"
    locSubJobDFTIDList=()
    locSubJobDFTOutList=()
    for idft in $(seq 0 $((${#toDFTStateLabels[@]}-1)) )
    do
        sjLab="${toDFTStateLabels[idft]}"
        sjId="$taskId.${sjLab}#"
        # WARNING the "_outSub" string is assumed in downstread scripts
        sjOut="$wrkDir/${molNum}_outSubDFT-$sjLab.sdf"
        locInpSDF="${allSdfToDFT[idft]}"
        submitSubJob "$locInpSDF" "$sjOut" "$wrkDir" "$runDFTScript" "$sjId"
        locSubJobDFTIDList+=("$sjId")
        locSubJobDFTOutList+=("$sjOut")
        allDFTSubJobsLabelList+=("$sjLab")
        allDFTSubJobsIDList+=("$sjId")
        allDFTSubJobsOutList+=("$sjOut")
    done

    # Wait for completion of parallel jobs
    echo "Start waiting for completion of DFT jobs..."
    allDFTDone=0
    for i in $(seq 1 "$maxWaitForHPC")
    do
        if [ "$i" == $maxWaitForHPC ]
        then
            sleep "$stepFirstWaitForHPC$unitFirstWaitForHPC"
        else
            sleep "$stepWaitForHPC$unitWaitForHPC"
        fi
        dateTime=$(date)
        running=0
        for j in $(seq 0 $((${#toDFTStateLabels[@]}-1)) )
        do
            tcl="$wrkDir/tcl-${locSubJobDFTIDList[$j]}"
            if [ ! -f "$tcl" ]
            then
                running=$((running+1))
                echo "$dateTime Still waiting for $tcl"
                if [ "$i" == $maxWaitForHPC ]
                then
                    errMsg="#WaitingDFT: time limit reached (task abandoned)"
                    abandon "$locInpSDF" "$E_OPTERROR"
                fi
            else
                if [ "$j" == $((${#toDFTStateLabels[@]}-1)) ] && [ "$running" == 0 ]
                then
                    allDFTDone=1
                    break;
                fi
            fi
        done
        if [ "$allDFTDone" == 1 ]
        then
            echo "$dateTime All done: stop waiting."
            break;
        fi
    done

    # Check outcome, and collect results
    echo "Checking output files of DFT jobs..."
    for j in $(seq 0 $((${#toDFTStateLabels[@]}-1)) )
    do
        # Check output is found and not containing error
        out="${locSubJobDFTOutList[$j]}"
        if [ ! -f "$out" ]
        then
            errMsg="#RetrieveDFTResults: output file $out not found."
            abandon "$locInpSDF" "$E_OPTERROR"
        else
            # NB: here one could want to keep running if some, specific subjob
            # does not terminate properly.
            if grep -q "MOL_ERROR" "$out" > /dev/null
            then
                abandonDueToChild "$out"
            fi
        fi

        mv "$wrkDir/${locSubJobDFTIDList[$j]}.log" "$wrkDir/${molNum}_DFT_${locSubJobDFTIDList[$j]}.log"
        if [ "$cleanup" == 0 ]
        then
            rm -f "$wrkDir/tcl-${locSubJobDFTIDList[$j]}"
            rm -f "$wrkDir/${locSubJobDFTIDList[$j]}.nho"
        fi
    done
}


#
# Function to submit a batch of DFT jobs to be executed sequentially by a single
# child job, and wait for completion.
# No argument required: this function uses global variables as input. Namely:
# 'allSdfToDFT' and 'toDFTStateLabels'
#
# Many variable names refer to DFT because the functionality was originally meant
# for DFT, but your should ignore this labeling.
#

function submitDFTAllInOnce(){
    # Submit all parallel sub-jobs
    echo "Submitting batch of DFT jobs: ${toDFTStateLabels[@]}"
    locSubJobDFTOutList=()
    inpFilesList=""
    outFilesList=""
    for idft in $(seq 0 $((${#toDFTStateLabels[@]}-1)) )
    do
        sjLab="${toDFTStateLabels[idft]}"
        # WARNING the "_outSub" string is assumed in downstread scripts
        sjOut="$wrkDir/${molNum}_outSubDFT-$sjLab.sdf"
        inpFilesList="$inpFilesList,${allSdfToDFT[idft]}"
        outFilesList="$outFilesList,$sjOut"
        allDFTSubJobsOutList+=("$sjOut")
        locSubJobDFTOutList+=("$sjOut")
    done

    tcl="$wrkDir/tcl-${taskId}_runDFT"
    submitSubJob "$inpFilesList" "$outFilesList" "$wrkDir" "$runManyDFTScript" "${taskId}_runDFT"

    # Wait for completion of the job
    echo "Start waiting for completion of DFT master job..."
    for i in $(seq 1 "$maxWaitForHPC")
    do
        if [ "$i" == $maxWaitForHPC ]
        then
            sleep "$stepFirstWaitForHPC$unitFirstWaitForHPC"
        else
            sleep "$stepWaitForHPC$unitWaitForHPC"
        fi
        dateTime=$(date)
        if [ ! -f "$tcl" ]
        then
            echo "$dateTime Waiting for $tcl"
            if [ "$i" == $maxWaitForHPC ]
            then
                errMsg="#WaitingDFT: time limit reached (task abandoned)"
                abandon "$locInpSDF" "$E_OPTERROR"
            fi
        else
            echo "$dateTime All done: stop waiting."
            break;
        fi
    done

    # Check outcome, and collect results
    echo "Checking output files of DFT jobs..."
    for j in $(seq 0 $((${#toDFTStateLabels[@]}-1)) )
    do
        # Check output is found and not containing error
        out="${locSubJobDFTOutList[$j]}"
        if [ ! -f "$out" ]
        then
            errMsg="#RetrieveDFTResults: output file $out not found."
            abandon "$locInpSDF" "$E_OPTERROR"
        else
            # NB: here one could want to keep running if some, specific subjob
            # does not terminate properly.
            if grep -q "MOL_ERROR" "$out" > /dev/null
            then
                abandonDueToChild "$out"
            fi
        fi

        if [ "$cleanup" == 0 ]
        then
            rm -f "$wrkDir/tcl-${taskId}_runDFT"
            rm -f "$wrkDir/${taskId}_runDFT.nho"
        fi
    done

    mv "$wrkDir/${taskId}_runDFT.log" "$wrkDir/${molNum}_run_DFT.log"
}


# Function to submit a batch of pre opt DFT jobs for pseudo transition states 
# to be executed sequentially by a single child job, and wait for completion.
# No argument required: this function uses global variables as input. Namely:
# 'allSdfToDFT' and 'toDFTStateLabels'
#
# Many variable names refer to DFT because the functionality was originally meant
# for DFT, but your should ignore this labeling.
#
function submitDFTPseudoTS(){
    # Submit all parallel sub-jobs
    echo "Submitting batch of DFT jobs: ${preDFTOptLabels[@]}"
    locSubJobDFTOutList=()
    inpFilesList=""
    outFilesList=""
    for idft in $(seq 0 $((${#preDFTOptLabels[@]}-1)) )
    do
        sjLab="${preDFTOptLabels[idft]}"
        # WARNING the "_outSub" string is assumed in downstread scripts
        sjOut="$wrkDir/${molNum}_outSubPreDFT-$sjLab.sdf"
        inpFilesList="$inpFilesList,${allSdfToPreDFT[idft]}"
        outFilesList="$outFilesList,$sjOut"
        allPreDFTSubJobsOutList+=("$sjOut")
        locSubJobPreDFTOutList+=("$sjOut")
    done

    tcl="$wrkDir/tcl-${taskId}_runDFT"
    submitSubJob "$inpFilesList" "$outFilesList" "$wrkDir" "$runPseudoTSDFTPreOpt" "${taskId}_runDFT"

    # Wait for completion of the job
    echo "Start waiting for completion of DFT master job..."
    for i in $(seq 1 "$maxWaitForHPC")
    do
        if [ "$i" == $maxWaitForHPC ]
        then
            sleep "$stepFirstWaitForHPC$unitFirstWaitForHPC"
        else
            sleep "$stepWaitForHPC$unitWaitForHPC"
        fi
        dateTime=$(date)
        if [ ! -f "$tcl" ]
        then
            echo "$dateTime Waiting for $tcl"
            if [ "$i" == $maxWaitForHPC ]
            then
                errMsg="#WaitingDFT: time limit reached (task abandoned)"
                abandon "$locInpSDF" "$E_OPTERROR"
            fi
        else
            echo "$dateTime All done: stop waiting."
            break;
        fi
    done

    # Check outcome, and collect results
    echo "Checking output files of DFT jobs..."
    for j in $(seq 0 $((${#preDFTOptLabels[@]}-1)) )
    do
        # Check output is found and not containing error
        out="${locSubJobPreDFTOutList[$j]}"
        if [ ! -f "$out" ]
        then
            errMsg="#RetrieveDFTResults: output file $out not found."
            abandon "$locInpSDF" "$E_OPTERROR"
        else
            # NB: here one could want to keep running if some, specific subjob
            # does not terminate properly.
            if grep -q "MOL_ERROR" "$out" > /dev/null
            then
                abandonDueToChild "$out"
            fi
        fi

        if [ "$cleanup" == 0 ]
        then
            rm -f "$wrkDir/tcl-${taskId}_runDFT"
            rm -f "$wrkDir/${taskId}_runDFT.nho"
        fi
    done

    mv "$wrkDir/${taskId}_runDFT.log" "$wrkDir/${molNum}_run_PreDFT.log"
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
# Function to remove a property from an SDF file.
# @param $1 the name of the property to be removed
# @param $2 the SDF file to alter
#
function removePropertyFromSDF() {
    propName="$1"
    file="$2"
    awk -v nlines=2 -v ptrn=$propName '$0~ptrn {for (i=0; i<nlines; i++) {getline}; next} 1' "$file" > "${file}_tmp"
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
    res=$?
    if [ "$res" != 0 ]
    then
        cp "$latestSDF" "$outSDF"
        addPropertyToSingleMolSDF "MOL_ERROR" "$errMsg" "$outSDF"
    fi
    echo " "
    echo "ERROR MESSAGE: "
    echo "$errMsg"
    echo " "
    #NB: the exit code is detected by DENOPTIM
    exit "$es" 
    #NB: we trap the EXIT signal 
}

#
# Abandon this script due to some error in child processes. The child process
# has already reported the error in the output SDF file.
# @param $1 the pathname of the latest result (i.e., an SDF file)
#
function abandonDueToChild {
    latestSDF="$1"
    errMsg=$(grep -A1 "MOL_ERROR" "$latestSDF" | tail -n 1)
    es=$(grep -A1 "EXIT_STATUS" "$latestSDF" | tail -n 1)
    cp "$latestSDF" "$outSDF"
    echo " "
    echo "ERROR MESSAGE from: $latestSDF (Status: $es)"
    echo "'""$errMsg""'"
    echo " "
    #NB: the exit code is detected by DENOPTIM
    exit "$es" 
    #NB: we trap the EXIT signal 
}

#
# Function to submit the modeling of a single state
# @param $1 pathname of the SDF with the chemical system used as input
# @param $2 pathname of the SDF file containing the output chemical system
# @param $3 pathname of work space
# @param $4 pathname of the BASH script to be called for this task
# @param $5 identification string for the child job
#
function submitSubJob() {
    in="$1"
    out="$2"
    wd="$3"
    script="$4"
    id="$5"
    cd "$wd"
    nohup "$script" "$in" "$out" "$wd" "$id" > "$id.nho" 2>&1&
    sleep 1
    errLines=0
    if [ -f "$id.nho" ]; then
        errLines=$(wc -l "$id".nho | awk '{print $1}')
    fi
    if [ "$errLines" != 0 ]; then
        errMsg="#SubmitNoHup: error message from nohup"
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    cd "$locDir"
}

#
# Perform all tasks to be done when exiting the script for whatever reason
#
function finish {
    finishTime=$(date +%s)
    runTime=$(($finishTime - $beginTime))
    date
    echo "Finished in $runTime seconds" 
    #Add here any mandatory task to be run on termination
}
trap finish EXIT


#
# Check N-Ru-Carbene angle to detect outliers
#
function detectNRuCOutlier() {
structureId=$1
inFile=$wrkDir/${molNum}_outSubXTB-$structureId.sdf
tolerance="55.0"
echo "Starting ACC..."
smallestAngle=$( autocompchem -p -t MeasureGeomDescriptors --verbosity 1 --infile "$inFile" --smarts "ANG [\$([#7](~[#44])=[#6]=[#8,#16])] [#44] [\$([#6;X3](~[#7;X3])~[#6,#7]),\$([#6;X3]1~[#6;X3]~[#6;X3]~[#7;X3]~[#6;X3]~[#6;X3]~1)]" --onlybonded true | grep ANG | grep -Eo "= [0-9]{1,3}\.[0-9]{1,16}" | awk '{print $2}' | sort -n | head -n 1 )
echo "ACC completed: smallest angle was: $smallestAngle degrees"
if (( $(echo "$smallestAngle < $tolerance" | bc -l) ))
    then 
    errMsg="#outlierDetection: N-Ru-Carbene angle is too small."
    abandon "$inpSDF" "$E_OPTERROR"
fi
}

#
# Check shortest non-bonded (and non 1,3 bonded) Ru,H distance to detect outliers.
#
function closeContact() {
structureId=$1
inFile=$wrkDir/${molNum}_outSubPreDFT-$structureId.sdf
echo -e "VERBOSITY: 1\nTASK: AnalyzeVDWClashes\nINFILE: "$inFile"\n\$TARGETSMARTS: [#1] [#44]\nCUTOFF: 0.1\nALLOWANCE: -0.2\nALLOWANCE13: 5.0" > $wrkDir/${molNum}_${structureId}_clash.param.tmp
echo "Checking $structureId for Ru,H close contacts:"
echo "Starting ACC..."
closestContact=$( autocompchem -p "$wrkDir/${molNum}_${structureId}_clash.param.tmp" | grep AtomClash | grep -E "H[0-9]{1,3}:Ru[0-9]{1,3}|Ru[0-9]{1,3}:H[0-9]{1,3}" | grep -oE "d: [0-9]\.[0-9]{1,16}" | awk '{print $2}' | sort -n | head -n 1 | while read n; do echo "scale=10; $n / 1" | bc -l; done )
echo "ACC completed: closest non bonded Ru,H contanct for $structureId: $closestContact Angstrom"
rm $wrkDir/${molNum}_${structureId}_clash.param.tmp
}

#
# Check shortest non-bonded (and non 1,3 bonded) Ru,H distance of imported candidate.
#
function importedCloseContact() {
inFile=$1
echo -e "VERBOSITY: 1\nTASK: AnalyzeVDWClashes\nINFILE: "$inFile"\n\$TARGETSMARTS: [#1] [#44]\nCUTOFF: 0.1\nALLOWANCE: -0.2\nALLOWANCE13: 5.0" > $wrkDir/${molNum}_clash.param.tmp
echo "Checking for Ru,H close contacts:"
echo "Starting ACC..."
importedClosestContact=$( autocompchem -p "$wrkDir/${molNum}_clash.param.tmp" | grep AtomClash | grep -E "H[0-9]{1,3}:Ru[0-9]{1,3}|Ru[0-9]{1,3}:H[0-9]{1,3}" | grep -oE "d: [0-9]\.[0-9]{1,16}" | awk '{print $2}' | sort -n | head -n 1 | while read n; do echo "scale=10; $n / 1" | bc -l; done )
echo "ACC completed: closest non bonded Ru,H contanct: $importedClosestContact Angstrom"
rm $wrkDir/${molNum}_clash.param.tmp
}

#
# Recaluclate fitness of an old FIT file. This function reads the energies from an old file,
# and recalculates the fitness according to the fitness definition of the function. $1 = Path to fitness file. 
#
function recalculateFitness() {
fitFile=$1
# Setting energy variables from original fitness file.
freeEnergyA="$( grep -A1 "<freeEnergyA>" $fitFile | tail -n 1 )"
freeEnergyF="$( grep -A1 "<freeEnergyF>" $fitFile | tail -n 1 )"
freeEnergyE="$( grep -A1 "<freeEnergyE>" $fitFile | tail -n 1 )"
freeEnergyC="$( grep -A1 "<freeEnergyC>" $fitFile | tail -n 1 )"
freeEnergyL="$( grep -A1 "<freeEnergyL>" $fitFile | tail -n 1 )"
freeEnergyX="$( grep -A1 "<freeEnergyX>" $fitFile | tail -n 1 )"
freeEnergyZ="$( grep -A1 "<freeEnergyZ>" $fitFile | tail -n 1 )"
freeEnergyD="$( grep -A1 "<freeEnergyD>" $fitFile | tail -n 1 )"
G_Propene="-117.738020895"
G_Ethene="-78.4716324751"
G_HoveydaProd="-502.204109499"
G_SIMes="-924.27249048"
G_HG_RuCl2_SIMes="-2402.13728523"
DG_referenceProductionBarrier=".0222529298"
hartree_to_kcalmol="627.49467516"
# Calculating new fitness components.
# Descriptor 1: (C-Z) - (C-X).
coef1="4"
desc1=$( echo "$coef1 * ( $hartree_to_kcalmol * ( ( $freeEnergyZ - $freeEnergyD ) - ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene ) ) )" | bc -l )
# Descriptor 2: A relative to C (Exponential, Low values)
gain2="0.7"
coef2="-1"
cutoff2="9"
predesc2=$( echo "( $coef2 * ( e( ( $hartree_to_kcalmol * ( $freeEnergyC + $G_HoveydaProd - $freeEnergyA - 2*${G_Ethene} ) + $cutoff2 ) * l( $gain2 ) ) ) )" | bc -l )
if [ "$( echo "$predesc2 <= -100" | bc -l )" == "1" ]
then
    desc2="-100"
else
    desc2="$predesc2"
fi
# Descriptor 3: A relative to C (Exponential, High values)
gain3="1.5"
coef3="-1"
cutoff3="-9"
predesc3=$( echo "( $coef3 * ( e( ( $hartree_to_kcalmol * ( $freeEnergyC + $G_HoveydaProd - $freeEnergyA - 2*${G_Ethene} ) + $cutoff3 ) * l( $gain3 ) ) ) )" | bc -l )
if [ "$( echo "$predesc3 <= -100" | bc -l )" == "1" ]
then
    desc3="-100"
else
    desc3="$predesc3"
fi
 # Descriptor 4: Exchange of L ligand with SIMes
gain4="2.0"
coef4="-3"
cutoff4="-14.0"
predesc4=$( echo "( $coef4 * ( e( ( $hartree_to_kcalmol * ( $freeEnergyE + $G_SIMes - $G_HG_RuCl2_SIMes - $freeEnergyL ) + $cutoff4 ) * l( $gain4 ) ) ) )" | bc -l )
if [ "$( echo "$predesc4 <= -100" | bc -l )" == "1" ]
then
    desc4="-100"
else
    desc4="$predesc4"
fi
# Descriptor 5: Production Barrier relative to HG Ru SIMes
coef5="-0.5"
desc5=$( echo "( $coef5 * ( $hartree_to_kcalmol * ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene - $DG_referenceProductionBarrier ) ) )" | bc -l )
# Descritor 6: disfavouring non trans precursors
gain6="0.5"
coef6="-1"
cutoff6="-3.5"
predesc6=$( echo "( $coef6 * ( e( ( $hartree_to_kcalmol * ( $freeEnergyF - $freeEnergyA ) + $cutoff6 ) * l( $gain6 ) ) ) )" | bc -l )
if [ "$( echo "$predesc6 <= -100" | bc -l )" == "1" ]
then
    desc6="-100"
else
    desc6="$predesc6"
fi
# Calculating overall fitness
fitness=$( echo " $desc1 + $desc2 + $desc3 + $desc4 + $desc5 + $desc6 " | bc -l )
#echo "done"
export desc1 desc2 desc3 desc4 desc5 desc6
echo "$fitness"
}

###############################################################################
# Main
###############################################################################

#
# Parse arguments
#
if [ "$#" -lt 5 ]
then
    echo " "
    echo "Usage: `basename $0` required number of arguments not supplied"       
    echo "5 parameters must be supplied (in this order):"
    echo " <input.sdf>  <output.sdf> <workingDirectory> <taskID> <UIDFile>"
    echo " "
    exit 1
fi
inpSDF="$1"
outSDF="$2"
wrkDir="$3"
taskId="$4"
UIDFILE="$5"
locDir="$(pwd)"
if [ "$#" -eq 6 ]
then 
    echo " "
    echo " Usage of previous data not implemented! "
    exit 1;
    #echo " Importing data to recover information from previous runs. "
    #readPrevDataSettings $6
fi

molName=`basename "$inpSDF" .sdf`
molNum=`basename "$molName" "$inputLabelExt"`
preOutSDF="$wrkDir/preOut_${molNum}.sdf"
MOLUID="$( grep -A1 UID "$wrkDir/$inpSDF" | tail -n 1)"

# Initiates experiment shut down. shutDownRun="1" -> Do it; shutDownRun="0" -> Don't. 
# Shut down will cause all new jobs to abandon, and UIDs of abandoned jobs will
# be stored in a text file.
shutDownRun="0"
if [ "$shutDownRun" -eq "1" ]
then
    echo "$MOLUID" >> "$wrkDir/shutDownAbandonUIDs.txt"
    errMsg="#ShutDown: abandoned due to experiment shut down."
    abandon "$inpSDF" "$E_OPTERROR"
fi

#
# Setup Log file
#
log="$wrkDir/${molNum}_FProvider.log"
# From here redirect stdout and stderr to log file
exec > $log
exec 2>&1
echo "Starting fitness calculation (ID:$taskId) at $beginTime" 


#
# Hold trap: humans can hold all machine-controlled fitness provider runs
#
for i in $(seq 1 $holdMaxWait)
do
    if [ -f "$flagHoldFitnessProvider" ]
    then
        echo "HOLDING PID:$$ => Fitness provider for $molName " >> "$listHeldPIDs"
        dateTime=$(date)
        echo "Fitness provider held by user - $dateTime"
    else
        break
    fi
    sleep "$holdStep$holdTimeUnit"
done


#
# If there are partial results available, skip conf. search and PM6 refinement,
# and recover data from previous run
#
allSubJobsLabelList=()
allSubJobsIDList=()
allSubJobsOutList=()
selectedValues=()
selectedLabels=()
if [ "$isRestartRun" == 0 ]
then
    #
    # Fake fitness can be used for debugging and testing purposes
    #
    if [ "$useFakeFitness" == 0 ]
    then
        echo "Fake fitness! Edit script $0 to do real computations"
        cp "$inpSDF" "$preOutSDF"
        fitness=$(date +%N)
        addPropertyToSingleMolSDF "FITNESS" "$fitness" "$preOutSDF"
        mv "$preOutSDF" "$outSDF"
        exit 0
    fi

    #
    # Run real task
    #
    echo "Recovering data from previous run..."
    oldDir="${pathCompletedStatesDir[0]}"
    if [ ! -d "$oldDir" ]
    then
        echo "ERROR! Cannot find directory with old results."
        exit 1
    fi
    cp -r "$oldDir" "$wrkDir/data_from_previous_run"
    if [ $? != 0 ]
    then
        echo "ERROR! while copying folder $oldDir"
        exit 1
    fi
    oldLog="${pathCompletedStatesDir[0]}/${molNum}_FProvider.log"
    if [ ! -f "$oldLog" ]
    then
        echo "ERROR! Cannot find old log file for restarting the run."
        exit 1
    fi
    for rootLab in "${sortedLabelRoots[@]}"
    do
        # 
        # Identify the refined model with lowest energy for this root label
        #
        echo "Recovering lowest GOSE for state $rootLab..."
        if ! grep -q "Lowest energy among $rootLab" "$oldLog" > /dev/null
        then
            echo "ERROR! Cannot find data for '$rootLab' in old log file."
            exit 1
        fi
        lowestJobLab=$(grep "Lowest energy among $rootLab" "$oldLog" | awk '{print $6}')
        subJobLog=$(find "$oldDir" -name "*$lowestJobLab.log")
        if [ $? != 0 ]
        then
            echo "ERROR! Cannot find *$lowestJobLab.log file"
            exit 1
        fi
        lowestGOSE="GOSE0"
        if grep -q "New lowest PM6 energy conformation found" "$subJobLog"
        then
            lowestGOSE=$(grep -B1 "New lowest PM6 energy conformation found" "$subJobLog" | tail -n 2 | head -n 1 | sed 's/)/ /g' | awk '{print $8}')
        fi

        #
        # Extract the geometry
        #
        echo "Extracting geometry from previous ${lowestJobLab}_${lowestGOSE}..."
        # Setting new params
        accSprGOSEOutParFile="$wrkDir/${molName}_Sprt${lowestGOSE}Out.par"
        accSprGOSEOutLog="$wrkDir/${molName}_Sprt${lowestGOSE}Out.log"
        sprtGOSEOut="$wrkDir/${molNum}_outSub-$lowestJobLab.sdf"
        sprtGOSEInp="$oldDir/${molNum}_${lowestJobLab}_${lowestGOSE}.spartan"
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
        # Update variables that will be used later
        #
        selectedLabels+=("$lowestJobLab")
        allSubJobsLabelList+=("$lowestJobLab")
        allSubJobsOutList+=("$sprtGOSEOut")
        gm3dfInp="$wrkDir/${molNum}_GM3DF.sdf"
        cp "$oldDir/${molNum}_GM3DF.sdf" "$gm3dfInp"
    done
else
    #
    # Fake fitness can be used for debugging and testing purposes
    #
    if [ "$useFakeFitness" == 0 ]
    then
        echo "Fake fitness! Edit script to do real computations"
        cp "$inpSDF" "$preOutSDF"
        fitness=$(date +%N)
        addPropertyToSingleMolSDF "FITNESS" "$fitness" "$preOutSDF"
        mv "$preOutSDF" "$outSDF"
        exit 0
    fi
     
    #Checking for the existance of matching #_out.sdf files (completed molecules).
    echo "Searching for UID: $MOLUID"
    if [ importOld == "1" ]
    then
        knownMol="$( find $pathToOld -name "*_out.sdf" -exec grep -l "$MOLUID" {} \; | xargs grep -l "FITNESS" | xargs grep -L "IMPORTED" | head -n 1 )"
    fi

    if [[ ! -z $MOLUID ]]
    then
        if [[ -n $knownMol ]] && [[ ! $knownMol == *"standard input"* ]]
        then
            echo "Molecule found in archive: $knownMol"
            knownMolNum=$( echo $knownMol | grep -Eo "M[0-9]{8}" )
            #zSdf=$( echo "/volume/jek003/rutsxtb/evolutionary_design/results/Raw_data/$( echo $knownMol | awk -F"jek003/" '{print $2}' | awk -F"_FIT" '{print $1}' )/${knownMolNum}_outSubXTB-Z.sdf" )
            #if [[ -f "$zSdf" ]]
            #then
            #importedCloseContact "$zSdf"
            newGCode=$( grep -A1 "<GCODE>" "$wrkDir/${molNum}_out.sdf" | tail -n 1 )
            newGraphMsg=$( grep -A1 "<GraphMsg>" "$wrkDir/${molNum}_out.sdf" | tail -n 1 )
            #newCdkTitle=$( grep -A1 "<cdk:Title>" "$wrkDir/${molNum}_out.sdf" | tail -n 1 )
            newFirstLine=$( head -n 1 "$wrkDir/${molNum}_out.sdf" )
            importMessage="Imported from: $knownMol"
            cp "$knownMol" "$wrkDir/${molNum}_out.sdf"
            oldFirstLine=$( head -n 1 "$wrkDir/${molNum}_out.sdf" )
            newFitness=$( recalculateFitness "$knownMol" )
            #Replacing old GCode, GraphMsg and CdkTitle on imported #_out.sdf with new ones from #_out.sdf
            #if grep "DESCRIPTOR_8" $knownMol >/dev/null 2>/dev/null
            #then
            #    removePropertyFromSDF "DESCRIPTOR_8" "$wrkDir/${molNum}_out.sdf"
            #fi
            removePropertyFromSDF "GCODE" "$wrkDir/${molNum}_out.sdf"
            #removePropertyFromSDF "cdk:Title" "$wrkDir/${molNum}_out.sdf"
            removePropertyFromSDF "GraphMsg" "$wrkDir/${molNum}_out.sdf"
            portablesed "s/$oldFirstLine/$newFirstLine/g" "$wrkDir/${molNum}_out.sdf"
            # Descriptor 8: Close contact (Ru,H)
            #if [[ $( echo "$closestContact <= 1.95" | bc -l ) -eq "1" ]]
            #then
            #    desc8="-100"
            #else
            #    coef8="-1"
            #    desc8=$( echo "$coef8 * e( l( 10 ) * ( 21 - 10 * $importedClosestContact )  )" | bc -l )
            #fi
            #addPropertyToSingleMolSDF "DESCRIPTOR_8" "$desc8" "$wrkDir/${molNum}_out.sdf"
            addPropertyToSingleMolSDF "IMPORTED" "$importMessage" "$wrkDir/${molNum}_out.sdf"
            addPropertyToSingleMolSDF "GCODE" "$newGCode" "$wrkDir/${molNum}_out.sdf"
            #addPropertyToSingleMolSDF "cdk:Title" "$newCdkTitle" "$wrkDir/${molNum}_out.sdf"
            addPropertyToSingleMolSDF "GraphMsg" "$newGraphMsg" "$wrkDir/${molNum}_out.sdf"
            #fitComponents=$( grep -oEA1 "DESCRIPTOR_[0-9]" "$wrkDir/${molNum}_out.sdf" | grep DESCRIPTOR | while read f; do grep -A1 "$f" "$wrkDir/${molNum}_out.sdf" | tail -n 1; done )
            #newFitness=$( echo "$fitComponents" | awk '{total += $NF} END { print total }' )
            removePropertyFromSDF "DESCRIPTOR_1" "$wrkDir/${molNum}_out.sdf"
            removePropertyFromSDF "DESCRIPTOR_2" "$wrkDir/${molNum}_out.sdf"
            removePropertyFromSDF "DESCRIPTOR_3" "$wrkDir/${molNum}_out.sdf"
            removePropertyFromSDF "DESCRIPTOR_4" "$wrkDir/${molNum}_out.sdf"
            removePropertyFromSDF "DESCRIPTOR_5" "$wrkDir/${molNum}_out.sdf"
            removePropertyFromSDF "DESCRIPTOR_6" "$wrkDir/${molNum}_out.sdf"
            removePropertyFromSDF "FITNESS" "$wrkDir/${molNum}_out.sdf"
            recalculateFitness "$knownMol" >/dev/null
            addPropertyToSingleMolSDF "DESCRIPTOR_1" "$desc1" "$wrkDir/${molNum}_out.sdf"
            addPropertyToSingleMolSDF "DESCRIPTOR_2" "$desc2" "$wrkDir/${molNum}_out.sdf"
            addPropertyToSingleMolSDF "DESCRIPTOR_3" "$desc3" "$wrkDir/${molNum}_out.sdf"
            addPropertyToSingleMolSDF "DESCRIPTOR_4" "$desc4" "$wrkDir/${molNum}_out.sdf"
            addPropertyToSingleMolSDF "DESCRIPTOR_5" "$desc5" "$wrkDir/${molNum}_out.sdf"
            addPropertyToSingleMolSDF "DESCRIPTOR_6" "$desc6" "$wrkDir/${molNum}_out.sdf"
            addPropertyToSingleMolSDF "FITNESS" "$newFitness" "$wrkDir/${molNum}_out.sdf"
            exit 0
            #else
            #errMsg="#knownMol: Molecul already exists in the database."
            #abandon "$inpSDF" "$E_OPTERROR"
            #fi
        fi
    fi
    #
    # Run real task starting with the First parallel batch: 
    # build selected states and choose ligand conformation
    #
    submitParallelBatchMMSE "1st" "$inpSDF" firstSubJobsLabelList[@] "$build3DConfScript"

    #
    # Extract conformation of the ligand
    #
    echo "Extracting L-ligand conformation..."
    absMinSDF="$wrkDir/${molNum}_outSub-$absMinLabel.sdf"
    gm3dfInp="$wrkDir/${molNum}_GM3DF.sdf"
    gm3dfPar="$wrkDir/${molNum}_GM3DF.par"
    gm3dfLog="$wrkDir/${molNum}_GM3DF.log"
    cp "$absMinSDF" "$gm3dfInp"
    echo "FRG-STRUCTURESFILE=$gm3dfInp" > "$gm3dfPar"
    echo "FRG-VERBOSITY=1" >> "$gm3dfPar"
    echo "FRG-CUTTINGRULESFILE=$coreLigandCutRules" >> "$gm3dfPar"
    echo "FRG-MAXFRAGSIZE=500" >> "$gm3dfPar"
    echo "FRG-REJECTELEMENT=Ru" >> "$gm3dfPar"
    # Submit GM3DFragmenter
    cd "$wrkDir"
    denoptim -r FRG "$gm3dfPar" > "$gm3dfLog" 2>&1 
    cd "$locDir"
    rm -f "$wrkDir/MolFrag-ratio_${molNum}_"*
    rm -f "$wrkDir/CPMap_${molNum}_"*
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
    # Prepare input for generation of other states
    #
    #NB: filename and extension are defined also in other scripts!
    libFragsAndLigand="$wrkDir/${molNum}_libFragsAndLigand.sdf"
    mv "$fragmenterDir/Fragments.sdf"  "$libFragsAndLigand"
    echo "Collecting fragments for building other states (1)"
    wc -l "$libFragsAndLigand"
    coreLigGraph="$wrkDir/${molNum}_core-ligGraph.sdf"
    cp "$coreLigGraphTmpl" "$coreLigGraph"
    portablesed "s/MOLNAME/$molNum/g" "$coreLigGraph"


    #
    # Extract conformation of the X-ligand
    #
    echo "Extracting X-ligand conformation..."
    gm3dfPar="$wrkDir/${molNum}_GM3DF-X.par"
    gm3dfLog="$wrkDir/${molNum}_GM3DF-X.log"
    echo "FRG-STRUCTURESFILE=$gm3dfInp" > "$gm3dfPar"
    echo "FRG-VERBOSITY=1" >> "$gm3dfPar"
    echo "FRG-CUTTINGRULESFILE=$coreXLigandCutRules" >> "$gm3dfPar"
    echo "FRG-MAXFRAGSIZE=500" >> "$gm3dfPar"
    echo "FRG-REJECTELEMENT=Ru" >> "$gm3dfPar"
    # Submit GM3DFragmenter
    cd "$wrkDir"
    denoptim -r FRG "$gm3dfPar" > "$gm3dfLog" 2>&1
    cd "$locDir"
    rm -f "$wrkDir/MolFrag-ratio_${molNum}_"*
    rm -f "$wrkDir/CPMap_${molNum}_"*
    # Check outcome
    fragmenterDir="$( grep "Output" "$gm3dfLog" | awk '{print $11}' )"
    if [ ! -f "$gm3dfLog" ]; then
        errMsg="#GM3DFragmenter: $gm3dfLog not found."
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    if ! grep -q " Completed Fragmenter" "$gm3dfLog" ; then
        errMsg="#GM3DFragmenter: non-zero exit status"
        abandon "$inpSDF" "$E_OPTERROR"
    fi

    # Keep only the first fragment
    xligFrags="$fragmenterDir/Fragments.sdf" 
    nEndFirst=$(grep -n -m 1 "\$\$\$\$" "$xligFrags" | awk -F":" '{print $1}')
    head -n "$nEndFirst" "$xligFrags"  >> "$libFragsAndLigand"
    echo "Collecting fragments for building other states (2)"
    wc -l "$libFragsAndLigand"
    # Reformatting the new fragments file
    cat "$libFragsAndLigand" | sed '/<ATTACHMENT_POINT>/,+2 d' | sed 's/<CLASS>/<ATTACHMENT_POINTS>/' > "${libFragsAndLigand}.tmp"
    mv "${libFragsAndLigand}.tmp" "$libFragsAndLigand"

    #
    # Second parallel batch:
    # build intermediates/transition states with the chosen ligand conformation
    #
    submitParallelBatchMMSE "2nd" "$coreLigGraph" secondSubJobsLabelList[@] "$build3DAndRefine"
fi    


#
# Prepare the preliminary output
# NB: preOutSDF will be overwritten once DFT results are available
#
cp "$gm3dfInp" "$preOutSDF"
addPropertyToSingleMolSDF "UID" "$MOLUID" "$preOutSDF"


#
# Get the proper conformer/stereoisomer for pseudo TS pre DFT opt
#
allSdfToPreDFT=()
for rootLab in "${preDFTOptLabels[@]}"
do
    #
    # Get the proper conformer/stereoisomer
    #
    found=0
    for subJobLab in "${selectedLabels[@]}"
    do
        rootS="${subJobLab:0:1}"
        if [ "$rootS" == "$rootLab" ]
        then
            found=1
            break;
        fi
    done
    if [ "$found" -ne 1 ]
    then
        errMsg="#RecoverSDFForXTB: unable to find SDF for $rootLab"
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    echo "Preparing pre DFT calculation for species $rootLab (from $subJobLab)"
    for j in $(seq 0 $((${#allSubJobsLabelList[@]}-1)) )
    do
        if [ "$subJobLab" == "${allSubJobsLabelList[$j]}" ]
        then
            break;
        fi
    done
    sdfToPreDFT="${allSubJobsOutList[$j]}"


    #
    # Store link to SDF file that will be submitted to XTB
    #
    allSdfToPreDFT+=("$sdfToPreDFT")
    echo "Geometry to pre DFT: $sdfToPreDFT"
done


# Submit job
submitDFTPseudoTS

#
# Get the proper conformer/stereoisomer
#
allSdfToXTB=()
for rootLab in "${toXTBStateLabels[@]}"
do
    #
    # Get the proper conformer/stereoisomer
    #
    found=0
    for subJobLab in "${selectedLabels[@]}"
    do
        rootS="${subJobLab:0:1}"
        if [ "$rootS" == "$rootLab" ]
        then
            found=1
            break;
        fi
    done
    if [ "$found" -ne 1 ]
    then
        errMsg="#RecoverSDFForXTB: unable to find SDF for $rootLab"
        abandon "$inpSDF" "$E_OPTERROR"
    fi
    echo "Preparing XTB calculation for species $rootLab (from $subJobLab)"
    for j in $(seq 0 $((${#allSubJobsLabelList[@]}-1)) )
    do
        if [ "$subJobLab" == "${allSubJobsLabelList[$j]}" ]
        then
            break;
        fi
    done
    sdfToXTB="${allSubJobsOutList[$j]}"

    #
    # Store link to SDF file that will be submitted to XTB
    #
    if  [ "$rootS" == "X" ] || [ "$rootS" == "Z" ]
    then
        continue
        #allSdfToXTB+=("$wrkDir/${molNum}_outSubPreDFT-$rootS.sdf")
    else
        allSdfToXTB+=("$sdfToXTB")
    fi
    echo "Geometry to XTB: $sdfToXTB"
done

#
# Submit xTB to HPC
#

# Old way submits one job for each label
#allXTBSubJobsLabelList=()
#allXTBSubJobsIDList=()
#allXTBSubJobsOutList=()
# WARNING! uses allSdfToXTB[@], toDFTStateLabels[@] and sortedLabelRoots[@] 
# and I really meand those whith the 'DFT' string in the name!
#submitParallelBatchXTB

# Submit one job that does all XTB tasks
allXTBSubJobsOutList=()
submitXTBAllInOnce

# Check xtb output for outliers
if grep -A1 "GraphENC" "$wrkDir/${molNum}_out.sdf" | tail -n 1 | awk '{print  $2}' | grep -E "[0-9]{1,6}_83_1_[0-9]{1,5}" >/dev/null 2>/dev/null
then
    echo "Running Outlier detection Z/X"
    detectNRuCOutlier "Z"
    detectNRuCOutlier "X"
elif grep -A1 "GraphENC" "$wrkDir/${molNum}_out.sdf" | tail -n 1 | awk '{print  $2}' | grep -E "[0-9]{1,6}_92_1_[0-9]{1,5}" >/dev/null 2>/dev/null
then
    echo "Running Outlier detection Z/X"
    detectNRuCOutlier "Z"
    detectNRuCOutlier "X"
else
    echo "No NCX found. Outlier detection will not be run."
fi

#
# Get the proper conformer/stereoisomer
#
allSdfToDFT=()
for sdfFromXTB in "${allXTBSubJobsOutList[@]}"
do
    allSdfToDFT+=("$sdfFromXTB")
    echo "Geometry to DFT: $sdfFromXTB"
done
# Adding X, Z and D(from C) from DFT optimization.
mv "$wrkDir/${molNum}_outSubPreDFT-C.sdf" "$wrkDir/${molNum}_outSubPreDFT-D.sdf"
for sdfFromDft in "${preDFTtoSPLabels[@]}"
do
    allSdfToDFT+=("$wrkDir/${molNum}_outSubPreDFT-$sdfFromDft.sdf")
    echo "Geometry to DFT: $wrkDir/${molNum}_outSubPreDFT-$sdfFromDft.sdf"
done

#
# Submit DFT to HPC
#

# Old submission of one job for each molecule
#allDFTSubJobsLabelList=()
#allDFTSubJobsIDList=()
#allDFTSubJobsOutList=()
# WARNING! uses allSdfToDFT[@], toDFTStateLabels[@] and sortedLabelRoots[@] as global arrays
#submitParallelBatchDFT

# New submission of one job running all DFT tasks
allDFTSubJobsOutList=()
# WARNING! uses allSdfToDFT[@], toDFTStateLabels[@] and sortedLabelRoots[@] as global arrays
submitDFTAllInOnce

# Measuring closest non bonded Ru,H distance (stored in var: 'closestContact')
closeContact "Z" 
closeContact "X" 

# Variables containing the final energies in Hartree
freeEnergyA=$( grep -A1 "DFT-ENERGY" "$wrkDir/${molNum}_outSubDFT-A.sdf" | tail -n 1 )
freeEnergyC=$( grep -A1 "DFT-ENERGY" "$wrkDir/${molNum}_outSubDFT-C.sdf" | tail -n 1 )
freeEnergyE=$( grep -A1 "DFT-ENERGY" "$wrkDir/${molNum}_outSubDFT-E.sdf" | tail -n 1 )
freeEnergyX=$( grep -A1 "DFT-ENERGY" "$wrkDir/${molNum}_outSubDFT-X.sdf" | tail -n 1 )
freeEnergyZ=$( grep -A1 "DFT-ENERGY" "$wrkDir/${molNum}_outSubDFT-Z.sdf" | tail -n 1 )
freeEnergyL=$( grep -A1 "DFT-ENERGY" "$wrkDir/${molNum}_outSubDFT-L.sdf" | tail -n 1 )
freeEnergyF=$( grep -A1 "DFT-ENERGY" "$wrkDir/${molNum}_outSubDFT-F.sdf" | tail -n 1 )
freeEnergyD=$( grep -A1 "DFT-ENERGY" "$wrkDir/${molNum}_outSubDFT-D.sdf" | tail -n 1 )

#Values in Hatree
G_Propene="-117.738020895"
G_Ethene="-78.4716324751" 
G_HoveydaProd="-502.204109499"

#Reference for synthesizability (Ligand exchange) in Hartree
G_SIMes="-924.27249048"
G_HG_RuCl2_SIMes="-2402.13728523"

#Reference for production barrier ( "D -> X" barrier relative to Ru SIMes) in Hartree
DG_referenceProductionBarrier=".0222529298"

hartree_to_kcalmol="627.49467516"

# Descriptor 1: (A-Z) - (A-X). 
coef1="4"
desc1=$( echo "$coef1 * ( $hartree_to_kcalmol * ( ( $freeEnergyZ - $freeEnergyD ) - ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene ) ) )" | bc -l )

# Descriptor 2: A relative to C (Exponential, Low values)
# Identity: e( a * l(b)) = b^a
gain2="0.7"
coef2="-1"
cutoff2="9"
predesc2=$( echo "( $coef2 * ( e( ( $hartree_to_kcalmol * ( $freeEnergyC + $G_HoveydaProd - $freeEnergyA - 2*${G_Ethene} ) + $cutoff2 ) * l( $gain2 ) ) ) )" | bc -l )
if [ "$( echo "$predesc2 <= -100" | bc -l )" == "1" ]
    then
    desc2="-100"
    else
    desc2="$predesc2"
fi

# Descriptor 3: A relative to C (Exponential, High values)
gain3="1.5"
coef3="-1"
cutoff3="-9"
predesc3=$( echo "( $coef3 * ( e( ( $hartree_to_kcalmol * ( $freeEnergyC + $G_HoveydaProd - $freeEnergyA - 2*${G_Ethene} ) + $cutoff3 ) * l( $gain3 ) ) ) )" | bc -l )
if [ "$( echo "$predesc3 <= -100" | bc -l )" == "1" ]
    then
    desc3="-100"
    else
    desc3="$predesc3"
fi

# Descriptor 4: Exchange of L ligand with SIMes  
gain4="2.0"
coef4="-3"
cutoff4="-14.0"
predesc4=$( echo "( $coef4 * ( e( ( $hartree_to_kcalmol * ( $freeEnergyE + $G_SIMes - $G_HG_RuCl2_SIMes - $freeEnergyL ) + $cutoff4 ) * l( $gain4 ) ) ) )" | bc -l )
if [ "$( echo "$predesc4 <= -100" | bc -l )" == "1" ]
    then
    desc4="-100"
    else
    desc4="$predesc4"
fi

# Descriptor 5: Production Barrier relative to HG Ru SIMes
coef5="-0.5"
desc5=$( echo "( $coef5 * ( $hartree_to_kcalmol * ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene - $DG_referenceProductionBarrier ) ) )" | bc -l )

# Descritor 6: disfavouring non trans precursors
gain6="0.5"
coef6="-1"
cutoff6="-3.5"
predesc6=$( echo "( $coef6 * ( e( ( $hartree_to_kcalmol * ( $freeEnergyF - $freeEnergyA ) + $cutoff6 ) * l( $gain6 ) ) ) )" | bc -l )
if [ "$( echo "$predesc6 <= -100" | bc -l )" == "1" ]
    then
    desc6="-100"
    else
    desc6="$predesc6"
fi

#FITNESS (Sum of fitness from descriptors 1-8)
fitness=$( echo " $desc1 + $desc2 + $desc3 + $desc4 + $desc5 + $desc6 " | bc -l )

#Preparing final sdf file
addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_1" "desc1=( echo '(  4 * ( hartree_to_kcalmol * ( ( freeEnergyZ - freeEnergyD ) - ( freeEnergyX + 2*G_Ethene - freeEnergyD - 2*G_Propene ) ) )' | bc -l )" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_2" "desc2=( echo '( -1 * ( e( ( hartree_to_kcalmol * ( freeEnergyC + G_HoveydaProd - freeEnergyA - ( 2 * G_Ethene ) ) + 9 ) * l( 0.7 ) ) ) )' | bc -l )" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_3" "desc3=( echo '( -1 * ( e( ( hartree_to_kcalmol * ( freeEnergyC + G_HoveydaProd - freeEnergyA - ( 2 * G_Ethene ) ) + -9 ) * l( 1.5 ) ) ) )' | bc -l )" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_4" "desc4=( echo '( -3 * ( e( ( hartree_to_kcalmol * ( freeEnergyE + G_SIMes - G_HG_RuCl2_SIMes - freeEnergyL ) + -3 ) * l( 1.7 ) ) ) )' | bc -l )" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_5" "desc5=( echo '( -0.5 * ( hartree_to_kcalmol * ( freeEnergyX + G_HoveydaProd - freeEnergyA - ( 2 * G_Propene ) - DG_referenceProductioniBarrier ) ) )' | bc -l )" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_6" "desc6=( echo '( -1 * ( e( ( hartree_to_kcalmol * ( freeEnergyF - freeEnergyA ) + -3.5 ) * l( 0.5 ) ) ) )' | bc -l )" "$preOutSDF"
addPropertyToSingleMolSDF "CALCULATION_OF_OVERALL_FITNESS" "fitness=( echo ' desc1 + desc2 + desc3 + desc4 + desc5 + desc6 ' | bc -l )" "$preOutSDF"
addPropertyToSingleMolSDF "freeEnergyA" "${freeEnergyA}" "$preOutSDF"
addPropertyToSingleMolSDF "freeEnergyF" "${freeEnergyF}" "$preOutSDF"
addPropertyToSingleMolSDF "freeEnergyE" "${freeEnergyE}" "$preOutSDF"
addPropertyToSingleMolSDF "freeEnergyC" "${freeEnergyC}" "$preOutSDF"
addPropertyToSingleMolSDF "freeEnergyL" "${freeEnergyL}" "$preOutSDF"
addPropertyToSingleMolSDF "freeEnergyX" "${freeEnergyX}" "$preOutSDF"
addPropertyToSingleMolSDF "freeEnergyZ" "${freeEnergyZ}" "$preOutSDF"
addPropertyToSingleMolSDF "freeEnergyD" "${freeEnergyD}" "$preOutSDF"
addPropertyToSingleMolSDF "G_Propene" "-117.738020895" "$preOutSDF"
addPropertyToSingleMolSDF "G_Ethene" "-78.4716324751" "$preOutSDF"
addPropertyToSingleMolSDF "G_HoveydaProd" "-502.204109499" "$preOutSDF"
addPropertyToSingleMolSDF "G_SIMes" "-924.27249048" "$preOutSDF"
addPropertyToSingleMolSDF "G_HG_RuCl2_SIMes" "-2402.13728523" "$preOutSDF"
addPropertyToSingleMolSDF "DG_referenceProductionBarrier" "0.01468968" "$preOutSDF"
addPropertyToSingleMolSDF "hartree_to_kcalmol" "627.49467516" "$preOutSDF"
addPropertyToSingleMolSDF "DFT_IDENTIFIER" "$( grep "Current timestamp (ms):" "$wrkDir/${molNum}_run_DFT.log" | awk '{print $4}' )" "$preOutSDF"
addPropertyToSingleMolSDF "DFT_MACHINE" "$( grep machine "$wrkDir/${molNum}_run_DFT.log" | awk '{print $2}' )" "$preOutSDF" 
addPropertyToSingleMolSDF "XTB_IDENTIFIER" "$( grep "Current timestamp (ms):" "$wrkDir/${molNum}_run_XTB.log" | awk '{print $4}' )" "$preOutSDF"
addPropertyToSingleMolSDF "XTB_MACHINE" "$( grep machine "$wrkDir/${molNum}_run_XTB.log" | awk '{print $2}' )" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_1" "$desc1" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_2" "$desc2" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_3" "$desc3" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_4" "$desc4" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_5" "$desc5" "$preOutSDF"
addPropertyToSingleMolSDF "DESCRIPTOR_6" "$desc6" "$preOutSDF"
addPropertyToSingleMolSDF "FITNESS" "$fitness" "$preOutSDF"

#
# Cleanup
#
if [ "$cleanup" == 0 ]
then
#    rm  "$inpSDF"
    rm -f "$gm3dfInp"
    rm -f "$gm3dfPar"
    rm -f "$gm3dfLog"
    rm -f "$libFragsAndLigand"
    rm -f "$graphEdTaskMol"
    rm -f "$GrEdPar"
    rm -f "$GrEdOut"
    rm -f "$GrEdLog"
    rm -f "$coreLigGraph"
    rm -f "$wrkDir/tcl-${molNum}"*
    rm -f "$wrkDir/${molNum}_libFragsAndLigand.sdf"
    rm -f "$wrkDir/Fragments_${molNum}_"*
    rm -f "$wrkDir/${molNum}_"*"Sprt"*
    rm -f "$wrkDir/${molNum}_MMSE2nd"*
    rm -f "$wrkDir/${molNum}_MMSE1st"*
    rm -f "$wrkDir/${molNum}_"*"_toXTB.sdf"
    rm -f "$wrkDir/${molNum}_"*"_toDFT.sdf"
    rm -f "$wrkDir/${molNum}_GM3DF.log"
    rm -f "$wrkDir/${molNum}_GM3DF.par"
    rm -f "$wrkDir/${molNum}_outSubDFT-"*".sdf"
    rm -f "$wrkDir/${molNum}_outSub-"*".sdf"
    rm -f "$wrkDir/${molNum}_run_DFT.log"
    rm -f "$wrkDir/${molNum}_run_XTB.log"
    rm -fr "$wrkDir/${molNum}_"*".spartan"
    rm -fr "$wrkDir/${molNum}_"*"CSSPR"*
    rm -fr "$wrkDir/xtb_${molNum}"
    
    #for subJobOutSDF in "${allSubJobsOutList[@]}"
    #do
    #    rm -f "$subJobOutSDF"
    #done
    mkdir "$wrkDir/${molNum}"
    mv "$wrkDir/${molNum}_"* "$wrkDir/${molNum}"
    #cp "$wrkDir/${molNum}/"*"outSub"* "$wrkDir/"
    echo "Making archive with log files and co." 
    tar -cpzf "$wrkDir/${molNum}.tar.gz" -C "$wrkDir" "${molNum}"
    rm -r "$wrkDir/${molNum}"
fi
mv "$preOutSDF" "$outSDF"


#
# Exit
# 
echo "Fitness evaluation completed for molecule $molNum!"
exit 0
