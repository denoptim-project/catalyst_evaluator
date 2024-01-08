#!/bin/bash
###############################################################################
#
# Script to start the evaluation of a candidate catalyst.
# For instructions run this script without any argument (e.g., ./scriptname.sh)
#
############################################################################### 

#
# Non-conda evnironment
#

# Pathname to the folder containing Tinker's executables
export TINKERBIN="/usr/local/TinkerLFMM_RCP_2/bin"

# Pathname to the Spartan executable
export SPARTANEXE="/usr/bin/spartan20"

# Set the following to 1 to make every waiting step check for completion with high frequency. This is meant for running tests that may fail quickly. Do not use it for production runs.
export HIGHFREQUENCY=0

# Set the following to 0 to allow running XTB jobs locally, if CPU loading allows so, or set it to 1 to force the sending of XTB jobs to a remote worker regardless of local CPU loading.
export RUNXTBLOCALLY=0

# Set this to the maximum number of main processes to run. Typically this is set to the number of cores of the machine hosting the evaluator.
export MAXPARALLELPROCESSES=32

# Version of 'sed' (specify either GNU or BSD)
export SEDSYNTAX="GNU"

################################################################################
#
# Settings you most likely do not need to change
#
################################################################################

export SHELL="/bin/bash"
export FFPARAMFILE="../data/params.MMFF94_CS-2.2"
myDir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
export REMOTEWORKERBRIDGE="$myDir/tools/RemoteWorkersBridge/submit_tool/submit.py"
export REMOTEWORKERBRIDGEHOME="$myDir/tools/RemoteWorkersBridge"

configurationFile="$myDir/tools/RemoteWorkersBridge/configuration"
if [ ! -f "$configurationFile" ]; then echo "ERROR: missing $configurationFile"; exit -1 ; fi
# WARNING: assumption that the 'wdirOnRemote' is the same on all the workers
# that need basis sets.
basisSetDir="$(grep wdirOnRemote "$configurationFile" | head -n 1 | awk -F'=' '{print $2}')/basisset"


################################################################################
#
# Utilities and functions
#
################################################################################

#
# Prints the usage instrution of this script
#
function printUsage() {
    cat <<EOF

 Usage:
 =====

 ./scriptname.sh -i <input> [-d <workspace>]

 where 

 -i <input> is a pathname to an SDF file containing the definition of
            a catalyst in terms of Denoptim graph
            (see https://github.com/denoptim-project/DENOPTIM).
 
 -d <workspace> is the pathname to a file system location that will be used
                as work space, i.e., we'll start processes from that location.
                Default value is '.'.

EOF
}

#
# Function that checks if the arguments of a command line option is present.
#
# @param $1 the argument index (0-based integer) of the argument right before
# the one we are checking the existence of.
# @param $2 a string identifying what $1 is supposed to be (used only for logging).
# @param $3 total number of arguments given to the script that calles this function.
# @param $4-$# original arguments of the script that calls this function.
# @return 1 when no argument is found for a command line option.
#

function checkArg() {
   local ii="$1"
   local option="$2"
   local tot="$3"
   local args=("${@:4}")
   ii=$((ii+1))
   if [ "$ii" -ge "$tot" ]
   then
      echo "ERROR! Missing argument for option '$option'."
      return 1
   fi
   local argument=${args[$i+1]}
   if [[ "$argument" == "-"* ]]
   then
      echo "ERROR! Missing argument for option '$option'."
      return 1
   fi
   return 0
}

#
# Function that checks whether the arguments of a command line option it valid.
# If not, it exits with an error.
#

function ensureGoodArg() {
   local option="$2"
   if ! checkArg $@
   then
     printUsage
     exit 1
   fi
}

if [ "$#" -lt 1 ]
then
    echo "ERROR! Missing input file option."
    printUsage
    exit -1
fi

inpPathName=""
baseWrkDir="$(pwd)"
args=("$@")
for ((i=0; i<$#; i++))
do
  arg="${args[$i]}"
  case "$arg" in
    "-i") ensureGoodArg "$i" "$arg" "$#";
          inpPathName=${args[$i+1]};;
    "-d") ensureGoodArg "$i" "$arg" "$#";
          baseWrkDir=${args[$i+1]};;
    *);;
  esac
done

if [ ! -f "$inpPathName" ]
then
    echo "ERROR! Input file '$inpPathName' not found!"
    printUsage
    exit -1
fi

inpFile="$(basename $inpPathName)"
runName=${inpFile%.*}
export WORKDIR="$baseWrkDir/$runName"
logFile="$runName.log"

################################################################################
#
# Main
#
################################################################################

if [ -d "$WORKDIR" ]; then
    echo " "
    echo "ERROR! $WORKDIR exists already!"
    echo " "
    exit
fi
mkdir "$WORKDIR"
if [ "$?" -ne 0 ]; then
    echo " "
    echo "ERROR! Cannot make work directory '$WORKDIR'!"
    echo " "
    exit
fi
chmod o-rwx "$WORKDIR"
chmod g-w "$WORKDIR"

cp -r "$myDir"/data/* "$WORKDIR"
cp "$myDir"/src/*sh "$WORKDIR"
cp "$inpPathName" "$WORKDIR"

# WARNING! We place the FF file in the $HOME
cp "$FFPARAMFILE" "$HOME/PARAMS.MMFF.DENOPTIM"
if [ ! -f "$HOME/PARAMS.MMFF.DENOPTIM" ]
then
    echo " "
    echo "ERROR! Could not copy modified force field file to \$HOME"
    echo " "
    exit
fi

cd "$WORKDIR"

filesToModify=$(find . -type f | xargs grep -l "UNSET_BASISSETDIR")
for f in $filesToModify
do
    echo "Updating file $f"
    if [ "$SEDSYNTAX" == "GNU" ]
    then
        sed -i "s|UNSET_BASISSETDIR|$basisSetDir|g" $f
    elif [ "$SEDSYNTAX" == "BSD" ]
    then
        sed -i '' "s|UNSET_BASISSETDIR|$basisSetDir|g" $f
    fi
done

echo "Now calling the fitness provider..."
./fitness_provider.sh "$inpFile" "$runName"_out.sdf "$WORKDIR" 1 MOLUID.txt
exitValue="$?"

if [ 0 -ne "$exitValue" ]; then
  echo Returning non-zero exit status: $exitValue
  exit $exitValue
fi
echo "All done!"
exit 0
