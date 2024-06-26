#!/bin/bash
#
# Script that runs the evaluation of a small test catalyst
#

# Setup environment
echo "Loading the environment..."
eval "$(conda shell.bash hook)"
conda activate RuCatEvaluator
if [ 0 -ne "$?" ]; then
  echo "ERROR! Could not use Conda to activate the RuCatEvaluator environment!"
  exit -1
fi

# Run evaluation
testDir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo "Starting catalyst evaluation in $testDir..."
cd "$testDir"
../evaluate_catalyst.sh -i Mol1.sdf $@
if [ 0 -ne "$?" ]; then
  echo "ERROR! Non-zero exit status from catalyst evaluation"
  exit -1
fi

# Assert consistency of results
echo "Running assertions..."
if [ ! -f Mol1/Mol1.tar.gz ]; then echo ERROR: missing archive ; exit -1 ; fi
if [ ! -f Mol1/Mol1_out.sdf ]; then echo ERROR: missing output file ; exit -1 ; fi
if ! grep -q '<FITNESS>' Mol1/Mol1_out.sdf ; then echo ERROR: fitness value not found in output; echo -1 ; fi
actualFitness=$(grep -A 1 '<FITNESS>' Mol1/Mol1_out.sdf | tail -n 1 )
actualFitnessRounded=$( echo "$actualFitness" | xargs printf "%.*f\n" "1" )
expectedFitness="6.4" # NB: we compare only these 3
if [ "$expectedFitness" != "$actualFitnessRounded" ]; then echo "ERROR: result is '$actualFitnessRounded' instead of '$expectedFitness'" ; exit -1; fi

echo "Test PASSED!"
exit 0
