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
../evaluate_candidate.sh -i Mol1.sdf

# Assert consistency of results
echo "Running assertions..."
if [ ! -f Mol1/Mol1.tar.gz ]; then echo ERROR: missing archive ; exit -1 ; fi
if [ ! -f Mol1/Mol1_out.sdf ]; then echo ERROR: missing output file ; exit -1 ; fi
if ! grep -q '<FITNESS>' Mol1/Mol1_out.sdf ; then echo ERROR: fitness value not found in output; echo -1 ; fi
actualFitness=$(grep -A 1 '<FITNESS>' Mol1/Mol1_out.sdf | tail -n 1 | cut -c 1-5)
expectedFitness="34.4" # NB: we compare only these 4 digits
if [ "$expectedFitness" != "$actualFitness" ]; then echo "ERROR: result is '$actualFitness' instead of '$expectedFitness'" ; exit -1; fi

echo "Test PASSED!"
exit 0
