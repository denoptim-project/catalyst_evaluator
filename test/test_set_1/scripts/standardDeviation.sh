#!/bin/bash

candidates=("H0 HI HII_prime HII HC1Ph")
tags=("_1 _2 _3 _4 _5")
numCandidates="$( echo $candidates | wc -w )"
allStdDev=""
myDir="$( dirname "$0" )"

for candidate in $candidates
do
    buildString=""
    if [ -f "$candidate/${candidate}_out.sdf" ]; then
        if grep -q "<FITNESS>" "$candidate/${candidate}_out.sdf"; then
            buildString+="$candidate/${candidate}_out.sdf "
        fi
    fi
    for tag in $tags
    do
        if [ -f "$candidate$tag/$candidate${tag}_out.sdf" ]; then
            if grep -q "<FITNESS>" "$candidate$tag/$candidate${tag}_out.sdf"; then
                buildString+="$candidate$tag/$candidate${tag}_out.sdf "
            fi
        fi
    done
    echo "Candidate group: $candidate"
    stdDev="$( $myDir/standardDeviation.py $buildString )"
    stdDevD1="$( $myDir/standardDeviationD1.py $buildString )"
    stdDevD2="$( $myDir/standardDeviationD2.py $buildString )"
    stdDevD3="$( $myDir/standardDeviationD3.py $buildString )"
    stdDevW4="$( $myDir/standardDeviationW1.py $buildString )"
    stdDevW5="$( $myDir/standardDeviationW2.py $buildString )"
    stdDevW6="$( $myDir/standardDeviationW3.py $buildString )"
    stdDevW7="$( $myDir/standardDeviationW4.py $buildString )"
    mean="$( $myDir/mean.py $buildString )"
    meanD1="$( $myDir/mean_D1.py $buildString )"
    meanD2="$( $myDir/mean_D2.py $buildString )"
    meanD3="$( $myDir/mean_D3.py $buildString )"
    meanW4="$( $myDir/mean_W1.py $buildString )"
    meanW5="$( $myDir/mean_W2.py $buildString )"
    meanW6="$( $myDir/mean_W3.py $buildString )"
    meanW7="$( $myDir/mean_W4.py $buildString )"
    echo "Candidates in group:  $buildString"
    echo "mean FIT: $mean"
    echo "mean D1: $meanD1"
    echo "mean D2: $meanD2"
    echo "mean D3: $meanD3"
    echo "mean W4: $meanW4"
    echo "mean W5: $meanW5"
    echo "mean W6: $meanW6"
    echo "mean W7: $meanW7"
    echo "StdDev FIT: $stdDev"
    echo "stdDev D1: $stdDevD1"
    echo "stdDev D2: $stdDevD2"
    echo "stdDev D3: $stdDevD3"
    echo "stdDev W4: $stdDevW4"
    echo "stdDev W5: $stdDevW5"
    echo "stdDev W6: $stdDevW6"
    echo "stdDev W7: $stdDevW7"
    echo ""
    allStdDev+="$stdDev + "
done
averageStdDev="$( echo "scale=3; ( $allStdDev 0 ) / $numCandidates" | bc -l )"
echo "Average over all candidates: $averageStdDev"
