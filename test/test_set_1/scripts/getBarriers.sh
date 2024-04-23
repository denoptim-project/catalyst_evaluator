#!/bin/bash

candidates=("H0 HI HII_prime HII HC1Ph")
tags=("_1 _2 _3 _4 _5")
G_Propene="-117.738020895"
G_Ethene="-78.4716324751"
G_HoveydaProd="-502.204109499"
hartree_to_kcalmol="627.49467516"
myDir="$( dirname "$0" )"

# Main

buildOut=""
buildOut="$( echo "Candidate Barrier_metathesis Barrier_B-HE DDE_Met/BE" )"

for candidate in $candidates
do 
    counter="1"
    allXBarriers=""
    allZBarriers=""
    allDDE=""
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
                counter=$( echo "$counter + 1" | bc )
            fi
        fi
    done
    for candidate in $buildString
    do
        energyD="$( grep -A1 "<freeEnergyD>" "$candidate" | tail -n 1 )"
        energyX="$( grep -A1 "<freeEnergyX>" "$candidate" | tail -n 1 )"
        energyZ="$( grep -A1 "<freeEnergyZ>" "$candidate" | tail -n 1 )"
        barrierX="$( echo "scale=1; ( $hartree_to_kcalmol *( $energyX + 2*$G_Ethene - $energyD - 2*$G_Propene ) )" | bc -l )"
        barrierZ="$( echo "scale=1; ( $hartree_to_kcalmol *( $energyZ - $energyD ) )" | bc -l )"
        DDE="$( echo "scale=1; $barrierZ - $barrierX" | bc -l )"
        buildOut+="$( echo -e "\n$candidate $barrierX $barrierZ $DDE" )"
        allXBarriersSTD+="$barrierX "
        allZBarriersSTD+="$barrierZ "
        allXBarriers+="$barrierX + "
        allZBarriers+="$barrierZ + "
        allDDE+="$DDE + "
    done
    stdX="$( $myDir/standardDeviation_numbers.py $allXBarriersSTD)"
    stdZ="$( $myDir/standardDeviation_numbers.py $allZBarriersSTD)"
    allXBarriers+="0"
    allZBarriers+="0"
    allDDE+="0"
    meanX="$( echo "scale=3; ( $allXBarriers ) / $counter" | bc -l )"
    meanZ="$( echo "scale=3; ( $allZBarriers ) / $counter" | bc -l )"
    meanDDE="$( echo "scale=3; ( $allDDE ) / $counter" | bc -l )"
    buildOut+="$( echo -e "\nMEAN: $meanX $meanZ $meanDDE" )"
    buildOut+="$( echo -e "\nSTDV: $stdX $stdZ NaN" )"
    buildOut+="$( echo -e "\n------------------------------- ------------- ------------ -------------" )"
    allXBarriersSTD=""
    allZBarriersSTD=""
done
echo "$buildOut" | column -t
