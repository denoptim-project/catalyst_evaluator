#!/bin/bash

#
#   This script updates the fitness of desired files according to the fitness definition defined bellow. 
#   Make alterations to the current fitness as is pleased, and run the script see the effect on test set 1.
#  
#   Usage: UpdateFitness.sh paths_to_output_files_to_be_altered.txt
#
#

fitFileList="$1"

if ! [ -f "$fitFileList" ]
then
    echo "ERROR: Text file with fit file paths was not found."
    echo "Usage: updateFitness.sh <arg1>"
    echo "<arg1>: Text file containing paths to fitness files (sdf)"
    exit 0
fi

# Set to "1" to activate HQF instead of AF (This changes the constants to match the respective methods) 
HQF="0"  #  <---{ 0 = AF, 1 = HQF }

# FUNCTIOINS ########################################################################################################

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

# MAIN ##############################################################################################################
# Looping over original fitness files listed in the input.
cat $fitFileList | while read f
do
    # Setting energy variables from original fitness file.
    freeEnergyA="$( grep -A1 "<freeEnergyA>" $f | tail -n 1 )"
    freeEnergyF="$( grep -A1 "<freeEnergyF>" $f | tail -n 1 )"
    freeEnergyE="$( grep -A1 "<freeEnergyE>" $f | tail -n 1 )"
    freeEnergyC="$( grep -A1 "<freeEnergyC>" $f | tail -n 1 )"
    freeEnergyL="$( grep -A1 "<freeEnergyL>" $f | tail -n 1 )"
    freeEnergyX="$( grep -A1 "<freeEnergyX>" $f | tail -n 1 )"
    freeEnergyZ="$( grep -A1 "<freeEnergyZ>" $f | tail -n 1 )"
    freeEnergyD="$( grep -A1 "<freeEnergyD>" $f | tail -n 1 )"

    ########### AF/HQF Shared values  ######################
    
    hartree_to_kcalmol="627.49467516" # (kcal/mol*hartree)
    boltzmann_constant="0.0000000000000000000000138066" # (J/K)
    avogadro_number="602214000000000000000000" # (1/mol)
    jmol_to_kcalmol="0.00023901" # (kcal/J)
    temp="313.15" # in K ( 40 degrees C to match laboratory experiments)
    kT=$( echo "$boltzmann_constant * $avogadro_number * $jmol_to_kcalmol * $temp" | bc -l ) # in kcal/mol
    magnitude3="4.29872427" # magnitude<n> represents n=1,2,3 orders of magnitude increase of rate constant at 40 C according to Eyring in kcal/mol
    magnitude2="2.86581618"
    magnitude1="1.43290809"
    magnitude4="5.73163236"

    ########### HQF values below ############################

    if [[ $HQF == "1" ]]
    then
        DDG_reference_HGII="0.019174109"
        G_Propene="-117.734910132"
        G_Ethene="-78.477866552"
        G_HoveydaProd="-502.210801644"
        G_SIMes="-924.245711409"
        G_HG_RuCl2_SIMes="-2402.142783425"
        DG_referenceProductionBarrier="0.024552309"
        DG_referencePrecursorStabilityHII="0.0100161"
        DG_referencePrecursorStabilityHI="0.02165023"
        DDG_referencePrecursorStability="$( echo "$DG_referencePrecursorStabilityHI - $DG_referencePrecursorStabilityHII" | bc -l )"
        DG_referenceSynt="0.0162160110"

    ############ AF values below ##########################

    else
        DDG_reference_HGII="0.0222312922"
        G_Propene="-117.738020895"
        G_Ethene="-78.4716324751"
        G_HoveydaProd="-502.204109499"
        G_SIMes="-924.27249048"
        G_HG_RuCl2_SIMes="-2402.13728523"
        DG_referenceProductionBarrier="0.0222714248"
        DG_referencePrecursorStabilityHII="-0.0052201621"
        DG_referencePrecursorStabilityHI="0.0117775312"
        DDG_referencePrecursorStability="$( echo "$DG_referencePrecursorStabilityHI - $DG_referencePrecursorStabilityHII" | bc -l )"
        DG_referenceSynt="0.0226072217"
    fi

    ##########################################################
    # Descriptor 1: Boltzmann factor (C-Z) - (C-X) - DDG_ref(HG2).
    coef1="1"
    desc1=$( echo "$coef1 * ( e( $hartree_to_kcalmol * ( ( $freeEnergyZ - $freeEnergyD ) - ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene ) - $DDG_reference_HGII ) / $kT ) ) " | bc -l )
    DESCRIPTOR_DEFINITION_1='desc1=( echo "1 * ( e( $hartree_to_kcalmol * ( ( $freeEnergyZ - $freeEnergyD ) - ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene ) - $DDG_reference_HGII ) / $kT ) ) " | bc -l )'
 
    # Descriptor 2: Linear: (C-Z) - (C-X).
    coef2="1"
    desc2=$( echo "$coef2 * ( $hartree_to_kcalmol * ( ( $freeEnergyZ - $freeEnergyD ) - ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene ) ) )" | bc -l )
    DESCRIPTOR_DEFINITION_2='desc2=$( echo "1 * ( $hartree_to_kcalmol * ( ( $freeEnergyZ - $freeEnergyD ) - ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene ) ) )" | bc -l )'

    # Descriptor 3: Production Barrier relative to HG Ru SIMes
    coef3="1"
    threshold3=$( echo "( $hartree_to_kcalmol * $DG_referenceProductionBarrier )" | bc -l )
    DG_cat=$( echo "$hartree_to_kcalmol * ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene )" | bc -l )
    if [ "$( echo "$DG_cat >= $threshold3" | bc -l )" == "1" ]
    then
        desc3="0.0000000001"
    else
        desc3=$( echo "$coef3 * ( $threshold3 - $DG_cat )" | bc -l )
    fi
    DESCRIPTOR_DEFINITION_3='desc3=$( echo "1 * $hartree_to_kcalmol * ( $DG_referenceProductionBarrier - ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene ) )" | bc -l )'

    # Dynamic weight 1 (sigmoid): Catalyst activity.
    threshold4="$( echo "( $DG_referenceProductionBarrier * $hartree_to_kcalmol ) + $magnitude4" | bc -l )"
    w1=$( echo "1 / ( 1 + e( ( ( $hartree_to_kcalmol * ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene ) ) - $threshold4 ) / $kT ) )" | bc -l )
    WEIGHT_DEFINITION_1='w1=$( echo "1 / ( 1 + e( $hartree_to_kcalmol * ( ( ( $freeEnergyX + 2*$G_Ethene - $freeEnergyD - 2*$G_Propene ) ) - ( $DG_referenceProductionBarrier + $magnitude4 ) ) / $kT ) )" | bc -l )'
 
    # Dynamic weight 2 (sigmoid): Stability of precursor wrt. MCB.
    threshold5=$( echo "( $hartree_to_kcalmol*$DDG_referencePrecursorStability ) + $magnitude1" | bc -l )
    DDG_stability=$( echo "$hartree_to_kcalmol * ( $freeEnergyC + $G_HoveydaProd - $freeEnergyA - 2*${G_Ethene} - $DG_referencePrecursorStabilityHII )" | bc -l )
    if [ "$( echo "$DDG_stability >= 0" | bc -l )" == "1" ]
    then 
        abs_DDG_Stability="$DDG_stability"
    else
        abs_DDG_Stability="$( echo "-1 * $DDG_stability" | bc -l )"
    fi
    w2=$( echo "1 / ( 1 + e( ( $abs_DDG_Stability - $threshold5 ) / $kT ) )" | bc -l )
    WEIGHT_DEFINITION_2='w2=$( echo "1 / ( 1 + e( $hartree_to_kcalmol * ( ABS( $freeEnergyC + $G_HoveydaProd - $freeEnergyA - 2*${G_Ethene} - $DG_referencePrecursorStabilityHII ) - ( $DDG_referencePrecursorStability + $magnitude1 ) ) / $kT ) )" | bc -l )' 
 
    # Dynamic weight 3 (sigmoid): Ligand exchange ( H2IMes  --exchange--> L ).
    coef6="1"
    threshold6=$( echo "( $hartree_to_kcalmol * $DG_referenceSynt ) + $magnitude1" | bc -l )
    DG_synt=$( echo "$hartree_to_kcalmol * ( $freeEnergyE + $G_SIMes - $G_HG_RuCl2_SIMes - $freeEnergyL )" | bc -l )
    w3=$( echo "1 / ( 1 + e( ( $DG_synt - $threshold6 ) / $kT ) )" | bc -l )
    WEIGHT_DEFINITION_3='w3=$( echo "1 / ( 1 + e( $hartree_to_kcalmol * ( ( $freeEnergyE + $G_SIMes - $G_HG_RuCl2_SIMes - $freeEnergyL ) - ( $DG_referenceSynt  + $magnitude1 ) ) / $kT ) )" | bc -l )'
 
    # Dynamic weight 4 (sigmoid): disfavouring non trans precursors.
    threshold7="$magnitude2"
    DG_stereo=$( echo "$hartree_to_kcalmol * ( $freeEnergyF - $freeEnergyA )" | bc -l )    
    w4=$( echo "1 / ( 1 + e( ( - $DG_stereo + $threshold7 ) / $kT ) )" | bc -l )
    WEIGHT_DEFINITION_4='w4=$( echo "1 / ( 1 + e( $hartree_to_kcalmol * ( - ( $freeEnergyF - $freeEnergyA ) + $magnitude2 ) / $kT ) )" | bc -l )'

    #echo "Calculating overall fitness"
    fitness=$( echo "( $desc1 + $desc2 + $desc3 ) * $w1 * $w2 * $w3 * $w4" | bc -l )
    FITNESS_DEFINITION='fitness=$( echo "( $desc1 + $desc2 + $desc3 ) * $w1 * $w2 * $w3 * $w4" | bc -l )'
 
    removePropertyFromSDF "DESCRIPTOR_1" "$f"
    removePropertyFromSDF "DESCRIPTOR_2" "$f"
    removePropertyFromSDF "DESCRIPTOR_3" "$f"
    removePropertyFromSDF "DESCRIPTOR_4" "$f"
    removePropertyFromSDF "DESCRIPTOR_5" "$f"
    removePropertyFromSDF "DESCRIPTOR_6" "$f"
    removePropertyFromSDF "WEIGHT_1" "$f"
    removePropertyFromSDF "WEIGHT_2" "$f"
    removePropertyFromSDF "WEIGHT_3" "$f"
    removePropertyFromSDF "WEIGHT_4" "$f"
    removePropertyFromSDF "FITNESS" "$f"
    removePropertyFromSDF "DESCRIPTOR_DEFINITION_1" "$f"
    removePropertyFromSDF "DESCRIPTOR_DEFINITION_2" "$f"
    removePropertyFromSDF "DESCRIPTOR_DEFINITION_3" "$f"
    removePropertyFromSDF "DESCRIPTOR_DEFINITION_4" "$f"
    removePropertyFromSDF "DESCRIPTOR_DEFINITION_5" "$f"
    removePropertyFromSDF "DESCRIPTOR_DEFINITION_6" "$f"
    removePropertyFromSDF "WEIGHT_DEFINITION_1" "$f"
    removePropertyFromSDF "WEIGHT_DEFINITION_2" "$f"
    removePropertyFromSDF "WEIGHT_DEFINITION_3" "$f"
    removePropertyFromSDF "WEIGHT_DEFINITION_4" "$f"
    removePropertyFromSDF "FITNESS_DEFINITION" "$f"
    addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_1" "$DESCRIPTOR_DEFINITION_1" "$f" 
    addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_2" "$DESCRIPTOR_DEFINITION_2" "$f" 
    addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_3" "$DESCRIPTOR_DEFINITION_3" "$f" 
    addPropertyToSingleMolSDF "WEIGHT_DEFINITION_1" "$WEIGHT_DEFINITION_1" "$f" 
    addPropertyToSingleMolSDF "WEIGHT_DEFINITION_2" "$WEIGHT_DEFINITION_2" "$f" 
    addPropertyToSingleMolSDF "WEIGHT_DEFINITION_3" "$WEIGHT_DEFINITION_3" "$f" 
    addPropertyToSingleMolSDF "WEIGHT_DEFINITION_4" "$WEIGHT_DEFINITION_4" "$f" 
    addPropertyToSingleMolSDF "FITNESS_DEFINITION" "$FITNESS_DEFINITION" "$f"
    addPropertyToSingleMolSDF "DESCRIPTOR_1" "$desc1" "$f"
    addPropertyToSingleMolSDF "DESCRIPTOR_2" "$desc2" "$f"
    addPropertyToSingleMolSDF "DESCRIPTOR_3" "$desc3" "$f" 
    addPropertyToSingleMolSDF "WEIGHT_1" "$w1" "$f" 
    addPropertyToSingleMolSDF "WEIGHT_2" "$w2" "$f" 
    addPropertyToSingleMolSDF "WEIGHT_3" "$w3" "$f" 
    addPropertyToSingleMolSDF "WEIGHT_4" "$w4" "$f" 
    addPropertyToSingleMolSDF "FITNESS" "$fitness" "$f"
    echo -e "$f successfully edited."
done
