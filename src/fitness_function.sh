#!/bin/bash
#
# This file collects functions that are needed to compute the fitness and write
# it into an SDF file with minimal dependencies (effectively only bash commands)
# for maximal portability.
#
# @author Marco Foscato
# @author Jonas Brattebø Ekeli
#

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
# Function that takes the raw energy values and computes the
# fitness according to the hard-coded formulation. The raw energies, the fitness
# score, and many other details defining the calcualtion of the fitness score
# are reported (i.e., appended as properties) to the SDF file given in among
# the command-line arguments.
#
# @param $1 pathname of SDF file where to save the information about the
#        calculation of the fitness (and the fitness itself).
# @param $2 energy value for TS of productive metathesis
# @param $3 energy value for TS of bH elimination
# @param $4 energy value for MCB intermediate (DFT-modelled)
# @param $5 energy value for catalyst precursor
# @param $6 energy value for MCB intermediate (xTB-modelled)
# @param $7 energy value for catalyst recursor with Cl replacing both X ligands
# @param $8 energy value for catalyst recursor with cis-anionic configuration
# @param $9 energy value for free ligand L
#
function computeFitness() {
    # General use constants
    hartree_to_kcalmol="627.49467516" # (kcal/mol*hartree)
    boltzmann_constant="0.0000000000000000000000138066" # (J/K)
    avogadro_number="602214000000000000000000" # (1/mol)
    jmol_to_kcalmol="0.00023901" # (kcal/J)
    temp="313.15" # (K) ( 40 degrees C to match laboratory experiments)
    kT=$( echo "$boltzmann_constant * $avogadro_number * $jmol_to_kcalmol * $temp" | bc -l ) # kcal/mol

    # Reference values in Hartree
    G_Propene="-117.738020895"
    G_Ethene="-78.4716324751"
    G_HoveydaProd="-502.204109499"
    G_SIMes="-924.27249048"
    G_PCy3="-1046.11440274"
    G_HG_RuCl2_SIMes="-2402.13728523"
    G_HI_precursor="-2523.95659027"
    DG_referenceProductionBarrier=".0222529298"
    DDG_reference_HGII="0.0222312922"
    DG_referenceProductionBarrier="0.0222714248"
    DG_referencePrecursorStabilityHII="-0.0052201621"
    DG_referencePrecursorStabilityHI="0.0117775312"
    DDG_referencePrecursorStability="$( echo "$DG_referencePrecursorStabilityHI - $DG_referencePrecursorStabilityHII" | bc -l )"

    # Magnitude<n> represents n=1,2,3 orders of magnitude increase of rate constant at 40 C according to the Eyring equation ( kcal/mol )
    magnitude3="4.29872427"
    magnitude2="2.86581618"
    magnitude1="1.43290809"
    magnitude4="5.73163236"

    # Parse CLI arguments
    if [ "$#" -ne 9 ]
    then
        echo "ERROR(Fitness function): wrong number of arguments ($#)"
        return 1
    fi
    sdfFile="$1"
    if [ ! -f "$sdfFile" ]
    then
      echo "ERROR(Fitness function): file '$sdfFile' not found!"
      return -1
    fi

    # DFT energies on DFT geometries
    freeEnergyX="$2" # TS productive metathesis
    freeEnergyZ="$3" # TS bH elimination
    freeEnergyD="$4" # MCB intermediate (DFT-modelled)
    # DFT energies on xTB geometries
    freeEnergyA="$5" # Precursor
    freeEnergyC="$6" # MCB intermediate (xTB-modelled)
    freeEnergyE="$7" # Precursor with Cl replacing both X ligands
    freeEnergyF="$8" # Precursor with cis-anionic configuration
    freeEnergyL="$9" # Free ligand L


    # Descriptor 1 -Boltzmann factor: (Z) - (X) - DDG_ref(HG2).
    coef1="1"
    DDG=$( echo "( ( $freeEnergyZ ) - ( $freeEnergyX + 2*$G_Ethene - 2*$G_Propene ) )" | bc -l )
    desc1=$( echo "$coef1 * ( e( $hartree_to_kcalmol * ( $DDG - $DDG_reference_HGII ) / $kT ) ) " | bc -l )
    DESCRIPTOR_DEFINITION_1='desc1=( echo "1 * ( e( $hartree_to_kcalmol * ( ( $freeEnergyZ ) - ( $freeEnergyX + 2*$G_Ethene - 2*$G_Propene ) - $DDG_reference_HGII ) / $kT ) ) " | bc -l )'


    # Descriptor 2 - Linear: (Z) - (X) (or 0 if lower than 0)
    coef2="1"
    desc2=$( echo "$coef2 * ( $hartree_to_kcalmol * ( $DDG ))" | bc -l )
    if [ "$( echo "$DDG < 0.0" | bc -l )" == "1" ]
    then
        desc2="0.0000000001"
    fi
    DESCRIPTOR_DEFINITION_2='desc2=$( echo "1 * ( $hartree_to_kcalmol * ( ( $freeEnergyZ ) - ( $freeEnergyX + 2*$G_Ethene - 2*$G_Propene ) ) )" | bc -l )'


    # Descriptor 3 - linear: Production Barrier relative to HG Ru SIMes.
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
    DG_synt=$( echo "$hartree_to_kcalmol * ( $freeEnergyE + $G_PCy3 - $G_HI_precursor - $freeEnergyL )" | bc -l )
    w3=$( echo "1 / ( 1 + e( ( $DG_synt - $magnitude1 ) / $kT ) )" | bc -l )
    WEIGHT_DEFINITION_3='w3=$( echo "1 / ( 1 + e( $hartree_to_kcalmol * ( ( $freeEnergyE + $G_PCy3 - $G_HI_precursor - $freeEnergyL ) - $magnitude1 ) / $kT ) )" | bc -l )'


    # Dynamic weight 4 (sigmoid): disfavouring non trans precursors.
    threshold7="$magnitude2"
    DG_stereo=$( echo "$hartree_to_kcalmol * ( $freeEnergyF - $freeEnergyA )" | bc -l )
    w4=$( echo "1 / ( 1 + e( ( - $DG_stereo + $threshold7 ) / $kT ) )" | bc -l )
    WEIGHT_DEFINITION_4='w4=$( echo "1 / ( 1 + e( $hartree_to_kcalmol * ( - ( $freeEnergyF - $freeEnergyA ) + $magnitude2 ) / $kT ) )" | bc -l )'

    #FITNESS (Sum of fitness from descriptos 1-3 weighted by dynamic weights 1-4).
    fitness=$( echo "( $desc1 + $desc2 + $desc3 ) * $w1 * $w2 * $w3 * $w4 " | bc -l )

    removePropertyFromSDF "DESCRIPTOR_1" "$sdfFile"
    removePropertyFromSDF "DESCRIPTOR_DEFINITION_1" "$sdfFile"
    removePropertyFromSDF "DESCRIPTOR_2" "$sdfFile"
    removePropertyFromSDF "DESCRIPTOR_DEFINITION_2" "$sdfFile"
    removePropertyFromSDF "DESCRIPTOR_3" "$sdfFile"
    removePropertyFromSDF "DESCRIPTOR_DEFINITION_3" "$sdfFile"
    removePropertyFromSDF "WEIGHT_1" "$sdfFile"
    removePropertyFromSDF "WEIGHT_DEFINITION_1" "$sdfFile"
    removePropertyFromSDF "WEIGHT_2" "$sdfFile"
    removePropertyFromSDF "WEIGHT_DEFINITION_2" "$sdfFile"
    removePropertyFromSDF "WEIGHT_3" "$sdfFile"
    removePropertyFromSDF "WEIGHT_DEFINITION_3" "$sdfFile"
    removePropertyFromSDF "WEIGHT_4" "$sdfFile"
    removePropertyFromSDF "WEIGHT_DEFINITION_4" "$sdfFile"
    removePropertyFromSDF "FITNESS" "$sdfFile"
    removePropertyFromSDF "freeEnergyX" "$sdfFile"
    removePropertyFromSDF "freeEnergyZ" "$sdfFile"
    removePropertyFromSDF "freeEnergyD" "$sdfFile"
    removePropertyFromSDF "freeEnergyA" "$sdfFile"
    removePropertyFromSDF "freeEnergyC" "$sdfFile"
    removePropertyFromSDF "freeEnergyE" "$sdfFile"
    removePropertyFromSDF "freeEnergyF" "$sdfFile"
    removePropertyFromSDF "freeEnergyL" "$sdfFile"
    removePropertyFromSDF "hartree_to_kcalmol" "$sdfFile"
    removePropertyFromSDF "boltzmann_constant" "$sdfFile"
    removePropertyFromSDF "avogadro_number" "$sdfFile"
    removePropertyFromSDF "jmol_to_kcalmol" "$sdfFile"
    removePropertyFromSDF "temp" "$sdfFile"
    removePropertyFromSDF "kT" "$sdfFile"
    removePropertyFromSDF "G_Propene" "$sdfFile"
    removePropertyFromSDF "G_Ethene" "$sdfFile"
    removePropertyFromSDF "G_HoveydaProd" "$sdfFile"
    removePropertyFromSDF "G_SIMes" "$sdfFile"
    removePropertyFromSDF "G_PCy3" "$sdfFile"
    removePropertyFromSDF "G_HG_RuCl2_SIMes" "$sdfFile"
    removePropertyFromSDF "G_HI_precursor" "$sdfFile"
    removePropertyFromSDF "DG_referenceProductionBarrier" "$sdfFile"
    removePropertyFromSDF "DDG_reference_HGII" "$sdfFile"
    removePropertyFromSDF "DG_referenceProductionBarrier" "$sdfFile"
    removePropertyFromSDF "DG_referencePrecursorStabilityHII" "$sdfFile"
    removePropertyFromSDF "DG_referencePrecursorStabilityHI" "$sdfFile"
    removePropertyFromSDF "DDG_referencePrecursorStability" "$sdfFile"
    removePropertyFromSDF "magnitude3" "$sdfFile"
    removePropertyFromSDF "magnitude2" "$sdfFile"
    removePropertyFromSDF "magnitude1" "$sdfFile"
    removePropertyFromSDF "magnitude4" "$sdfFile"
    addPropertyToSingleMolSDF "hartree_to_kcalmol" "$hartree_to_kcalmol" "$sdfFile"
    addPropertyToSingleMolSDF "boltzmann_constant" "$boltzmann_constant" "$sdfFile"
    addPropertyToSingleMolSDF "avogadro_number" "$avogadro_number" "$sdfFile"
    addPropertyToSingleMolSDF "jmol_to_kcalmol" "$jmol_to_kcalmol" "$sdfFile"
    addPropertyToSingleMolSDF "temp" "$temp" "$sdfFile"
    addPropertyToSingleMolSDF "kT" "$kT" "$sdfFile"
    addPropertyToSingleMolSDF "G_Propene" "$G_Propene" "$sdfFile"
    addPropertyToSingleMolSDF "G_Ethene" "$G_Ethene" "$sdfFile"
    addPropertyToSingleMolSDF "G_HoveydaProd" "$G_HoveydaProd" "$sdfFile"
    addPropertyToSingleMolSDF "G_SIMes" "$G_SIMes" "$sdfFile"
    addPropertyToSingleMolSDF "G_PCy3" "$G_PCy3" "$sdfFile"
    addPropertyToSingleMolSDF "G_HG_RuCl2_SIMes" "$G_HG_RuCl2_SIMes" "$sdfFile"
    addPropertyToSingleMolSDF "G_HI_precursor" "$G_HI_precursor" "$sdfFile"
    addPropertyToSingleMolSDF "DG_referenceProductionBarrier" "$DG_referenceProductionBarrier" "$sdfFile"
    addPropertyToSingleMolSDF "DDG_reference_HGII" "$DDG_reference_HGII" "$sdfFile"
    addPropertyToSingleMolSDF "DG_referenceProductionBarrier" "$DG_referenceProductionBarrier" "$sdfFile"
    addPropertyToSingleMolSDF "DG_referencePrecursorStabilityHII" "$DG_referencePrecursorStabilityHII" "$sdfFile"
    addPropertyToSingleMolSDF "DG_referencePrecursorStabilityHI" "$DG_referencePrecursorStabilityHI" "$sdfFile"
    addPropertyToSingleMolSDF "DDG_referencePrecursorStability" "$DDG_referencePrecursorStability" "$sdfFile"
    addPropertyToSingleMolSDF "magnitude3" "$magnitude3" "$sdfFile"
    addPropertyToSingleMolSDF "magnitude2" "$magnitude2" "$sdfFile"
    addPropertyToSingleMolSDF "magnitude1" "$magnitude1" "$sdfFile"
    addPropertyToSingleMolSDF "magnitude4" "$magnitude4" "$sdfFile"
    addPropertyToSingleMolSDF "freeEnergyX" "$freeEnergyX" "$sdfFile"
    addPropertyToSingleMolSDF "freeEnergyZ" "$freeEnergyZ" "$sdfFile"
    addPropertyToSingleMolSDF "freeEnergyD" "$freeEnergyD" "$sdfFile"
    addPropertyToSingleMolSDF "freeEnergyA" "$freeEnergyA" "$sdfFile"
    addPropertyToSingleMolSDF "freeEnergyC" "$freeEnergyC" "$sdfFile"
    addPropertyToSingleMolSDF "freeEnergyE" "$freeEnergyE" "$sdfFile"
    addPropertyToSingleMolSDF "freeEnergyF" "$freeEnergyF" "$sdfFile"
    addPropertyToSingleMolSDF "freeEnergyL" "$freeEnergyL" "$sdfFile"
    addPropertyToSingleMolSDF "DESCRIPTOR_1" "$desc1" "$sdfFile"
    addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_1" "$DESCRIPTOR_DEFINITION_1" "$sdfFile"
    addPropertyToSingleMolSDF "DESCRIPTOR_2" "$desc2" "$sdfFile"
    addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_2" "$DESCRIPTOR_DEFINITION_2" "$sdfFile"
    addPropertyToSingleMolSDF "DESCRIPTOR_3" "$desc3" "$sdfFile"
    addPropertyToSingleMolSDF "DESCRIPTOR_DEFINITION_3" "$DESCRIPTOR_DEFINITION_3" "$sdfFile"
    addPropertyToSingleMolSDF "WEIGHT_1" "$w1" "$sdfFile"
    addPropertyToSingleMolSDF "WEIGHT_DEFINITION_1" "$WEIGHT_DEFINITION_1" "$sdfFile"
    addPropertyToSingleMolSDF "WEIGHT_2" "$w2" "$sdfFile"
    addPropertyToSingleMolSDF "WEIGHT_DEFINITION_2" "$WEIGHT_DEFINITION_2" "$sdfFile"
    addPropertyToSingleMolSDF "WEIGHT_3" "$w3" "$sdfFile"
    addPropertyToSingleMolSDF "WEIGHT_DEFINITION_3" "$WEIGHT_DEFINITION_3" "$sdfFile"
    addPropertyToSingleMolSDF "WEIGHT_4" "$w4" "$sdfFile"
    addPropertyToSingleMolSDF "WEIGHT_DEFINITION_4" "$WEIGHT_DEFINITION_4" "$sdfFile"
    addPropertyToSingleMolSDF "CALCULATION_OF_OVERALL_FITNESS" 'fitness=$( echo "( $desc1 + $desc2 + $desc3 ) * $w1 * $w2 * $w3 * $w4 " | bc -l )' "$sdfFile"
    addPropertyToSingleMolSDF "FITNESS" "$fitness" "$sdfFile"
    return 0
}
