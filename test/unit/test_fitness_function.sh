#!/bin/bash
#
# Here we test the job submissin script for AutoCompChem jobs
#

function mustHaveString()
{
    str="$1"
    file="$2"
    refLineNo="$3"
    if ! grep -q "$str" "$file" ; then
        echo "ERROR! Test failed at '$refLineNo' of '$0'"
        exit -1
    fi
}

function mustNotHaveString()
{
    str="$1"
    file="$2"
    refLineNo="$3"
    if grep -q "$str" "$file" ; then
        echo "ERROR! Test failed at '$refLineNo' of '$0'"
        exit -1
    fi
}

function valueIsDistantFromZero()
{
    valueIsDistantFromReference "$1" '0.000' "$2" "$3"
}

function valueIsDistantFromReference()
{
    propertyName="$1"
    reference="$2"
    file="$3"
    refLineNo="$4"
    v=$(grep -A 1 "<$propertyName>" "$file" | tail -n 1)
    diff=$(echo "$v - $reference" | bc -l)
    if (( $(echo "$diff < 0 " | bc -l) )) ; then
        diff=$(echo "-1.0 * $diff " | bc -l)
    fi
    if (( $(echo "$diff > 0.001" | bc -l) )); then
        echo "ERROR! Test failed at '$refLineNo' of '$0'"
        exit -1
    fi
}


log=test_fitness_function.log
tmpFile=dummy_tmp.sdf

rm -f "$log" "$tmpFile" 2> /dev/null

cp dummy.sdf "$tmpFile"

source ../../src/fitness_function.sh

if [[ $(type -t computeFitness ) != function ]]
then
  echo "Not passed: function 'computeFitness' not defined!"
fi

computeFitness non-existing-file 2 3 4 5 6 7 8 9 > "$log"
mustHaveString 'file.*not found' "$log" $LINENO

computeFitness "$tmpFile" > "$log"
mustHaveString 'wrong number of arguments' "$log" $LINENO

computeFitness "$tmpFile" 2 3 4 5 6 7 8 > "$log"
mustHaveString 'wrong number of arguments' "$log" $LINENO

x="-2134.22385645"
z="-2055.67028667"
d="-2055.71274937"
a="-2400.94412028"
c="-2055.68864266"
e="-2400.94412028"
f="-2400.92847469"
l="-923.082037925"

computeFitness "$tmpFile" "$x" "$z" "$d" "$a" "$c" "$e" "$f" "$l" > "$log"
mustHaveString "${x//-/\\-}" "$tmpFile" $LINENO
mustHaveString "${z//-/\\-}" "$tmpFile" $LINENO
mustHaveString "${d//-/\\-}" "$tmpFile" $LINENO
mustHaveString "${a//-/\\-}" "$tmpFile" $LINENO
mustHaveString "${c//-/\\-}" "$tmpFile" $LINENO
mustHaveString "${e//-/\\-}" "$tmpFile" $LINENO
mustHaveString "${f//-/\\-}" "$tmpFile" $LINENO
mustHaveString "${l//-/\\-}" "$tmpFile" $LINENO
valueIsDistantFromReference 'DESCRIPTOR_1' '0.23449108' "$tmpFile" $LINENO
valueIsDistantFromReference 'DESCRIPTOR_2' '13.0474593' "$tmpFile" $LINENO
valueIsDistantFromReference 'DESCRIPTOR_3' '0.37754158' "$tmpFile" $LINENO
valueIsDistantFromReference 'WEIGHT_1' '1.0' "$tmpFile" $LINENO
valueIsDistantFromReference 'WEIGHT_2' '1.0' "$tmpFile" $LINENO
valueIsDistantFromReference 'WEIGHT_3' '1.0' "$tmpFile" $LINENO
valueIsDistantFromReference 'WEIGHT_4' '1.0' "$tmpFile" $LINENO
valueIsDistantFromReference 'FITNESS' '13.65855' "$tmpFile" $LINENO
mustNotHaveString 'WRONG' "$tmpFile" $LINENO
mustHaveString 'not_touched_value' "$tmpFile" $LINENO

computeFitness "$tmpFile" "$x" "$z" '-2055.80' "$a" "$c" "$e" "$f" "$l" > "$log"
valueIsDistantFromZero 'FITNESS' "$tmpFile" $LINENO
valueIsDistantFromZero 'WEIGHT_1' "$tmpFile" $LINENO

computeFitness "$tmpFile" "$x" "$z" "$d" "$a" '-2055.75' "$e" "$f" "$l" > "$log"
valueIsDistantFromZero 'FITNESS' "$tmpFile" $LINENO
valueIsDistantFromZero 'WEIGHT_2' "$tmpFile" $LINENO

computeFitness "$tmpFile" "$x" "$z" "$d" "$a" '-2055.60' "$e" "$f" "$l" > "$log"
valueIsDistantFromZero 'FITNESS' "$tmpFile" $LINENO
valueIsDistantFromZero 'WEIGHT_2' "$tmpFile" $LINENO

computeFitness "$tmpFile" "$x" "$z" "$d" "$a" "$c" '-2400.90' "$f" "$l" > "$log"
valueIsDistantFromZero 'FITNESS' "$tmpFile" $LINENO
valueIsDistantFromZero 'WEIGHT_3' "$tmpFile" $LINENO

computeFitness "$tmpFile" "$x" "$z" "$d" "$a" "$c" "$e" '-2400.95' "$l" > "$log"
valueIsDistantFromZero 'FITNESS' "$tmpFile" $LINENO
valueIsDistantFromZero 'WEIGHT_4' "$tmpFile" $LINENO

rm -f "$log" "$tmpFile"

echo "Test passed!"
