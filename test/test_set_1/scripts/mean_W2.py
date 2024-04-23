#!/usr/bin/env python
import math
import sys
from sys import argv

def calculate_mean(numbers_str):

    # Convert the string of numbers to a list of floats
    numbers = [float(x) for x in numbers_str.split() if float(x) != 0]
  
    # Check if there are non-zero numbers to avoid division by zero
    if len(numbers) == 0:
        return 0

    # Calculate the mean
    mean = sum(numbers) / len(numbers)

    return mean

fitness = ''
for file in argv[1:]:
	readfile = open(file,'r')
	for line in readfile:
		if '> <WEIGHT_2>' in line:
			thisline = readfile.read()
			fitness+=thisline[:7] + ' '
#print('Fitness values: ', fitness)
result = round(calculate_mean(fitness), 2)
print(result)
