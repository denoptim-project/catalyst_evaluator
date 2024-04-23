#!/usr/bin/env python
import math
import sys
from sys import argv

def calculate_standard_deviation(numbers_str):
    # Convert the string of numbers to a list of floats
    numbers = [float(x) for x in numbers_str.split()]

    # Calculate the mean
    mean = sum(numbers) / len(numbers)

    # Calculate the squared differences from the mean
    squared_diff = [(x - mean) ** 2 for x in numbers]

    # Calculate the variance
    variance = sum(squared_diff) / ( len(numbers) - 1 )

    # Calculate the standard deviation
    std_deviation = math.sqrt(variance)

    return std_deviation

fitness = ''
for file in argv[1:]:
	readfile = open(file,'r')
	for line in readfile:
		if '> <WEIGHT_4>' in line:
			thisline = readfile.read()
			fitness+=thisline[:7] + ' '
#print('Fitness values: ', fitness)
result = round(calculate_standard_deviation(fitness), 2)
print(result)
