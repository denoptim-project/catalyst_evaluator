#!/usr/bin/env python
import glob
import os
import re
import pandas as pd
from optparse import OptionParser


def calculate_mean(numbers_str):
    numbers = [float(x) for x in numbers_str.split()]
    mean = sum(numbers) / len(numbers)
    return mean


def is_float(string: any) -> bool:
    if string is None:
        return False
    try:
        float(string)
        return True
    except ValueError:
        return False


def extract_properties_from_sdf(pathname):
    if not os.path.isfile(pathname):
        print("ERROR: file '%s' is missing!" % pathname)
        quit(1)
    pattern = re.compile("^> *<.*>$")
    all_values = {}
    with open(pathname, "r") as input_file:
        lines = input_file.readlines()
        for index, line in enumerate(lines):
            if pattern.match(line):
                property_name = re.split('<|>', line)[2]
                value_str = lines[index + 1].rstrip()
                if is_float(value_str):
                    all_values[property_name] = float(value_str)
                else:
                    all_values[property_name] = value_str
        print(f"Parsed {len(all_values)} properties from {pathname}")
    return all_values


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f", type="string", dest="pathnames_list",
                      help="The pathname/s to the file/s to read, which are "
                           "expected to be SDF filed produced by DENOPTIM's "
                           "fitness provider (e.g., *_out.sdf). Use comma as "
                           "separator.")
    parser.add_option("-r", type="string", dest="file_finding_regex",
                      help="The regex to be used to find files to read, "
                           "which are expected to be SDF filed produced by "
                           "DENOPTIM's fitness provider (e.g., *_out.sdf). ")
    parser.add_option("-p", type="string", dest="property_to_analyze",
                      help="The string identifying the SDF property for which "
                           "to show basics statistics.",
                      default='FITNESS')

    (options, args) = parser.parse_args()
    if not options.pathnames_list and not options.file_finding_regex:
        parser.error('ERROR: No list of pathnames of regex given!')
    pathnames = []
    if options.pathnames_list:
        pathnames = options.pathnames_list.split(',')
    if options.file_finding_regex:
        regex = str(options.file_finding_regex).strip('\'')
        pathnames.extend(sorted(glob.glob(regex)))

    rows = []

    for pathname in pathnames:
        rows.append(extract_properties_from_sdf(pathname))


    df = pd.DataFrame(rows)
    if len(df.columns.tolist()) > 0:
        print('Retrieved properties:', df.columns.tolist())

    pd.set_option("display.precision", 8)

    if options.property_to_analyze:
        column_name = options.property_to_analyze.strip('\'')
        print(f'Values of {column_name}')
        print(df.loc[:, column_name])
        print('Statistics')
        print(df.agg({column_name: ['count',
                                    'min',
                                    'max',
                                    'mean',
                                    'std',
                                    'median',
                                    'skew']}))

