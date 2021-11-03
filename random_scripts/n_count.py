#!/usr/bin/env python

import argparse
import re
import os.path

def get_args():
    parser = argparse.ArgumentParser(description = " A program to check a consensus fasta for N-content")
    parser.add_argument("-i", "--input", help ="Input consensus fasta file", required = True)
    parser.add_argument("-s", "--save_path", help ="The path you want files saved to(default: current directory)", required = False, default= "./")
    return parser.parse_args()

args = get_args()
input_file = args.input
save_path = args.save_path

with open(input_file, "r") as inputfile:
    '''Small function that finds the name of the file and counts the percentage of N's in a multiline FASTA file'''
    n_count = 0
    bp_count = 0
    for line in inputfile:
        line = line.strip()
        if line.startswith(">"):
            name_components = line.split(re.search("\_(?!.*\_)", string = line)[0])
            map_number = str(name_components[0])
            map_number = str(map_number.split(">")[1])
            special_number = str(name_components[1])
            output_name = (map_number + "_" + special_number)
        else:
            for character in line:
                if character == "N":
                    n_count += 1
                    bp_count += 1
                else:
                    bp_count += 1
    n_percentage = (n_count / bp_count)*100

output_file = os.path.join(save_path, (output_name + ".csv"))
with open (output_file, "w") as output_file:
    output_file.write("sample,percent_N" + "\n")
    output_file.write(output_name + "," + str(round(n_percentage, 1)) + "%")


