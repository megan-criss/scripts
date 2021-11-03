#!/usr/bin/env python

import argparse
import re
import os.path


def get_args():
    parser = argparse.ArgumentParser(description = " A program to deduplicate a SAM file (check for PCR duplicates)")
    parser.add_argument("-f", "--file", help ="SAM file", required = True)
    parser.add_argument("-m", "--MAPnumber", help ="The MAP testing number", required = True)
    parser.add_argument("-s", "--save_path", help ="The path you want files saved to(default: current directory)", required = False, default= "./")
    return parser.parse_args()

args = get_args()
input_file = args.file
map_num = args.MAPnumber
save_path = args.save_path

def set_position(line):
    sam_line_strip = line.strip()
    elements = sam_line_strip.split("\t")
    cigar_string = elements[5]
    original_position = int(elements[3])
    bit_flag = int(elements[1])
    if ((bit_flag&2) == 0):                 #determines if properly paired or not, if == 0 then not mapped in proper pair
        boolean = False
        return boolean, "None"
    if ((bit_flag&2) != 0): 
        boolean = True
    ######################################## forward strand ##################################################################
        if ((bit_flag&16) == 0): #this means 0x10 is unchecked, means it is the forward strand
            ########################### R1 ####################################################################
            if ((bit_flag&64) != 0): #this means 0x40 is checked, R1 
                ####################### First in Pair #########################################################
                if ((bit_flag&128) == 0): #this means 0x80 is unchecked, first segment in template
                    if "S" in cigar_string:
                        clip_position = int(re.search('[0-9][0-9]*', string = cigar_string)[0])
                        clip_letter = re.search('[A-Z]', string = cigar_string)[0]
                        if clip_letter == "S":
                            new_position = abs(original_position - clip_position)
                            #print("A")
                            return boolean, new_position
                        else:
                            new_position = abs(original_position)
                            #print("B")
                            return boolean, new_position
                    else:
                        new_position = abs(original_position)
                        #print("C")
                        return boolean, new_position
            #################### R2 ############################################################
            elif ((bit_flag&64) != 0): #this means 0x40 is checked, R1
                #################### first in pair ############################################################
                if ((bit_flag&128) == 0): #this means 0x80 is UNchecked, FIRST segment in template
                    if "S" in cigar_string:
                        clip_position = int(re.search('[0-9][0-9]*', string = cigar_string)[0])
                        clip_letter = re.search('[A-Z]', string = cigar_string)[0]
                        if clip_letter == "S":
                            new_position = abs(original_position - clip_position)
                            #print("A")
                            return boolean, new_position
                        else:
                            new_position = abs(original_position)
                            #print("B")
                            return boolean, new_position
                    else:
                        new_position = abs(original_position)
                        #print("C")
                        return boolean, new_position         
                else:
                    new_position = "None"
                    return boolean, new_position
            else:
                new_position = "None"
                return boolean, new_position           
    ########################### Reverse Strand  #################################################################
        elif ((bit_flag&16) != 0): #this means 0x10 is checked, means it is the reverse strand
            ################### R2 #############################################################################
            if ((bit_flag&64) == 0): #this means 0x40 is unchecked, R2
                #################### Second in pair ############################################################
                if ((bit_flag&128) != 0): #this means 0x80 is checked, second segment in template
                    if "S" in cigar_string:
                        clip_letter = re.search('[A-Z]', string = cigar_string)[0]
                        if clip_letter == "S":
                            new_string = cigar_string.split(re.search('[A-Z]', string = cigar_string)[0])[1]
                            cigar_values = re.findall("[0-9][0-9]*", new_string)
                            cigar_ints = [int(i) for i in cigar_values]
                            cigar_sum = sum(cigar_ints)
                            new_position = abs(original_position + cigar_sum)
                            return boolean, new_position
                        else:
                            cigar_values = re.findall("[0-9][0-9]*", cigar_string)
                            cigar_ints = [int(i) for i in cigar_values]
                            cigar_sum = sum(cigar_ints)
                            new_position = abs(original_position + cigar_sum)
                            return boolean, new_position
                    else:
                        cigar_values = re.findall("[0-9][0-9]*", cigar_string)
                        cigar_ints = [int(i) for i in cigar_values]
                        cigar_sum = sum(cigar_ints)
                        new_position = abs(original_position + cigar_sum)
                        #print("O")
                        return boolean, new_position
                else:
                    new_position = "None"
                    return boolean, new_position
            #################### R1 ############################################################
            elif ((bit_flag&64) == 0): #this means 0x40 is ununchecked, R2
                #################### first in pair ############################################################
                if ((bit_flag&128) != 0): #this means 0x80 is checked, second segment in template
                    if "S" in cigar_string:
                        clip_letter = re.search('[A-Z]', string = cigar_string)[0]
                        if clip_letter == "S":
                            new_string = cigar_string.split(re.search('[A-Z]', string = cigar_string)[0])[1]
                            cigar_values = re.findall("[0-9][0-9]*", new_string)
                            cigar_ints = [int(i) for i in cigar_values]
                            cigar_sum = sum(cigar_ints)
                            new_position = abs(original_position + cigar_sum)
                            return boolean, new_position
                        else:
                            cigar_values = re.findall("[0-9][0-9]*", cigar_string)
                            cigar_ints = [int(i) for i in cigar_values]
                            cigar_sum = sum(cigar_ints)
                            new_position = abs(original_position + cigar_sum)
                            return boolean, new_position
                    else:
                        cigar_values = re.findall("[0-9][0-9]*", cigar_string)
                        cigar_ints = [int(i) for i in cigar_values]
                        cigar_sum = sum(cigar_ints)
                        new_position = abs(original_position + cigar_sum)
                        #print("O")
                        return boolean, new_position
                else:
                    new_position = "None"
                    #print("E")
                    return boolean, new_position    
            else:
                new_position = "None"
                #print("E")
                return boolean, new_position 
        else:
            new_position = "None"
            return boolean, new_position

#writes output to file named via the input map number

paired_file = os.path.join(save_path, ("intermediate_" + map_num + '.sam'))
unpaired_file = os.path.join(save_path, ("unpaired_reads_" + map_num + '.sam'))
error_file = os.path.join(save_path, ("error_reads_" + map_num + '.sam'))


with open(input_file, "r") as inputfile, open(paired_file, "w") as paired_file, open(unpaired_file, "w") as unpaired_file, open(error_file, "w") as error_file:
    for line in inputfile:
            stripped_read = line.strip()
            if line.startswith("@"):
                error_file.write(line)
            else:
                return_quantities = set_position(line)
                boolean = return_quantities[0]
                position_read = return_quantities[1]
                if boolean == False:
                    unpaired_file.write(str(stripped_read)+ "\n")
                elif boolean == True:    
                    if position_read == "None":
                        error_file.write(stripped_read + "\n")
                    else:
                        full_read = str(stripped_read)
                        header = full_read.split('\t', maxsplit = 1)[0]
                        rest_of_line = full_read.split('\t', maxsplit = 1)[1]
                        paired_file.write(header + "@" + str(position_read) + '\t' + rest_of_line + "\n")
