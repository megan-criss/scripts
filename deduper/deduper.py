#!/usr/bin/env python

import argparse
import gzip
import re

def get_args():
    parser = argparse.ArgumentParser(description = " A program to deduplicate a SAM file (check for PCR duplicates)")
    parser.add_argument("-f", "--file", help ="SAM file", required = True)
    parser.add_argument("-p", "--paired", help = "Makes sure that file is not paired end (not single-end), type y or n", type = bool, required = True)
    parser.add_argument("-u", "--umi", help = "designates file containing the list of UMIs (unset if randomers instead of UMIs)", required = True)
    return parser.parse_args()

args = get_args()
input_file = args.file
UMIs = args.umi 
paired = args.paired

if paired == "True":
    print("Use an un-paired file for this program you silly goose")
    exit()

#import os
#os.system('samtools sort %s -o sam_sorted.sam' %input_file)

#input_file = "/home/mcriss/bioinformatics/Bi610/deduper/Deduper/test_file.sam"
#UMIs = '/home/mcriss/bioinformatics/Bi610/deduper/Deduper/STL96.txt'

with open(UMIs, "r") as UMIs:
    umi_list=[]
    LN = 0
    for line in UMIs:
        line = line.strip()
        #print(line)
        umi_list.append(line)
        LN += 1

#sam_line = "NS500451:154:HWKTMBGXX:1:11101:24260:1121:TCGTAGGT	0	2	76814284	36	71M12M24N32M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU"

def real_umi(sam_line):
    ''' A function that is designed to check an UMI against a stored list of UMI's provided by the user to make sure its real!'''
    sam_line_strip = sam_line.strip("\n")
    elements = sam_line.split("\t")
    ident = elements[0]
    identifier = ident.split(":")
    umi = identifier[-1]
    if umi in umi_list:
        return(True)
    else: 
        error.write(sam_line_strip + '\n')
        return(False)

def set_position(sam_line):
    '''if (-) strand add to position column value, if (+) strand subtract from position column value'''
    sam_line_strip = sam_line.strip()
    elements = sam_line_strip.split("\t")
    cigar_string = elements[5]
    og_position = int(elements[3])
    bit_flag = int(elements[1])
    if ((bit_flag&16)==0):
        strand = "forward"
        if "S" in cigar_string:
            #is there is soft clipping in this cigar string then it searches for the S, takes the value before it and subtracts it from the position (leftmost soft clipping)
            clip_pos = re.search('[0-9][0-9]*', string = cigar_string)
            clip_let = re.search('[A-Z]', string = cigar_string)
            clip_letter = clip_let[0]
            clip_position = int(clip_pos[0])
            if "S" in clip_letter:
                new_position = og_position - clip_position
                return new_position, strand
            else:
                new_position = og_position
                return new_position, strand
        else:
            new_position = og_position
            return new_position, strand
    elif ((bit_flag&16) == 16): 
        strand = "reverse"
        #this portion is for the reverse strand, i.e. you need to see if there is soft clipping on the end of the cigar string- if the strand is reverse then all other values in the cigar string must be added to the position to find the actual starting position
        if "S" in cigar_string:
            clip_pos = re.search("[0-9][0-9]*", cigar_string)
            clip_let = re.search("[A-Z]", cigar_string)
            clip_letter=clip_let[0]
            if "S" in clip_letter:
                cigar_vals = cigar_string.split("S")
                new_cigar = cigar_vals[1]
                cigar_values = re.findall("[0-9][0-9]*", new_cigar)
                cigar_ints = [int(i) for i in cigar_values]
                cigar_sum = sum(cigar_ints)
                new_position = og_position + cigar_sum
                return new_position, strand
            else:
                cigar_values = re.findall("[0-9][0-9]*", new_cigar)
                cigar_ints = [int(i) for i in cigar_values]
                cigar_sum = sum(cigar_ints)
                new_position = og_position + cigar_sum
                return new_position, strand
        else:
            cigar_values = re.findall("[0-9][0-9]*", cigar_string)
            cigar_ints = [int(i) for i in cigar_values]
            cigar_sum = sum(cigar_ints)
            new_position = og_position + cigar_sum
            return new_position, strand
    else:
        duplicates.write(sam_line_strip + '\n')
            #position is equal to all other values besides s plus the original position

umi_dict = {} 
for i in umi_list: 
    umi_dict[i] = None

#Start the actualbody of code here:
chromosome_dict = {}
#umi_counter = 0
total_reads = 0

##opening files for kept reads/unkept reads to write to
keep = open('kept_reads.sam', 'w')
duplicates = open('duplicate_reads.sam', 'w')
error =  open('error_reads.sam', 'w')

read_chromosome = str(0)
LN = 0
with open(input_file, "r") as inputfile:
    #go chromosome by chromosome: chrom number is column 3
    for sam_line in inputfile:
        stripped_read = sam_line.strip()
        if sam_line.startswith("@"):
            pass 
        else:            #real_umi(sam_line)
            if real_umi(sam_line) == False:
                pass
            else:
                full_read = str(stripped_read)
                elements = stripped_read.split("\t")
                chromosome = elements[2]
                strand = elements[1]
                if read_chromosome != chromosome:
                    chromosome_dict = {}
                    read_chromosome = chromosome
                    print("cleared ditionary for new chromosome")
                if read_chromosome == chromosome:
                    x = set_position(sam_line)
                    umi = elements[0].split(':')[-1]
                    umi = str(umi)
                    position = str(x[0])
                    chromosome = str(chromosome)
                    strand = str(x[1])
                    header_info = (umi + ":" + strand + ":" + chromosome + ":" + position)
                    if header_info in chromosome_dict.keys():
                        duplicates.write(full_read + '\n')
                        LN+=1
                        if LN == 1:
                                print(full_read)
                    else: 
                        chromosome_dict[header_info] = 0
                        keep.write(full_read + "\n")

keep.close()
duplicates.close()
error.close()
UMIs.close()
