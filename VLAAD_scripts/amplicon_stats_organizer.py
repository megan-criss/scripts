#!/usr/bin/env python

import argparse
import os.path
import re
###############################################################
def get_args():
    parser = argparse.ArgumentParser(description = " A program to find motifs within a fasta file")
    parser.add_argument("-f", "--file", help ="Amplicon Stats file", required = True)
    parser.add_argument("-s", "--save_path", help ="The path you want files saved to(default: current directory)", required = False, default= "./")
    parser.add_argument("-m", "--MAPnumber", help ="The MAP testing number", required = True, type = str )
    return parser.parse_args()
args = get_args()
amp_file = args.file
save_path = args.save_path
map_num = args.MAPnumber
###############################################################
out_file = os.path.join(save_path, (map_num + "_ampliconstats_organized" + '.txt'))
###############################################################
def line_getter1(line):
    new_line = line.split(re.search("(.*?\t)", string = line)[1])
    return new_line
def line_getter2(line):
    new_line = line.split(re.search("(FSS\t)", string = line)[1])
    return new_line
def line_getter3(line):
    new_line = line.split(re.search("(FREADS\t)", string = line)[1])
    return new_line
def line_getter4(line):
    new_line = line.split(re.search("(FDEPTH\t)", string = line)[1])
    return new_line
def line_getter5(line):
    new_line = line.split(re.search("(FTCOORD\t)", string = line)[1])
    return new_line
###############################################################
def rule_0(line):
    outfile.write( "\n" + "\n" + "Amplicon Statistics for " + map_num + ":" + "\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+ "\n")
    outfile.write( "Summary Statistics (used for scaling the plot)" + "\n" + "----------------------------------------------" + "\n")
    return()
def rule_1(line):
    outfile.write("\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    outfile.write("\n" + "Amplicon locations from BED file" + "\n" + "(LEFT/RIGHT are <start>-<end> format and comma-separated for alt-primers)" + "\n" + "\n")
    outfile.write("REFERENCE    AMPLICON #   LEFT/RIGHT Position" + "\n" + "---------------------------------------------------------------------" + "\n")
    return()
def rule_2(line):
    outfile.write("\n" + "\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    outfile.write("Summary Statistics" + "\n")
    outfile.write( "------------------------------------" + "\n")
    return()
def rule_3(line):
    outfile.write("\n" + "\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    outfile.write("Number of reads per amplicon" + "\n" + "\n")
    outfile.write("Amplicon" + "\t" + "Number of Reads" + "\n")
    outfile.write( "--------------------------------------" + "\n")
    return()
def rule_4(line):
    outfile.write("\n" + "\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    outfile.write("Average read depth per amplicon" + "\n" + "\n")
    outfile.write("Amplicon" + "\t" + "Depth of Reads" + "\n")
    outfile.write( "--------------------------------------" + "\n")
    return()
def rule_5(line):
    outfile.write("\n" + "\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    outfile.write("\n" + "Template start,end coordinate frequencies per amplicon" + "\n" )
    outfile.write("(Status key: 0 for OK, 1 for skipping amplicon, 2 for unknown location)" + "\n" + "\n")
    outfile.write("Amplicon" + "\t" + "start, end, frequency, and status" + "\n")
    outfile.write( "-----------------------------------------------------" + "\n")
    return()
###############################################################
amp_file = open(amp_file, "r")
outfile = open(out_file, "w")
###############################################################
for line in amp_file:
    line = line.strip()
    #### printing summary statistics #1 (Amplicon and file counts. Always comes first)
    if line.startswith("# Summary statistics, used for scaling the plots."):
        rule_0(line)
    elif line.startswith("SS"):
        ss_line = line_getter1(line)
        outfile.write(ss_line[1] + "\n")
    #### printing Amplicon primer locations
    elif line.startswith("# AMPLICON	REF	NUMBER	LEFT	RIGHT"):
        rule_1(line)
    elif line.startswith("AMPLICON"):
        amplicon_line = line_getter1(line)
        outfile.write(amplicon_line[1] + "\n")
    #### printing File specific: summary stats
    elif line.startswith("# Use 'grep ^FSS | cut -f 2-' to extract this part."):
        rule_2(line)
    elif line.startswith("FSS"):
        fss_line = line_getter2(line)
        outfile.write(fss_line[1] + "\n")
    #### printing File specific: numbers of reads per amplicon
    elif line.startswith("# Use 'grep ^FREADS | cut -f 2-' to extract this part."):
        rule_3(line)
    elif line.startswith("FREADS"):
        freads_line = line_getter3(line)[1]
        counter = 0
        freads_line_split = list(freads_line.split("\t"))
        for item  in freads_line_split:
            mytuple = ("MAP", "Twist", "NoRT")
            if item.startswith(mytuple):
                pass
            else:
                counter += 1
                outfile.write(str(counter) + "\t" + item + "\n")
    #### Printing File specific: average read depth per amplicon
    elif line.startswith("FDEPTH"):
        rule_4(line)
        fdepth_line = line_getter4(line)[1]
        counter = 0
        fdepth_line_split = list(fdepth_line.split("\t"))
        for item  in fdepth_line_split:
            mytuple = ("MAP", "Twist", "NoRT")
            if item.startswith(mytuple):
                pass
            else:
                counter += 1
                outfile.write(str(counter) + "\t" + item + "\n")
    #### printing File specific: template start,end coordinate frequencies per amplicon
    elif line.startswith("# Use 'grep ^FTCOORD | cut -f 2-' to extract this part."):
        rule_5(line)
    elif line.startswith("FTCOORD"):
        ftcoord_line = line_getter5(line)[1]
        outfile.write(ftcoord_line + "\n")
    else:
        pass
###############################################################
outfile.write("\n" + "Finished writing Amplicon-stats summary, have a nice day :)" + "\n")
###############################################################
amp_file.close()
outfile.close()