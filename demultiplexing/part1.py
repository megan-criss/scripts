#!/usr/bin/env python

import matplotlib.pyplot as plt
import argparse
import gzip
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description = " A program to open a zipped file and create a plot to show average phred score per position")
    parser.add_argument("-f", "--file", help ="Zipped fastq", required = True, type = str)
    parser.add_argument("-l", "--readlength", help = "untrimmed read length", type = int, default = 101 )
    parser.add_argument("-p", "--pngname", help = "create name for plot png", type = str, required = True)
    return parser.parse_args()


args = get_args()
rfile= args.file
read_length = args.readlength
plot_name = args.pngname




#file = gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz")


def convert_phred(char):
    """Converts a single ASCII character into a phred score"""
    qscore = (ord(char))-33
    return qscore


def init_list(value=0.0):
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    # I think this line is unneccessary....
    lst=[]
    for values in range(read_length):
        lst.append(value)
    #print(array)
    return(lst)


def fill_list(file):
    """Filling list with phred scores for base position"""
    p_scores = init_list()
    with gzip.open(file, "rt") as fh:
#starting to separate lines and strip them of \n character     
        LN=0
        for line in fh:
            LN+=1
            linestr = line.strip()
            if LN%4 == 0:
                counter=0
                #converting phred scores and creating mean_scores
                for char in linestr:
                    score = convert_phred(char) 
                    p_scores[counter] += score
                    counter+=1
        #LN is now the total number of lines in the file! 
        #starting loop to calculate means of the Quality scores
    records = LN/4
    return p_scores, records

#populating lists
means, num_records = fill_list(rfile)

# divide in place 
for char in range(len(means)):
    means[char] = means[char] / num_records

#indexes needed to plot
indexes = []
for x in range(len(means)):
    indexes.append(x)


# plotting the mean scores
plt.figure(figsize=[12,6])
#plot axes
plt.bar(indexes, means)
plt.xlabel( "Position", size=40)
plt.ylabel("Mean Quality Score", size=40)
plt.title("Mean Quality Scores at Per Position", size=40)
plt.savefig(plot_name)
