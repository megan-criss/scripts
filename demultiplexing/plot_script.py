#!/usr/bin/env python

import matplotlib.pyplot as plt
import argparse
import gzip

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

def convert_phred(char):
    qscore = (ord(char))-33
    return qscore

def init_list(value=0.0):
    lst=[]
    for values in range(read_length):
        lst.append(value)
        return(lst)

def fill_list(file):
    p_scores = init_list()
    with gzip.open(file, "rt") as fh:  
        LN=0
        for line in fh:
            LN+=1
            linestr = line.strip()
            if LN%4 == 0:
                counter=0
                for char in linestr:
                    score = convert_phred(char) 
                    p_scores[counter] += score
                    counter+=1
    records = LN/4
    return p_scores, records

means, num_records = fill_list(rfile)

for char in range(len(means)):
    means[char] = means[char] / num_records

indexes = []
for x in range(len(means)):
    indexes.append(x)


plt.figure(figsize=[12,6])
#plot axes
plt.bar(indexes, means)
plt.xlabel( "Position")
plt.ylabel("Mean Quality Score")
plt.title("Mean Quality Scores Per Position")
plt.savefig(plot_name)
