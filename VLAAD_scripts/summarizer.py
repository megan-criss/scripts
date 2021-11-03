#!/usr/bin/env python

from re import S
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os.path

#argparse expressions to take in user input
def get_args():
    parser = argparse.ArgumentParser(description = " A program to find motifs within a fasta file")
    parser.add_argument("-f", "--file", help ="SAM file", required = True)
    parser.add_argument("-m", "--MAPnumber", help ="The MAP testing number", required = True)
    parser.add_argument("-s", "--save_path", help ="The path you want files saved to(default: current directory)", required = False, default= "./")
    parser.add_argument("-i", "--intermediate_file", help ="The file that has been sorted and filtered for only properly aligned pairs (used for fragment length", required = True)
    parser.add_argument("-l", "--fragment_length", help ="The max fragment length that will still be included and not filtered out (default: 500)", required = False, default= 500, type = int)
    parser.add_argument("-p", "--plots", help ="Do you want plots of fragment length? (True or False, default = False)", required = False, default = False, type = bool)
    return parser.parse_args()
    
args = get_args()
sam_file = args.file
map_num = args.MAPnumber
save_path = args.save_path
intermediate_file = args.intermediate_file
max_frag_length = args.fragment_length
want_plot = args.plots

#writes output to file named via the input map number
out_file = os.path.join(save_path, (map_num + "_fragment_statistics" + '.txt'))
filtered_sam = os.path.join(save_path, (map_num + "_fragment_filtered_sam" + '.sam'))
big_reads = os.path.join(save_path, (map_num + "_fragments_too_big" + '.sam'))
small_reads = os.path.join(save_path, (map_num + "_fragments_too_small" + '.sam'))

fragment_dict = dict()
fragment_list = []
fragment_zero_list = []
good_reads = open(filtered_sam, "w")
big_reads = open(big_reads, 'w')
small_reads = open(small_reads, "w")

def fragment_length(intermediate_file):
    '''determines fragment length using the intermediate file created using position appender script'''
    frag_count = 0
    while True:
            line1 = intermediate_file.readline()
            line2 = intermediate_file.readline()
            if(line2 == ""):
                break
            header1 = line1.split("\t")[0]
            name1 = (header1.split("@")[0])
            pos_read = int(header1.split("@")[1])
            sam_line_1 = line1.split("\t", 1)[1]
            header2 = line2.split("\t")[0]
            name2 = (header2.split("@")[0])
            pos_pair = int(header2.split("@")[1])
            sam_line_2 = line2.split("\t", 1)[1]
            if name1 != name2:
                print(name1, name2)
                raise Exception("The read names do not match, cannot continue")
            else: 
                frag_length = abs( pos_pair - pos_read )
                if frag_length == 0:
                    frag_count += 1
                    fragment_zero_list.append(str(name1) + ":" + str(name2))
                elif frag_length >= int(max_frag_length):
                    sam_reads = (name1 + "\t" + sam_line_1 + name2 + "\t" + sam_line_2 )
                    fragment_list.append(frag_length)
                    big_reads.write(sam_reads)
                elif frag_length < int(max_frag_length):
                    sam_reads = (name1 + "\t" + sam_line_1 + name2 + "\t" + sam_line_2 )
                    fragment_list.append(frag_length)
                    good_reads.write(sam_reads)
                    if frag_length <= 200:
                        small_reads.write(sam_reads)
                    else:
                        pass
# this decides the length between read pairs/ how far away pairs are from one another

with open(sam_file,'r') as sam_file, open(out_file, "w") as out_file, open(intermediate_file, "r") as intermediate_file:
    unpaired_count = 0
    paired_count = 0
    LN = 0
    print(map_num)
    for line in sam_file:
        LN += 1
        if LN % 200000 == 0:
            print("200,000 lines parsed through")
        if line.startswith("@"):
            good_reads.write(line)
            big_reads.write(line)
        else:
            #checks if read is paired or single
            flag=int(line.split()[1])
            # will check if the read's pair is mapped or unmapped
            if((flag & 2) == 2):
                paired_count += 1
            else:
                unpaired_count += 1
    fragment_length(intermediate_file)
    ###########################################################################
    out_file.write('Summary statistics for ' + str(map_num) + "\n" + "-----" + "\n")
    out_file.write("Number of reads aligned as pairs: " + str(paired_count) + "\n" )
    out_file.write("Number of reads aligned as singles: " + str(unpaired_count) + "\n")
    out_file.write("Number of fragments with length of zero: " + str(len(fragment_zero_list)) + "\n")
    out_file.write("\n" + "-----------------------------------------------------" + "\n")
    ###########################################################################
    out_file.write("Fragment length Stats " + "\n" + "-----" + "\n")
    out_file.write("Maximum Fragment Length: " + str(max(fragment_list)) + "\n")
    out_file.write("Minimum Fragment Length: " + str(min(fragment_list)) + "\n")
    out_file.write("Average Fragment Length: " + str(np.average(fragment_list)) + "\n")
    out_file.write("Median Fragment Length: " + str(np.median(fragment_list)) + "\n" + "-----------------------------------------------------" + "\n")
    out_file.write("Fragment Length , Counts" + "\n" + "-----" + "\n")
    ###########################################################################
    fragment_count_dict = dict()
    for i in fragment_list:
        fragment_count_dict[i] = fragment_count_dict.get(i, 0)+1
    fragment_count_list = list(fragment_count_dict.items())
    fragment_count_list.sort()
    ###########################################################################
    for item in fragment_count_list:
        out_file.write(str(item) + "\n")
    out_file.write("-----------------------------------------------------" + "\n""Reads with fragment length of zero:" + "\n" + "-----" + "\n")
    for item in fragment_zero_list:
        out_file.write(str(item) + "\n")

##############################################################################
##############################################################################
##############################################################################
##############################################################################

if want_plot == True:
    fig=plt.figure()
    plt.style.use('ggplot')
    #plot axes
    plt.xlabel( "Fragment length (limit 1000 bp)")
    plt.ylabel("Number of fragments in this category")
    plt.title("Bar Chart of Fragment Length vs. Quantity for " + str(map_num))
    plt.xlim([0,1000])
    '''plt.rcParams['figure.figsize'] = (20-9,9)'''
    #data to plot
    x=(fragment_count_dict.keys())
    y=(fragment_count_dict.values())
    plt.bar(list(x),list(y), color = 'C2')
    plt.show()
    plt.savefig(os.path.join(save_path, (map_num + "_fragment_distribution")))
    plt.close()
    plt.cla()
    ##############################################################################
    fig=plt.figure()
    plt.style.use('ggplot')
    #plot axes
    plt.xlabel( "Fragment lengths zoomed in (limit 500 bp)")
    plt.ylabel("Number of fragments in this category")
    plt.title("Bar Chart of Fragment Length vs. Quantity for " + str(map_num))
    plt.xlim([0,500])
    '''plt.rcParams['figure.figsize'] = (20-9,9)'''
    #data to plot
    x=(fragment_count_dict.keys())
    y=(fragment_count_dict.values())
    plt.bar(list(x),list(y), color = 'C2')
    plt.show()
    plt.savefig(os.path.join(save_path, (map_num + "_zoomed_fragment_distribution")))
    plt.close()
    plt.cla()
    #####################################################################################
    print("Finished plots and summary statistics, have a nice day :)")
else:
    pass

