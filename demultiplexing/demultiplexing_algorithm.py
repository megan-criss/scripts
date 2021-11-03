#!/usr/bin/env python

#importing necessary things
import argparse
import itertools
import re
import gzip

#defining necessary arguments
def get_args():
    parser = argparse.ArgumentParser(description = " A program to demultiplex FASTQ files")
    parser.add_argument("-r1", "--read1", help ="Insert the name of the first read file you wish to use", required = True)
    parser.add_argument("-r2", "--read2", help = "Insert the name of the second read file you wish to use", required = True)
    parser.add_argument("-i1", "--index1", help = "Insert the name of the first index file you wish to use", required = True)
    parser.add_argument("-i2", "--index2", help = "Insert the name of the second index file you wish to use", required = True)
    parser.add_argument("-idx", "--indexfile", help = "Insert the name of the file containing the indexes", required = True)
    parser.add_argument("-q", '--qcutoff', help = 'Enter the Q-cutoff score you want to use', default= 30, type = int)
    return parser.parse_args()


args = get_args()
R1= args.read1
R2 = args.read2
I1 = args.index1
I2 = args.index2
idx = args.indexfile
q = args.qcutoff

#defining function that will both convert phred score and help sort due to quality of read
def phred_check(phred_line):
    """Converts a single character into a phred score"""
    for letter in phred_line:
        qualityscore = (ord(letter)) - 33
        if qualityscore < q:
            return False
    return True

#print(phred_check('  P3%kjsjf'))


def reverse_complement(phred_line):
#A function to find the reverse complement of a sequence of DNA and return it, even if there are unkown bases'''
    comp_DNA=""
    for bp in phred_line:
        if bp == "T":
            bp_comp = "A"
        elif bp == "G": 
            bp_comp="C"
        elif bp == "A":
            bp_comp = "T"
        elif bp == "C":
            bp_comp = "G"
        else:
            bp_comp = bp
        comp_DNA = comp_DNA+bp_comp
        rev_comp = comp_DNA[::-1]
    return rev_comp

# line1= "ATG"
# line2= "NCG"
# print(reverse_complement(line1))
# print(reverse_complement(line2))


#creating a dictionary of indexes, and reverse comlpements of indexes, where keys are indexes and values are additional information
def make_index_dictionary(idx):
    idx_dictionary = {}
    with open(idx) as idxfile:
        index_counts = 0
        for index in idxfile:
            if index_counts == 0:
                index_counts += 1
            else:
                index = index.strip() 
                column = index.split("\t")
                #keys(idx seq):values(sample, group, treatment, idx ID)
                idx_key = column[-1]
                idx_value = column[0:4]
                idx_dictionary[idx_key] = (idx_value)
        # for key, value in idx_dictionary.items():
        #     print(key, value)
        return idx_dictionary

#calling index dictionary maker to make dictionary
index_dictionary = make_index_dictionary(idx)
#print(index_dictionary)

#making a list of all indexes
total_indexes=[]
for key in index_dictionary.keys():
    #print(key)
    total_indexes.append(key)
#print(total_indexes)


reverse_complement_dict = {}
for key in index_dictionary:
    reverse_complement_dict[reverse_complement(key)] = ('')
#print(reverse_complement_dict)


#creates a dictionary of permutations for all indexes present for index-hopped files
#Keys are permutations and values will be counts
def perm_dict_maker(pair_of_indexes):
    permutations_dict = {}
    permutations = itertools.product(index_dictionary, repeat = 2)
    for x in permutations:
        index1, index2 = x
        permutations_dict.setdefault(index1 + '-' + index2, 0)
    return permutations_dict

#calls permutation dictionary maker to create a dictionary from list of all indexes
permutation_dict=perm_dict_maker(total_indexes)
#print(permutation_dict.keys())

#setting filenames for matched files that are to be made


matched_files = {} 
#print(index_dictionary.values())
for idx_seq, info in index_dictionary.items():
    fileR1 = info[3] + "_" + "matched" + '_' + info[1] + "_" + info[2] + '_' + info[0] + '_' + "_R1.fastq"
    fileR2 = info[3] + "_" + "matched" + '_' + info[1] + "_" + info[2] + '_' + info[0] + '_' + "_R2.fastq"
    matched_files.setdefault(idx_seq, [fileR1, fileR2])
#print(matched_names)

    
#opening matched files to write information to
for key in matched_files:
    for x in range(0,2):
        matched_files[key][x] = open(matched_files[key][x], 'w')

#print(matched_files.keys())

##opening index hopped files/unmatched files to write to
index_hop_R1 = open('index_hopped_R1.fq', 'w')
index_hop_R2 = open('index_hopped_R2.fq', 'w')

#opening files for unknown indicies and low-quality indicies to write to
unknown_badqual_R1 = open('unknown_lowqual_R1.fq', 'w')
unknown_badqual_R2 = open('unknown_lowqual_R2.fq', 'w')



#starting main code, wish me luck
#def demultiplex_files(permutation_dict, index_dictionary, reverse_complement_dict, r1, r2, i1, i2):
#setting counters
ln = 0
rec_numb = 0
match_count = 0
unmatched_count = 0
unknown_count = 0
badquality_count = 0
    #setting lists to fill with records
rec_r1 = []
rec_r2 = []
rec_i1 = []
rec_i2 = []
    ###opening files of info to be demultiplexed
    ###remember to add gzip in
r1 = gzip.open(R1, "rt")
r2 = gzip.open(R2, 'rt') 
i1 = gzip.open(I1, "rt") 
i2 = gzip.open(I2, 'rt')
for line in zip(r1, r2, i1, i2):
    ln += 1 
    r1_line_stripped = line[0].strip()
    r2_line_stripped = line[1].strip()
    i1_line_stripped = line[2].strip()
    i2_line_stripped = line[3].strip()
    if ln // 4 == rec_numb:
        rec_r1.append(r1_line_stripped)
        rec_r2.append(r2_line_stripped)
        rec_i1.append(i1_line_stripped)
        rec_i2.append(i2_line_stripped)
    else:
        rec_r1.append(r1_line_stripped)
        rec_r2.append(r2_line_stripped)
        rec_i1.append(i1_line_stripped)
        rec_i2.append(i2_line_stripped)
        rec_numb += 1
        #print(rec_i1, '\n' ,sep = '\n')

        pair_of_indexes = rec_i1[1] + '-' + reverse_complement(rec_i2[1])
        #print(rec_i1[1])
        #print(matched_files.keys())
        #first check index quality scores
        #print(phred_check(rec_i1[3]))
        #print(rec_i1[1])
        #print(reverse_complement(rec_i2[1]))
        if phred_check(rec_i1[3]) is True and phred_check(rec_i2[3]) is True:
            #if quality scores are good, check if the indexes are real
            if rec_i1[1] in index_dictionary.keys() and reverse_complement(rec_i2[1]) in index_dictionary.keys():
                #making sure reverse complements are real keys
                if reverse_complement(rec_i1[1]) in reverse_complement_dict.keys() and rec_i2[1] in reverse_complement_dict.keys():
                    #if indexes real, and quality scores are good, check if indexes match
                    permutation_dict[pair_of_indexes] += 1
                    if rec_i1[1] == reverse_complement(rec_i2[1]):
                        for idx_seq in matched_files:
                            if idx_seq == rec_i1[1]:
                                matched_files[idx_seq][0].write(rec_r1[0] + "\t" + pair_of_indexes + '\n')
                                matched_files[idx_seq][0].write(rec_r1[1] + '\n')
                                matched_files[idx_seq][0].write(rec_r1[2] + '\n')
                                matched_files[idx_seq][0].write(rec_r1[3] + '\n')
                                matched_files[idx_seq][1].write(rec_r2[0] + "\t" + pair_of_indexes + '\n')
                                matched_files[idx_seq][1].write(rec_r2[1] + '\n')
                                matched_files[idx_seq][1].write(rec_r2[2] + '\n')
                                matched_files[idx_seq][1].write(rec_r2[3] + '\n')
                                match_count += 1
                    else:
                        #indexes are real, quality scores are good, but indexes don't match
                        permutation_dict[pair_of_indexes] += 1 
                        index_hop_R1.write(rec_r1[0] + "\t" + pair_of_indexes + '\n')
                        index_hop_R1.write(rec_r1[1] + '\n')
                        index_hop_R1.write(rec_r1[2] + '\n')
                        index_hop_R1.write(rec_r1[3] + '\n')
                        index_hop_R2.write(rec_r1[0] + "\t" + pair_of_indexes + '\n')
                        index_hop_R2.write(rec_r1[1] + '\n')
                        index_hop_R2.write(rec_r1[2] + '\n')
                        index_hop_R2.write(rec_r1[3] + '\n')
                        unmatched_count += 1
                else:
                    #quality scores are good but indexes are unknown
                    unknown_badqual_R1.write(rec_r1[0] + "\t" +  pair_of_indexes + '\n')
                    unknown_badqual_R1.write(rec_r1[1] + '\n')
                    unknown_badqual_R1.write(rec_r1[2] + '\n')
                    unknown_badqual_R1.write(rec_r1[3] + '\n')
                    unknown_badqual_R2.write(rec_r1[0] + "\t" + pair_of_indexes + '\n')
                    unknown_badqual_R2.write(rec_r1[1] + '\n')
                    unknown_badqual_R2.write(rec_r1[2] + '\n')
                    unknown_badqual_R2.write(rec_r1[3] + '\n')
                    unknown_count += 1
            else:
                #quality scores are good but indexes are unknown
                unknown_badqual_R1.write(rec_r1[0] + "\t" +  pair_of_indexes + '\n')
                unknown_badqual_R1.write(rec_r1[1] + '\n')
                unknown_badqual_R1.write(rec_r1[2] + '\n')
                unknown_badqual_R1.write(rec_r1[3] + '\n')
                unknown_badqual_R2.write(rec_r1[0] + "\t" + pair_of_indexes + '\n')
                unknown_badqual_R2.write(rec_r1[1] + '\n')
                unknown_badqual_R2.write(rec_r1[2] + '\n')
                unknown_badqual_R2.write(rec_r1[3] + '\n')
                unknown_count += 1
        else:
            #quality scores are bad
            unknown_badqual_R1.write(rec_r1[0] + "\t" + pair_of_indexes + '\n')
            unknown_badqual_R1.write(rec_r1[1] + '\n')
            unknown_badqual_R1.write(rec_r1[2] + '\n')
            unknown_badqual_R1.write(rec_r1[3] + '\n')
            unknown_badqual_R2.write(rec_r1[0] + "\t" + pair_of_indexes + '\n')
            unknown_badqual_R2.write(rec_r1[1] + '\n')
            unknown_badqual_R2.write(rec_r1[2] + '\n')
            unknown_badqual_R2.write(rec_r1[3] + '\n')
            badquality_count += 1

        rec_r1=[]
        rec_r2=[]
        rec_i1=[]
        rec_i2=[]

r1.close()
r2.close()
i1.close()
i2.close()

#print(index_dictionary)
#print(permutation_dict)

#creating stats file
stats = open('stats_file.txt', 'w')
number_matched = match_count
number_hopped = unmatched_count
number_unknown = unknown_count
number_lowqual = badquality_count
number_records = rec_numb

stats.write('Overall counts:'+ '\n')
stats.write("Total number of records:"+ ' ' + str(number_records) + '\n')
stats.write("Total number of records with matching indices:"+ '\t' + str(number_matched) + '\n')
stats.write("Total number of records with hopped indices:"+ '\t' + str(number_hopped) + '\n')
stats.write("Total number of records with unknown indicies:"+ '\t' + str(number_unknown) + '\n')
stats.write("Total number of records with bad quality scores:"+ '\t' + str(number_lowqual) + '\n')

stats.write('\n' +'Overall percentages:'+ '\n')
stats.write("Percent of records with matching indicies:" + '\t' + str(((number_matched/number_records)*100)) + "%" + '\n')
stats.write("Percent of records with index hopping:" + '\t' + str(((number_hopped/number_records)*100)) + '%' + '\n')
stats.write("Percent of records with unknown indicies:" + '\t' + str(((number_unknown/number_records)*100)) + '%' + '\n')
stats.write("Percent of records with quality scores below threshhold:" + '\t' + str(((number_lowqual/number_records)*100)) + '%' + '\n')

stats.write('\n' + "percentage of reads per sample:" + '\n')
for index in index_dictionary:
    for x in range(4):
        stats.write(str(index_dictionary[index][x]) + '\t')
    matching_indexes = (str(index) + "-" + str(index))
    for pairs in permutation_dict:
        if matching_indexes == pairs:
            stats.write(str((permutation_dict[pairs]/number_records)*100))
    stats.write('\n')

stats.write('\n' + "Permutation counts:" + '\n')
for key, value in permutation_dict.items():
    stats.write(str(key) + '\t' + str(value) + '\n')

for key in matched_files:
    for x in range(0,2):
        matched_files[key][x].close()

index_hop_R1.close()
index_hop_R2.close()
unknown_badqual_R1.close()
unknown_badqual_R2.close()
stats.close()

