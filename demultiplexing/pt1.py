#!/usr/bin/env python
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=meanscores     ### Job Name
#SBATCH --output=meanscores_%j.out.out       ### File in which to store job output
#SBATCH --error=meanscores_%j.out.err       ### File in which to store job error messages
#SBATCH --time=0-04:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1       ### Number of tasks to be launched
#SBATCH --account=bgmp          ### Account used for job submission

import gzip

file = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
#do not change this path - move your file to this location
#file = "./test.fastq"

def convert_phred(char):
    """Converts a single character into a phred score"""
    return (ord(char))-33



def init_list(array, value=0.0):
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    # I think this line is unneccessary....
    array=[]
    for values in range(102):
        array.append(value)
    #print(array)
    return(array)
mean_scores=[]
mean_scores=init_list(mean_scores)
real_means=[]
real_means=init_list(real_means)

with gzip.open(file, "r") as fh:
#starting to separate lines and strip them of \n character     
    LN=0
    for line in fh:
        LN+=1
        linestr = line.strip()
        if LN%4 == 0:
            counter=0
            #converting phred scores and creating mean_scores
            for char in linestr:
                score=convert_phred(char) 
                mean_scores[counter] += score
                counter+=1
    #LN is now the total number of lines in the file! 
#starting loop to calculate means of the Quality scores
            records=LN/4
            for char in range(102):
                real_means[char]=mean_scores[char]/records
    


#Use this cell to generate your plot
import matplotlib.pyplot as plt
fig=plt.figure()
plt.style.use('ggplot')

#plot axes
ax=fig.add_axes([0,0,1,1])
plt.xlabel( "Base Pair Number", size=40)
plt.ylabel("Mean Quality Score", size=40)
plt.title("Mean Quality Scores per Read", size=40)

# plt.rc('ytick', labelsize=40)
# plt.rc('xtick', labelsize=40)

#data to plot
x=(list(range(0,102)))
y=(real_means)

#setting axes
ax.bar(x,y)
plt.savefig("test_plot.png")



file.close()
