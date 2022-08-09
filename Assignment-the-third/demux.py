#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import bioinfo
import gzip
import argparse
import math
from itertools import permutations

barcode_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

# argparse to get the file name, # of records and the read length
def get_args():
    #defines variables
    parser = argparse.ArgumentParser(description="The variables to initialize the dictionary")
    parser.add_argument("-r1", help="read one")
    parser.add_argument("-r2", help="read two")
    parser.add_argument("-i1", help="index one")
    parser.add_argument("-i2", help="index two")
    parser.add_argument("-k", help="known barcodes list")
    parser.add_argument("-o", help="output plot")
    return parser.parse_args()
args=get_args()

r1 = gzip.open(args.r1,"rt")
r2 = gzip.open(args.r2, "rt")
i1 = gzip.open(args.i1, "rt")
i2 = gzip.open(args.i2, "rt")
k = open(args.k, "rt")

#writes files that contain each of these things
r1_unknown = open("r1_unknown.fq", "w")
r2_unknown = open("r2_unknown.fq", "w")
r1_hopped = open("r1_hopped.fq", "w")
r2_hopped = open("r2_hopped.fq", "w")

#gets the barcodes from the indexes.txt and puts them into a set
known_codes = set()
with open(args.k, 'r') as k:
    lines=k.readlines()[1:]
    for x in lines:
        known_codes.add((x).split('\t')[4].strip())

#writes the file names based on barcodes
file_dictionary = {}
for index in known_codes:
    file_dictionary[index] = [open(f'output/{index}_r1_matched.fq', 'w'), open(f'output/{index}_r2_matched.fq', 'w')]

#creates a dictionary with all permutations of the barcodes and sets the count to 0
hop_dictionary = {}
p = permutations(known_codes,2)
for j in list(p):
    hop_dictionary[j] = 0

#creates a dictionary with all matching sets of barcodes and sets the count to 0
match_dictionary = {}
for x in known_codes:
    match_dictionary[x, x] = 0

#sets the count to 0
unknown = 0

while True:
    #reads and labels header lines
    r1_header = r1.readline().strip()
    if r1_header == "":
        break
    r2_header = r2.readline().strip()
    i1_header = i1.readline().strip()
    i2_header = i2.readline().strip()
    #splits each file into their individual lines
    r1_array = np.array([r1_header, r1.readline().strip(), r1.readline().strip(), r1.readline().strip()])
    r2_array = np.array([r2_header, r2.readline().strip(), r2.readline().strip(), r2.readline().strip()])
    i1_array = np.array([i1_header, i1.readline().strip(), i1.readline().strip(), i1.readline().strip()])
    i2_array = np.array([i2_header, i2.readline().strip(), i2.readline().strip(), i2.readline().strip()])

    rev = bioinfo.rev_comp(i2_array[1])

#writes all the string concatenation stuff so I don't have to clog up the following code
    fq_r1_string = r1_header + ('_') + i1_array[1] + ('_') + rev \
            + ('\n') + r1_array[1] + ('\n') + r1_array[2] + ('\n') + r1_array[3] + ('\n')
    fq_r2_string = r2_header + ('_') + i1_array[1] + ('_') + rev \
            + ('\n') + r2_array[1] + ('\n') + r2_array[2] + ('\n') + r2_array[3] + ('\n')

#quality control, checks if phred in index quality score is lower than 30 and puts it in unknown if so
    written = False
    for letter in i1_array[3]:
        phred = bioinfo.convert_phred(letter)
        if int(phred) < 30:
            r1_unknown.write(fq_r1_string)
            r2_unknown.write(fq_r2_string)
            unknown += 1
            written = True
            break
        elif letter in i2_array[3]:
            phred = bioinfo.convert_phred(letter)
            if int(phred) < 30:
                r1_unknown.write(fq_r1_string)
                r2_unknown.write(fq_r2_string)
                unknown += 1
                written = True
                break
    if written == False:
        #checks if the reverse complement of the barcodes in i2 are in i1 and writes to a file
        if i1_array[1] in known_codes:
            if rev == i1_array[1]:
                file_dictionary[i1_array[1]][0].write(fq_r1_string)
                file_dictionary[i1_array[1]][1].write(fq_r2_string)
                match_dictionary[(i1_array[1]), (i1_array[1])] += 1
            else:
                #if the reverse comp of i2 is in the known barcodes, put in hopped
                if rev in known_codes:
                    r1_hopped.write(fq_r1_string)
                    r2_hopped.write(fq_r2_string)
                    hop_dictionary[(i1_array[1]), rev] += 1
                #if not in barcodes, put in unknown
                else: 
                    r1_unknown.write(fq_r1_string)
                    r2_unknown.write(fq_r2_string)
                    unknown += 1
                    
        #if not the reverse comp, put in unknown
        else:
            r1_unknown.write(fq_r1_string)
            r2_unknown.write(fq_r2_string)
            unknown += 1


Total_Records = sum(match_dictionary.values()) + sum(hop_dictionary.values()) + unknown

#write a loop for the matched dictionary's output
print("Matched Barcodes")
for a in match_dictionary:
    print(a, match_dictionary[a], sep=('\t'))

#write a loop for the hopped dictionary's output
print("Hopped Barcodes" + ('\n'))
for b in hop_dictionary:
    print(b, hop_dictionary[b], sep=('\t'))

print("------------------------------------------------------------")

#total records
print("Total Matched")
print(sum(match_dictionary.values()))
print("Total Hopped")
print(sum(hop_dictionary.values()))
print("Total Unknown")
print(unknown)
print("Total Records")
print(Total_Records)

Percentage_of_Matched = sum(match_dictionary.values())/Total_Records * 100
Percentage_of_Hopped = sum(hop_dictionary.values())/Total_Records * 100
Percentage_of_Unknown = unknown/Total_Records * 100

print("Percentage of Matched")
print(Percentage_of_Matched)
print("Percentage of Hopped")
print(Percentage_of_Hopped)
print("Percentage of Unknown")
print(Percentage_of_Unknown)

#makes a pie chart
labels = 'Matched', 'Hopped', 'Unknown'
sizes = [Percentage_of_Matched, Percentage_of_Hopped, Percentage_of_Unknown]
fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=labels, autopct='%1.1f%%', colors=('#002243', '#003F6C', '#006AE4'))
ax1.axis('equal')
plt.title('Proportion of Matched, Hopped, and Unknown Reads')
plt.savefig(args.o)