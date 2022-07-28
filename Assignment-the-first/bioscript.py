#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import bioinfo
import gzip
import argparse

# argparse to get the file name, # of records and the read length
def get_args():
    #defines variables
    parser = argparse.ArgumentParser(description="The variables to initialize the dictionary")
    parser.add_argument("-n", help="Number of records", type=int)
    parser.add_argument("-l", help="The length of the reads", type=int)
    parser.add_argument("-f", help="The file name")
    parser.add_argument("-o", help="Output file name")
    return parser.parse_args()
args=get_args()

my_list = [0]*args.l

# with gzip.open(file,"rt") as fh:
#     x = len(fh.readlines())
#     x = int(x)/4
# some_array = np.zeros((args.l,args.n),dtype=np.int8)

def pop_list(file: str) -> tuple[list, int]:
    with gzip.open(args.f,"rt") as fh:
        count = 0
        for line in fh:
            line = line.strip()
            count += 1
            if count % 4 == 0:
                for i, base in enumerate(line):
                    my_list[i] += bioinfo.convert_phred(base)
                    # some_array[i, int(count/4)-1] += bioinfo.convert_phred(base)
    return (my_list, count)

my_list, num_lines = pop_list(args.f)

for i, item in enumerate(my_list):
    my_list[i] = (item/(num_lines/4))

fig, ax = plt.subplots(figsize=(10, 6))
plt.bar(range(args.l), my_list, align='edge', color='#007030')
plt.xlabel('Base Position')
plt.ylabel('Mean Quality Score')
plt.title('Mean Quality Scores per Base Position')
plt.savefig(args.o)