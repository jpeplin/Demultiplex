Lab Notebook
Jack Peplinski
7/26/22

This is my lab notebook for the demultiplexing assignments. 

## ---- Assignment the first ----

## Part 1 – Quality Score Distribution per-nucleotide

1.
Initial Data Exploration: 

ls -lah
1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc
1   1   102
1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc
1   1   9
1294_S1_L008_R3_001.fastq.gz | head -2 | tail -1 | wc
1   1   9
1294_S1_L008_R4_001.fastq.gz | head -2 | tail -1 | wc
1   1   102

R1 and R4 contain the biological reads and R2 and R3 contain the indices. I knew this because the read lengths for the indices will be smaller as they contain the barcodes not the sequences. 

I looked through the quality scores to find the phred 33, this is because in phred 64, the bulk of the phred scores will be lowercase. 

2. For the histograms, I wrote a bash script to run the bioscript.py python script for taking the mean average of each base position's quality score. Initially I was using the code from PS9 which created a numpy array. This worked for R2 and R3 but even adding multiple nodes and creating a 8bit array was not enough to lower the amount of memory used on talapas. Because of this, I switched to PS4's list method instead. It took a bit to do but eventually worked. 

For the quality score cutoff, I looked up Illumina's manual to double check that quality score 30 was a good benchmark for cutting off scores. Several papers I've read used 30 as the cutoff and it made sense to me as most reads were not below 30 anyway. 

I did the following code to assess the number of records containing N's in the barcodes: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep -c 'N' 

## Part 2 – Develop an algorithm to de-multiplex the samples

For part 2, Jason went over the steps to create several folders full of matched, unknown, hopped, and low-quality reads. I included the low-quality folders because it felt wrong to put low-quality but matched files into an unknown folder. 

I created a record for each basic scenario that would put the records into each category of folder. Once this was all done, I pushed the files to github. 

## ---- Assignment the second ----

For the second assignment, I peer-reviewed two of my classmates psuedocode. 

## Peter's:

Line 39: I like what you have here, maybe include a function for grabbing the correct lines as well?

Line 88: What exactly are you concatenating and where are you putting it?

Line 92: What is your cutoff for the quality score?

Line 97: I like the way this is organized, really nicely laid out. Outlining the code like this seems painful but once you have it done I bet it will speed up the coding process!

## Logan's:

Line 20: Interesting... how/why did you come up with two different quality score cutoffs?

Line 46: The functions look really good, well-organized and laid out.

## Isis's (doesn't have line numbers):

What's your quality score cutoff? :D

       "If Index1 or Index2 contain any "N"

        Move the Read1 (record) to "Unknown_R1.fq" and move Read2 (record) to "Unknown_R2.fq".while while            redirecting change the header of Read1 and Read2 record to include Index1:Index2*"

Make sure the N you're talking about is in the sequence not anywhere else in the record.

Your pseudocode looks good. You'll be all ready for coding now!

## ---- Assignment the third ----

For assignment the third, we actually had to write out the code for demultiplexing the samples.

The first step for me was importing all the packages I needed for this undertaking. I knew I needed matplotlib for the graph at the end, I needed gzip for unzipping the files, bioinfo to call my functions for phred and the new functions I would write, argparse for assigning variables, etc. Later, I added itertools with the permutations package. 

Then I could actually get to coding. The first thing was to set my variables using argparse. I knew I'd need some for the reads, indices, and the known barcodes list as well as for the output plot. 

Then I needed to open the files using gzip and write files that I would put the results into: 

r1 = gzip.open(args.r1,"rt")
r2 = gzip.open(args.r2, "rt")
i1 = gzip.open(args.i1, "rt")
i2 = gzip.open(args.i2, "rt")
k = open(args.k, "rt")

r1_unknown = open("r1_unknown.fq", "w")
r2_unknown = open("r2_unknown.fq", "w")
r1_hopped = open("r1_hopped.fq", "w")
r2_hopped = open("r2_hopped.fq", "w")

--------------------------------------------------

Next, I set about parsing the barcodes file. I wrote a function that split the tsv into columns and pulled out the barcodes column and put those strings into a set to be called later. 

Next came the dictionary making. I wrote a for loop for writing the file names based on barcodes with the help of Rachel. 

Next, I created two more dictionaries for the hopped and matched barcodes, using the permutations package in itertools to match all the hopped barcodes with different ones. Now I had a place for everything to go.

--------------------------------------------------

Next I set a counter for the unknown reads to 0 so that I could get a value for every time an unknown read would be found. 

Just beneath that, I added a while True statement and used readline for all the header lines for each file (for adding the barcodes to) and made an array for each file so that I could reconstruct the files into FASTQ format later on. 

I wrote two functions that I put into bioinfo: a function that created a complement to a nucleotide and a function that created a reverse complement (using that complement function).

With Leslie's help, I set that reverse comp to a variable to save computing power. 

Next, I created two string concatenations for the r1 and the r2 strins that basically reconstructs the FASTQ files with the altered header lines. 

Now I could start writing the loops.

--------------------------------------------------

The first thing I did was write a for loop for identifying the quality score cutoff. If there was a base read in the indexes that was below 30 I would throw the read into the unknown files and add to the unknown counter. If it passed the quality check, it would move on. 

I first checked if the barcode in index one was in the barcodes set I made earlier. If not, it would go in unknown. Then I checked if the reverse complement of the barcodes in the second index are in index one. If the reverse comp was equal to the index, it would go in matched. If not, it would move on. If the reverse comp of index 2 is in the known barcodes, it would go in hopped, if not, unknown.

To run it, I used an sbatch script using the files taken from the shared folder on talapas. 

Basically, that was it. I wrote some simple for loops for the dictionary output and summing the values of the dictionaries to get the results. My graph didn't turn out because I forgot to write the output file so I made the graph again using the values I had outputted. Here's what that looked like:

#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

labels = 'Matched', 'Hopped', 'Unknown'
sizes = [79.40787630204026, 0.13918225582949836, 20.45294144213024]
explode = (0.05,0.05,0.05)
fig1, ax1 = plt.subplots()
ax1.pie(sizes, shadow=True, labels=labels, autopct='%1.1f%%', colors=('#2AAEC6', '#1CC9EB', '#33E0EA'), startangle=90, explode = explode)
ax1.axis('equal')
plt.tight_layout()
plt.subplots_adjust(top=0.88)
plt.title('Proportion of Matched, Hopped, and Unknown Reads')
plt.savefig("plot.png")

The final results matched up with what everyone else was finding. I had the same amount of total reads and about 80% matched with 0.1% hopped. 

Total Matched
288446518
Total Hopped
505575
Total Unknown
74294642
Total Records
363246735
Percentage of Matched
79.40787630204026
Percentage of Hopped
0.13918225582949836
Percentage of Unknown
20.45294144213024

If I had to change anything, I think I would add more elif statements so I didn't have to nest the for loops as much which probably slowed my program down a bit. Mine took about 2.5 hrs to run and SJ's took about 1 hr (they used the elif more).