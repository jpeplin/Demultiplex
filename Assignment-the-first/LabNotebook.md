Lab Notebook
Jack Peplinski
7/26/22

This is my lab notebook for the demultiplexing assignments. 

## ---- Assignment the first ----

## Part 1 – Quality Score Distribution per-nucleotide

1.
Initial Data Exploration: 

$ ls -lah
$ 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc
1   1   102
$ 1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc
1   1   9
$ 1294_S1_L008_R3_001.fastq.gz | head -2 | tail -1 | wc
1   1   9
$ 1294_S1_L008_R4_001.fastq.gz | head -2 | tail -1 | wc
1   1   102

R1 and R4 contain the biological reads and R2 and R3 contain the indices. I knew this because the read lengths for the indices will be smaller as they contain the barcodes not the sequences. 

I looked through the quality scores to find the phred 33, this is because in phred 64, the bulk of the phred scores will be lowercase. 

2. For the histograms, I wrote a bash script to run the bioscript.py python script for taking the mean average of each base position's quality score. Initially I was using the code from PS9 which created a numpy array. This worked for R2 and R3 but even adding multiple nodes and creating a 8bit array was not enough to lower the amount of memory used on talapas. Because of this, I switched to PS4's list method instead. It took a bit to do but eventually worked. 

For the quality score cutoff, I looked up Illumina's manual to double check that quality score 30 was a good benchmark for cutting off scores. Several papers I've read used 30 as the cutoff and it made sense to me as most reads were not below 30 anyway. 

I did the following code to assess the number of records containing N's in the barcodes: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep -c 'N' 

## Part 2 – Develop an algorithm to de-multiplex the samples

For part 2, Jason went over the steps to create several folders full of matched, unknown, hopped, and low-quality reads. I included the low-quality folders because it felt wrong to put low-quality but matched files into an unknown folder. 

I created a record for each basic scenario that would put the records into each category of folder. Once this was all done, I pushed the files to github. 







