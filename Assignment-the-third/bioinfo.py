# Author: Jack Peplinski jpeplin5@uoregon.edu

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.5"
             
import re
import math
import matplotlib.pylab as plt
import numpy as np
import argparse             

DNA_bases = "ACTGN"
RNA_bases = "ACUGN"



def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return (ord(letter) - 33)

def qual_score(phred_score: str) -> float:
    """Calculates the average quality score from a string of phred scores"""
    length = len(phred_score)
    x = 0
    sum = 0
    for item in phred_score:
        bases = convert_phred(item)
        x += 1
        sum = sum + bases
    return(sum/length)

DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')

def validate_base_seq(seq: str,RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("U" if RNAflag else "T") + seq.count("G") + seq.count("C")

def gc_content(seq):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq), "String contains invalid characters"
    seq = seq.upper()
    Gs = seq.count("G")
    Cs = seq.count("C")
    return (Gs+Cs)/len(seq)

def oneline_fasta():
    first = True
    with open('file', 'r') as f, open("newfile","w") as outfh:
        for line in f:
            line = line.strip("\n")
            if line[0] == ">":
                if not first:
                    outfh.write(header+'\n')
                    outfh.write(seq+'\n')
                first = False
                header = line
                seq = ""
            else:
                seq += line
        outfh.write(header+'\n')
        outfh.write(seq+'\n')

comp_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}

#making the complements
def comp(seq):
    bases = list(seq)
    bases = [comp_dict[base] for base in bases]
    return ''.join(bases)

#making the reverse complements
def rev_comp(s):
    return comp(s[::-1])
    
if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests") 

