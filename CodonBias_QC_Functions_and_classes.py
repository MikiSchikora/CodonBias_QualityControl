# this is a module that contains all the functions and classes used for CodonBias_QC

import numpy as np
import re
from Bio import pairwise2

def calc_tAI(sequence, type='tAI'):

    """ This function takes a  DNA sequence object or sequence that comes from a multifasta and outputs a tAI floats.
    This is for Saccharomyces cerevisiae."""

    if sequence.__class__.__name__=="SeqRecord":
        sequence = sequence.seq._data

    # define yeast tAI values taken from Table S2 in Tuller et al.
    tAI_weights = {
        'TTT': 0.2703,
        'TTC': 0.6158,
        'TTA': 0.4310,
        'TTG': 0.7537,
        'TCT': 0.6773,
        'TCC': 0.4877,
        'TCA': 0.1848,
        'TCG': 0.1207,
        'TAT': 0.2163,
        'TAC': 0.4926,
        'TAA': None,
        'TAG': None,
        'TGT': 0.1081,
        'TGC': 0.2463,
        'TGA': None,
        'TGG': 0.3695,
        'CTT': 0.0270,
        'CTC': 0.0616,
        'CTA': 0.1847,
        'CTG': 0.0591,
        'CCT': 0.1232,
        'CCC': 0.0887,
        'CCA': 0.6158,
        'CCG': 0.1970,
        'CAT': 0.1892,
        'CAC': 0.4310,
        'CAA': 0.5542,
        'CAG': 0.2389,
        'CGT': 0.3695,
        'CGC': 0.2660,
        'CGA': 0.000037,
        'CGG': 0.0616,
        'ATT': 0.8005,
        'ATC': 0.5764,
        'ATA': 0.1232,
        'ATG': 0.6158,
        'ACT': 0.6773,
        'ACC': 0.4877,
        'ACA': 0.2464,
        'ACG': 0.1404,
        'AAT': 0.2703,
        'AAC': 0.6158,
        'AAA': 0.4310,
        'AAG': 1.0000,
        'AGT': 0.0541,
        'AGC': 0.1232,
        'AGA': 0.6773,
        'AGG': 0.2783,
        'GTT': 0.8621,
        'GTC': 0.6207,
        'GTA': 0.1232,
        'GTG': 0.1626,
        'GCT': 0.6773,
        'GCC': 0.4877,
        'GCA': 0.3080,
        'GCG': 0.0985,
        'GAT': 0.4325,
        'GAC': 0.9852,
        'GAA': 0.8621,
        'GAG': 0.3990,
        'GGT': 0.4325,
        'GGC': 0.9852,
        'GGA': 0.1847,
        'GGG': 0.1823
    }

    # also nTE wights:
    nTE_weights = {
        'TAA': None,
        'TAG': None,
        'TGA': None,
        'TTT': 0.096,
        'TTC': 0.285,
        'TTA': 0.147,
        'TTG': 0.227,
        'TCT': 0.248,
        'TCC': 0.292,
        'TCA': 0.095,
        'TCG': 0.133,
        'TAT': 0.106,
        'TAC': 0.284,
        'TGT': 0.122,
        'TGC': 0.513,
        'TGG': 0.312,
        'CTT': 0.021,
        'CTC': 0.112,
        'CTA': 0.125,
        'CTG': 0.052,
        'CCT': 0.083,
        'CCC': 0.125,
        'CCA': 0.285,
        'CCG': 0.37,
        'CAT': 0.127,
        'CAC': 0.483,
        'CAA': 0.176,
        'CAG': 0.181,
        'CGT': 0.473,
        'CGC': 1.00,
        'CGA': 0.00013,
        'CGG': 0.363,
        'ATT': 0.23,
        'ATC': 0.28,
        'ATA': 0.069,
        'ATG': 0.261,
        'ACT': 0.281,
        'ACC': 0.318,
        'ACA': 0.131,
        'ACG': 0.166,
        'AAT': 0.072,
        'AAC': 0.216,
        'AAA': 0.094,
        'AAG': 0.263,
        'AGT': 0.035,
        'AGC': 0.119,
        'AGA': 0.265,
        'AGG': 0.281,
        'GTT': 0.318,
        'GTC': 0.421,
        'GTA': 0.097,
        'GTG': 0.134,
        'GCT': 0.247,
        'GCC': 0.31,
        'GCA': 0.171,
        'GCG': 0.145,
        'GAT': 0.101,
        'GAC': 0.41,
        'GAA': 0.159,
        'GAG': 0.185,
        'GGT': 0.139,
        'GGC': 0.889,
        'GGA': 0.159,
        'GGG': 0.275,
    }

    #codons and codon bias calculation
    codons = list(map(''.join, zip(*[iter(sequence)]*3)))

    if len(codons)>=10:
        if type=='tAI':
            cb = geo_mean([tAI_weights[c] for c in codons if tAI_weights[c]])

        elif type=='nTE':
            cb = geo_mean([nTE_weights[c] for c in codons if nTE_weights[c]])
    else:
        cb = 0

    return cb


def get_first_orf(chr_seq,tSS):

    """ This function takes a nucleotide sequence (chr_seq) and a tSS, which is the index in which a random transcript starts.
    It returns the first ORF possible."""

    sliced_chr = chr_seq[tSS:] # chromosomal region where to search

    pattern = re.compile(r'(?=(ATG(?:...)*?'r')(?=TAG|TGA|TAA))') #orf pattern

    try:
        first_orf = pattern.findall(sliced_chr)[0][:-3]
    except (AttributeError, IndexError):
        first_orf = ''

    return first_orf


def geo_mean(iterable):

    """ This calculates the geometric mean of an iterable object"""
    a = np.array(iterable)
    return a.prod()**(1.0/len(a))

def find_orf(seq):

    """ This takes a Seq object and finds the ORF that starts at position 1"""

    stop_pos = seq.translate()._data.find('*')*3
    return seq[0:stop_pos+3]


def calc_id(seqA, seqB):

    """ Takes two sequences and returns it's %identity"""

    global_align = pairwise2.align.globalxx(seqA._data, seqB._data)
    seq_length = max(len(seqA), len(seqB))
    matches = global_align[0][2]

    return (matches / seq_length) * 100