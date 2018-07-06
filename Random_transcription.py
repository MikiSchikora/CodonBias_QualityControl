# This script asks things about what happens if you have random transcription of the genome

import os
import sys
from Bio import SeqIO
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# add modules
sys.path.append('/Volumes/GoogleDrive/Mi unidad/CareyLab/Miki/Python_modules/')
import CodonBias_QC_Functions_and_classes as CB
from GFF3parser import *

# define the plot to do, from command-line
type_plot = sys.argv[1]
if type_plot not in ['length_ORFs','tAI_in_random','nTE_in_random','tAI_in_random_weight']:
    type_plot = 'tAI_in_random_weight'


########## compare tAI and CAI of native vs random genes ###############

orf_coding_path = '/Users/miki_schikora_tamarit/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/Yeast/S288C_reference_genome_R64-2-1_20150113/orf_coding_all_R64-2-1_20150113.fasta'
transcripts_path = '/Users/miki_schikora_tamarit/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/Yeast/S288C_reference_genome_R64-2-1_20150113/orf_trans_all_R64-2-1_20150113.fasta'
#RNAexpression_path = '~/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/vanDijk15/FitFlowRNAseq.tab'
#RNAexpression_path = '~/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/vanDijk15/FitFlow__Fast_Slow.tab'
RNAexpression_path = '~/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/Carey15/Carey15_RNAseq.tab'
genome_path = '/Users/miki_schikora_tamarit/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/Yeast/S288C_reference_genome_R64-2-1_20150113/genome.fasta'

# generate the tAI for random ORFs

# generate a set of orfs comming from random positions in the genome, through all chromosomes
n_orfs = 100 # number of orfs to generate of each chromosome
tAI_random  = []
nTE_random = []

for chr in SeqIO.parse(genome_path, "fasta"):

    print('generating random ORFs in chr %s'%(chr.id))

    chr_seq = chr.seq._data

    # generate random orfs n_orfs times:
    for n in range(0,n_orfs):
        tSS  = random.randint(0,len(chr_seq))# simulate a transcription start site. It's an index, so it starts with 0

        tAI_random.append(CB.calc_tAI(CB.get_first_orf(chr_seq,tSS)))# generate the random orf
        nTE_random.append(CB.calc_tAI(CB.get_first_orf(chr_seq, tSS), type='nTE'))

tAI_random = np.array(tAI_random)
tAI_random = tAI_random[tAI_random>0]

nTE_random = np.array(nTE_random)
nTE_random = nTE_random[nTE_random>0]

# record the distribution of ORF length
length_genome = np.array([len(o.seq._data)/3 for o in SeqIO.parse(orf_coding_path, "fasta")])

# assign an expression weight to each gene
#Expression = pd.read_table(RNAexpression_path).set_index('TranscriptID')['FastFPKM'].to_dict()
#Expression = pd.read_table(RNAexpression_path,names=['ORF','name','position','FastFPKM','SlowFPKM']).set_index('ORF')['FastFPKM'].to_dict()
Expression = pd.read_table(RNAexpression_path).groupby('target_id').sum().to_dict()['tpm']


# calc tAI for native ORFs
tAI_genome = np.array(list(map(CB.calc_tAI, SeqIO.parse(orf_coding_path, "fasta"))))
tAI_genome = tAI_genome[tAI_genome>0]

# calc weighted tAI
tAI_genome_weight = []
n_orfs_added=0
for orf in SeqIO.parse(orf_coding_path, "fasta"):

    id = orf.id
    try:
        expression_w = Expression[id]
        tAI_genome_weight += [CB.calc_tAI(orf)]*int(expression_w*1)
        n_orfs_added += 1
    except KeyError:
        pass
tAI_genome_weight = np.array(tAI_genome_weight)
tAI_genome_weight = tAI_genome_weight[tAI_genome_weight>0]


nTE_genome = np.array([CB.calc_tAI(s,type='nTE') for s in SeqIO.parse(orf_coding_path, "fasta")])
nTE_genome = nTE_genome[nTE_genome>0]


# plot genome length

if type_plot == 'length_ORFs':
    plt.close()
    sns.set(style="white", palette="muted", color_codes=True)
    sns.distplot(length_genome, color="r", bins=400, hist=True)
    plt.tight_layout()
    plt.xlabel('nuber of codons')
    plt.ylabel('gene occurence')
    plt.show()

if type_plot == 'tAI_in_random':
    #plot of tAI
    plt.close()
    sns.set(style="white", palette="muted", color_codes=True)
    sns.distplot(tAI_genome, color="m", label='native genes')
    sns.distplot(tAI_random, color="b", label='random ORFs')
    plt.tight_layout()
    plt.xlabel('tRNA adaptation index (tAI)')
    plt.ylabel('frequency density')
    plt.legend()
    plt.show()


if type_plot == 'tAI_in_random_weight':
    #plot of tAI
    plt.close()
    sns.set(style="white", palette="muted", color_codes=True)
    sns.distplot(tAI_genome_weight, color="m", label='native genes')
    sns.distplot(tAI_random, color="b", label='random ORFs')
    #plt.tight_layout()
    plt.xlabel('tRNA adaptation index (tAI)')
    plt.ylabel('frequency density')
    plt.legend()
    plt.show(format='pdf')
    plt.savefig("graph.pdf",format='pdf')

if type_plot == 'nTE_in_random':
    #plot of nTE
    plt.close()
    sns.set(style="white", palette="muted", color_codes=True)
    sns.distplot(nTE_genome, color="m", label='native genes')
    sns.distplot(nTE_random, color="b", label='random ORFs')
    plt.tight_layout()
    plt.xlabel('nTE')
    plt.ylabel('frequency density')
    plt.legend()
    plt.show()



########################################################################

