# This script is to ask about the relationship between tAI and transcript errors. Simulated

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from Bio.Alphabet import generic_dna

# add modules
sys.path.append('/Volumes/GoogleDrive/Mi unidad/CareyLab/Miki/Python_modules/')
import CodonBias_QC_Functions_and_classes as CB
import GFF3parser


# define paths
genome_path = '/Users/miki_schikora_tamarit/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/Yeast/S288C_reference_genome_R64-2-1_20150113/genome.fasta'
orf_coding_path = '/Users/miki_schikora_tamarit/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/Yeast/S288C_reference_genome_R64-2-1_20150113/orf_coding_all_R64-2-1_20150113.fasta'
RNAexpression_path = '~/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/Carey15/Carey15_RNAseq.tab'
utr3_path = '/Users/miki_schikora_tamarit/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/Nagalakshmi08/Nagalakshmi_2008_3UTRs_V64.gff3'

# assign an expression value to each gene
Expression = pd.read_table(RNAexpression_path).groupby('target_id').sum().to_dict()['tpm']
Expression_probabilities = (pd.read_table(RNAexpression_path).groupby('target_id').sum()['tpm']/sum(Expression.values())).to_dict()

# put the genome in a dictionary
genome = {}
for chromosome in SeqIO.parse(genome_path, "fasta"):
    genome[chromosome.id] = chromosome.seq._data

# put the ORFs in a dictionary
ORFs = {}
for orf in SeqIO.parse(orf_coding_path, "fasta"):
    ORFs[orf.name] = orf.seq


#make a dictionary that includes ORF and 3'UTR
transcripts = {}
for record in GFF3parser.parseGFF3(utr3_path):

    if record is None or record.seqid=="chrmt":
        continue
    chr = record.seqid
    gene_id = record.attributes['ID'][0:-5]

    # generate 3'UTR seq
    if record.strand=="+":
        utr3 = Seq(genome[chr][record.start-1:record.end],generic_dna)
    if record.strand=="-":
        utr3 = Seq(genome[chr][record.start-1:record.end], generic_dna).reverse_complement()

    try:
        transcripts[gene_id] = [ORFs[gene_id],utr3,Expression_probabilities[gene_id]]
    except KeyError:
        continue

# create a dataframe with all the information for a set of insertions. each gene gets as many insertions as TPM it has

#indels_ds = {'transcript_id':[], 'native_tAI':[], 'indel_tAI':[], 'native_after_indel_tAI':[], 'after_indel_tAI':[],
#             'pctg_identity_prot':[], 'pctg_identity_nuc':[], 'utr3_len_native':[], 'utr3_len_indel':[]
#             }

indels_ds = {'transcript_id':[], 'native_tAI':[], 'indel_tAI':[], 'native_after_indel_tAI':[], 'after_indel_tAI':[],
             'utr3_len_native':[], 'utr3_len_indel':[], 'indel_pos': [], 'orf_l':[]
             }

t=0
# loop through all genes
for trans_id, trans_info in transcripts.items():

    t+=1

    # generate native protein sequence, needed below
    orf = trans_info[0]
    #prot_seq = orf.translate()
    native_tAI = CB.calc_tAI(orf)
    utr3_len_native = len(trans_info[1])
    orf_len_native = len(orf)

    #print(t, trans_id)


    # generate an insertion for each tpm. The probability is multiplied by the simulated number of transcripts
    for I in range(0,int(trans_info[2]*100000)):

        print(t, trans_id, I)

        # add obvious things:
        indels_ds['transcript_id'].append(trans_id)
        indels_ds['native_tAI'].append(native_tAI)
        indels_ds['utr3_len_native'].append(utr3_len_native)
        indels_ds['orf_l'].append(orf_len_native)

        # generate the ORF with a deletion
        ins_pos = random.randint(0,len(orf))
        indel_orf = CB.find_orf(orf[:ins_pos] + orf[(ins_pos+1):] + trans_info[1])
        indels_ds['indel_tAI'].append(CB.calc_tAI(indel_orf))
        indels_ds['utr3_len_indel'].append(len(orf) + len(trans_info[1]) - len(indel_orf))
        indels_ds['indel_pos'].append(ins_pos)

        # align and add to ds
        #indel_prot_seq = indel_orf.translate()
        #indels_ds['pctg_identity_prot'].append(CB.calc_id(prot_seq,indel_prot_seq))
        #indels_ds['pctg_identity_nuc'].append(CB.calc_id(orf, indel_orf))

        # calculate the natice tAI from the insertion on (from the next codon on):
        while ins_pos%3 != 0: # add until starting position is multiple of 3
            ins_pos+=1
        indels_ds['native_after_indel_tAI'].append(CB.calc_tAI(orf[ins_pos:]))
        indels_ds['after_indel_tAI'].append(CB.calc_tAI(indel_orf[ins_pos:]))


#  convert into dataframe:
INDELS_df = pd.DataFrame(indels_ds)

# load TATA
TATA_info = pd.read_table("~/Google Drive File Stream/Mi unidad/CareyLab/ExternalData/Basehoar04/Basehoar04_TATA.csv", delimiter=";")



INDELS_df = pd.merge(TATA_info,INDELS_df, left_on="Gene", right_on="transcript_id")

#~/Google\ Drive\ File\ Stream/Mi\ unidad/CareyLab/ExternalData/Basehoar04/

idx = (INDELS_df.native_tAI>0) & (INDELS_df.indel_tAI>0) & (INDELS_df.index < 100000) & \
      (INDELS_df.after_indel_tAI>0) & (INDELS_df.native_after_indel_tAI>0) & \
      (INDELS_df.utr3_len_native>0) & (INDELS_df.utr3_len_indel>0) & (INDELS_df.orf_l>30) \
      & (INDELS_df['1=TATA; 0=TATA-less']!=np.nan)
ds = INDELS_df.loc[idx]

# plot
#X = ds['native_tAI']; Xlabel = 'native tAI'; Xlim = [0.25,0.7]
#X = ds['native_after_indel_tAI']; Xlabel = 'tAI after IN/DEL of native gene'; Xlim = [0.1,0.8]
#X = ds['native_after_indel_tAI'] - ds['after_indel_tAI']; Xlabel = "tAI reduced by IN/DEL"; Xlim = [0,0.5]
#X = ds['after_indel_tAI'] / ds['native_after_indel_tAI'] ; Xlabel = "tAI IN/DEL / native tAI"; Xlim = [0.2,0.8]
#X = ds['utr3_len_indel'] - ds['utr3_len_native'] ; Xlabel = "3'UTR increased by IN/DEL"; Xlim = [0,1000]
X = ds['utr3_len_native']; Xlabel = "native 3'UTR length"; Xlim = [0,900]

#Y = ds['indel_tAI']; Ylabel = 'with IN/DEL tAI'; Ylim = [0.25,0.7]
#Y = ds['after_indel_tAI']; Ylabel = 'tAI after IN/DEL'; Ylim = [0.1,0.8]
#Y = ds['utr3_len_indel'] - ds['utr3_len_native']; Ylabel = "3'UTR increased by IN/DEL"; Ylim = [0,650]
#Y = ds['utr3_len_indel'] / ds['utr3_len_native']; Ylabel = "3'UTR IN/DEL / native 3'UTR length"; Ylim = [0,10]

Y = ds['utr3_len_indel']; Ylabel = "3'UTR length with IN/DEL"; Ylim = [0,900]
#Y = (((ds.orf_l - ds.indel_pos)/ds.orf_l)*100).sort_values(); Ylabel = "% protein removed by IN/DEL"; Ylim = [0,100]

color = ds['utr3_len_indel'].sort_values(); color_label = "resulting 3'UTR length after IN/DEL"
#color = (((ds.orf_l - ds.indel_pos)/ds.orf_l)*100).sort_values(); color_label = "% protein removed by IN/DEL"

type="ratio"

# scatter
if type=="scatter":
    plt.close()
    sns.set(style="white", palette="colorblind", color_codes=True)
    plt.scatter(X, Y, alpha = 0.1, c = color, cmap = 'viridis')
    plt.plot(Xlim,Xlim,'k--')
    #plt.plot(X.sort_values(),pd.rolling_median(Y.sort_values(), window=100), 'k-')
    cbar = plt.colorbar(); cbar.set_label(color_label)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.xlim(Xlim)
    plt.ylim(Ylim)
    plt.title('How tAI varies with in/dels')
    plt.show()

if type=="density":
    #scatter density
    plt.close()
    sns.set(style="white")
    sns.kdeplot(X, Y, cmap="Reds", shade=True, shade_lowest=False)
    #g = sns.jointplot(X[(Y<800) & (X>0) & (Y>0)], Y[(Y<800) & (X>0) & (Y>0)], kind="kde", size=7, space=0)
    #g = sns.jointplot(X, Y, kind="kde", size=7, space=0)

    plt.plot(Xlim,Xlim,'k--')
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    #g.set_axis_labels(Xlabel,Ylabel)
    plt.xlim(Xlim)
    plt.ylim(Ylim)
    plt.title('Density plot')
    plt.show()

if type=="ratio":
    # ratio plots
    plt.close()
    plt.scatter(X, Y, alpha = 0.1, c = color, cmap = 'viridis')
    #plt.plot(Xlim,Xlim,'k--')
    #plt.plot(X.sort_values(),pd.rolling_median(Y.sort_values(), window=100), 'k-')
    cbar = plt.colorbar(); cbar.set_label(color_label)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.xlim(Xlim)
    plt.ylim(Ylim)
    plt.title('Looking at ratios')
    plt.show()

if type=="boxplot":

    plt.close()

    purple = (0.66, 0.486, 0.753)
    blue = (0.396, 0.529, 0.812)
    color_code_purple_blue = {"without FE": purple, "with FE":blue}

    sns.set_style("whitegrid")

    X_data = pd.Series(["without FE"]*len(ds) + ["with FE"]*len(ds))
    #Y_data = list(ds.native_after_indel_tAI.append(ds.after_indel_tAI))
    Y_data = list(np.log10(ds.utr3_len_native.append(ds.utr3_len_indel)))



    bp_ds = pd.DataFrame(data={"x_data":X_data, "y_data":Y_data, 'tata_type':ds['1=TATA; 0=TATA-less']})

    ax = sns.violinplot(x="x_data", y="y_data", data=bp_ds, palette=color_code_purple_blue)

    #ax = sns.violinplot(x="x_data", y="y_data", data=bp_ds, hue='tata_type', split=True)
    #hue
    #ax.set_yscale("log")

    plt.xlabel("")
    plt.ylabel("")
    #plt.ylim([-5,1000])
    plt.show()



