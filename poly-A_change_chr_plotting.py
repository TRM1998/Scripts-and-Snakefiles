import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Entrez
from scipy.stats import mannwhitneyu
pd.options.mode.chained_assignment = None  # default='warn'

def extract_transcript_reference(reference_id):
    if reference_id != None:
        ref_id = reference_id.split("_")
        if len(ref_id) == 6:
            transcript_reference = ref_id[3] +  "_" + ref_id[4]
        else:
            transcript_reference = "Failed"
        return transcript_reference

def extract_chrosomosome_reference(ref_id):
    """
    Fetches the chromosome RefSeq id the ref_id string
    Args:
        ref_id (string): Contains the RefSeq chromosome and transcript ids
    Returns:
        chromosome_reference (string): Chromosome id
    """
    ref_id = ref_id.split("|")
    ref_id = ref_id[1].split("_")
    chromosome_reference = ref_id[0] + "_" + ref_id[1]
    return chromosome_reference



def match(key, dictionary):
    """
    Returns the value in a supplied dictionary for a particular key
    Args:
        key (string / int): Key
        dictionary (dict): Dictionary
    Returns:
       matched (string / int): The matched value for the key from the dictionary
    """
    matched = dictionary.get(key)
    return matched

def assign_type(type):
    return type

input_stats_file_path = "../../data2/tom/RNA/poly-A-MWU-filtered.tsv"
input_WT_path = "../../data2/tom/RNA/poly-A-WT.tsv"
input_SSS_path = "../../data2/tom/RNA/poly-A-SSS.tsv"
input_position_path = "../../data2/tom/RNA/outnanocompore_percentages.tsv"
output_file_path = "../../data2/tom/RNA/poly-A-chr.tsv"
output_graph_path = "../../data2/tom/RNA/poly-A-chr_change.png"

#df_stats = pd.read_csv(input_stats_file_path, sep='\t', header=0)
#df_stats = df_stats.iloc[700:744]
df_position = pd.read_csv(input_position_path, sep = "\t", header = 0)
#df_WT = pd.read_csv(input_WT_path, sep='\t', header=0)
#df_WT["type"] = df_WT.apply(lambda row : assign_type("WT"), axis = 1)
#df_SSS = pd.read_csv(input_SSS_path, sep='\t', header=0)
#df_SSS["type"] = df_SSS.apply(lambda row : assign_type("SSS"), axis = 1)

#df = pd.DataFrame(columns=('transcript_reference', "type", 'polya_length'))
#for transcript in df_stats['transcript_reference']:
#    df_1 = df_WT.loc[df_WT['transcript_reference'] == transcript]
#    df_2 = df_SSS.loc[df_SSS['transcript_reference'] == transcript]
#    df_merged = pd.concat([df_1, df_2])
#    df_melted = df_merged.melt(id_vars = ["transcript_reference", "type", "polya_length"]) #, value_vars = ['polya_length'])
    #Need to merge df_1 and df_2 and specify which dataset each one came from
#    df = df.append(df_melted)

input_tails_file_path = "../../data2/tom/RNA/significant_tail_differences.tsv"
df = pd.read_csv(input_tails_file_path, sep='\t', header=0)
#df["transcript"] = df.apply(lambda row : extract_transcript_reference(row["transcript_reference"]), axis = 1)
#df = df.loc[(df['transcript'] != "Failed")]
df['chromosome_reference'] = df.apply(lambda row : extract_chrosomosome_reference(row['transcript_reference']), axis = 1)
#chr_dictionary = {}
#position_dictionary = {}
#locus_dictionary = {}
#for transcript in df["transcript_id"].unique():
#    data = df_position.loc[df_position['transcript_reference'] == str(transcript)]
#    position = data['exact_genomic_position'].values[0:]
#    chromosome = data['chr'].values[0:]
#    locus = data['locus'].values[0:]
#    print(locus)
#    chr_dictionary[transcript] = chromosome
    #position_dictionary[transcript] = position
#    locus_dictionary[transcript] = locus
def get_gene_id(ref_id):
    """
    Fetches the gene id of a matching RefSeq id from Entrez
    Args:
        ref_id (string): RefSeq id
    Returns:
       gene_id (string): Gene id 
    """
    Entrez.email = "zx20935@bristol.ac.uk"
    handle = Entrez.efetch(db="nucleotide", id=ref_id, rettype="gb", retmode="text")
    read = handle.read()
    matched_line = ([line for line in read.split('\n') if "GeneID" in line])[0]
    gene_id = (matched_line.split(":"))[1]
    gene_id = gene_id[:-1]
    return gene_id

def get_locus(ref_id):
    """
    Fetches the chromosomal locus of a matching gene id from Entrez
    Args:
        ref_id (string): Gene id
    Returns:
       locus (int): Chromosomal locus / The base the gene starts at on a particular chromosome
    """
    Entrez.email = "zx20935@bristol.ac.uk"
    handle = Entrez.efetch(db="gene", id=ref_id, rettype="gb", retmode="text")
    read = handle.read()
    lines = read.splitlines()
    line = lines[6].split("(")
    location = line[1].split(".", 1)
    locus = int(location[0])
    return locus
chromosome_refseq_dictionary = {"NC_001133.9":1,"NC_001134.8":2,"NC_001135.5":3,"NC_001136.10":4,"NC_001137.3":5,"NC_001138.5":6,"NC_001139.9":7,"NC_001140.6":8,"NC_001141.2":9,"NC_001142.9":10,"NC_001143.9":11,"NC_001144.5":12,"NC_001145.3":13,"NC_001146.8":14,"NC_001147.6":15,"NC_001148.4":16}
df['chromosome'] = df.apply(lambda row : match(row['chromosome_reference'], chromosome_refseq_dictionary), axis = 1)
df = df[df['transcript_id'].notna()]
position_dictionary = {}
for item in df['transcript_id'].unique():
    gene_id = get_gene_id(item)
    position_dictionary[item] = get_locus(gene_id)
df['locus'] = df.apply(lambda row : match(row['transcript_id'], position_dictionary), axis = 1)

def get_change(old, new):
    new = float(new)
    old = float(old)
    if old == new:
        return 100.0
    try:
        #Want to get the difference change -> to - will be a decrease from WT to SSS
        return ((new - old) / old) * 100.0 #(abs(new - old) / old) * 100.0 
    except ZeroDivisionError:
        return 0

df['tail_change'] = df.apply(lambda row : get_change(row["average_tail_WT_length"], row["average_tail_SSS_length"]), axis = 1)
#df['chromosome'] = df.apply(lambda row : match(row['transcript_id'], chr_dictionary), axis = 1)
#df['exact_position'] = df.apply(lambda row : match(row['transcript_id'], position_dictionary), axis = 1)
#df['locus'] = df.apply(lambda row : match(row['transcript_id'], locus_dictionary), axis = 1)
#I may end up generating a face_plot of chromosomes + loci

df.chromosome = df.chromosome.astype('category')
#df.chromosome = df.chromosome.cat.set_categories(['ch-%i' % i for i in range(16)], ordered=True) #This row apepars to wipe the chromosome stuff
df = df.sort_values('chromosome')
df['ind'] = range(len(df))
df_grouped = df.groupby(('chromosome'))
df.to_csv(output_file_path, sep='\t', encoding='utf-8')

# How to plot gene vs. -log10(pvalue) and colour it by chromosome?t
print(df_grouped.head())
fig = plt.figure()
ax = fig.add_subplot(111)
colors = ['red','green','blue', 'yellow']
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='tail_change', color=colors[num % len(colors)], ax=ax) # y = 'GMM_logit_pvalue', color=colors[num % len(colors)], ax=ax) 
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df)])
#ax.set_ylim([0, 3.5])
ax.set_xlabel('Chromosome')
ax.set_ylabel('Relative change in poly-A tail length')
ax.set_title("Relative change in poly-A tail length against genomic position")
fig.savefig(output_graph_path)
print("done")
