import sys
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez
import seaborn as sns

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

def exact_position(locus, position):
    """
    Calculates sum of the chromosomal base number position the gene starts at, position in the transcript the read was taken, and the genomic position of that chromosome
    Args:
        locus (int): The base the gene starts at on a particular chromosome
        position (int): The positon of the read within the transcript
        genomic_position (int): The position of that chromosome in the Saccharomyces cerevisiae genome
    Returns:
       exact_position (int): The exact base in the genome which the read was taken at
    """
    exact_position = locus + position
    return exact_position
def calculate_intensity_difference(c1_mean_intensity, c2_mean_intensity):
    difference = abs(c1_mean_intensity - c2_mean_intensity)
    return difference

input_file_path = "../../../../../data2/tom/RNA/results/outnanocompore_shift_stats.tsv"#"outnanocompore_results_short.tsv"#"../../../../../data2/tom/RNA/results/outnanocompore_shift_stats.tsv"
output_file_path = "../../../../data2/tom/RNA/signal_chr.tsv"#"test.tsv"# "../../../../data2/tom/RNA/signal_chr.tsv"
output_graph_path = "../../../../data2/tom/RNA/signal_chr_manhattan_plot_graph.png"#"test.png"#"../../../../data2/tom/RNA/signal_chr_manhattan_plot_graph.png"
#output_graph_path = "nanocompore_signal_genomic_location_graph.png"
df_original = pd.read_csv(input_file_path, sep='\t', header=0)
#df2 = df_original[df_original['GMM_logit_pvalue'].astype(str) != "nan"]
df = df_original.loc[df_original['ref_id'].str.contains("mrna")]
#df = df3.loc[(df3['GMM_logit_pvalue'].astype(float) <= 0.005)]
#Extract of transcript and chromosome information from the existing ref_id column
df = df.loc[df['ref_id'].str.contains("mrna")]
df["transcript_reference"] = df.apply(lambda row : extract_transcript_reference(row["ref_id"]), axis = 1)
df = df.loc[(df['transcript_reference'] != "Failed")]
df['chromosome_reference'] = df.apply(lambda row : extract_chrosomosome_reference(row['ref_id']), axis = 1)

#Change the each Saccharomyces cerevisiae chromosome RefSeq id into its actual chromosome number and then apply this to the datatframe
chromosome_refseq_dictionary = {"NC_001133.9":1,"NC_001134.8":2,"NC_001135.5":3,"NC_001136.10":4,"NC_001137.3":5,"NC_001138.5":6,"NC_001139.9":7,"NC_001140.6":8,"NC_001141.2":9,"NC_001142.9":10,"NC_001143.9":11,"NC_001144.5":12,"NC_001145.3":13,"NC_001146.8":14,"NC_001147.6":15,"NC_001148.4":16}
df['chromosome'] = df.apply(lambda row : match(row['chromosome_reference'], chromosome_refseq_dictionary), axis = 1)
chromosome_size_dictionary = {1:230218, 2:813184, 3:316620, 4:1531933, 5:576874, 6:270161, 7:1090940, 8:562643, 9:439888, 10:745751, 11:666816, 12:1078177, 13:924431, 14:784333, 15:1091291, 16:948066}
#df = df.dropna(how='all',axis=0, inplace=True) 
df = df[df['transcript_reference'].notna()]
#Use the RefSeq gene id to search for the gene id and then use that to get the chromosomal locus for every unique transcript_id, before matching these values this to the data frame. Use this locus, the transcript position, in addition to the cumulative genomic position of each chromosome to get the exact genomic location of every transcript position.
position_dictionary = {}
print(df.head())
for item in df['transcript_reference'].unique():
    gene_id = get_gene_id(item)
    position_dictionary[item] = get_locus(gene_id)
df['locus'] = df.apply(lambda row : match(row['transcript_reference'], position_dictionary), axis = 1)
df['exact_position'] = df.apply(lambda row : exact_position(row['locus'], row['pos']), axis = 1)
df['inter_condition_mean_intensity_difference'] = df.apply(lambda row : calculate_intensity_difference(row['c1_mean_intensity'], row['c2_mean_intensity']), axis = 1)
df = df.loc[(df['inter_condition_mean_intensity_difference'].astype(float) >= 2)]
df.to_csv(output_file_path, sep='\t', encoding='utf-8')
df.chromosome = df.chromosome.astype('category')
#df.chromosome = df.chromosome.cat.set_categories(['ch-%i' % i for i in range(16)], ordered=True) #This row apepars to wipe the chromosome stuff
df = df.sort_values('chromosome')
df['ind'] = range(len(df))
df_grouped = df.groupby(('chromosome'))
df.to_csv(output_file_path, sep='\t', encoding='utf-8')

# How to plot gene vs. -log10(pvalue) and colour it by chromosome?
print(df_grouped.head())
fig = plt.figure()
ax = fig.add_subplot(111)
colors = ['red','green','blue', 'yellow']
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='inter_condition_mean_intensity_difference', color=colors[num % len(colors)], ax=ax) # y = 'GMM_logit_pvalue', color=colors[num % len(colors)], ax=ax) 
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df)])
#ax.set_ylim([0, 3.5])
ax.set_xlabel('Chromosome')
ax.set_ylabel('Difference in mean signal intensity')
ax.set_title("Difference in mean signal intensity throughout the genome between conditions")
fig.savefig(output_graph_path)
#g = sns.FacetGrid(df, col="chr") #, height=3.5, aspect=.65)
#g.map_dataframe(sns.lmplot(x='exact_position', y='inter_condition_mean_intensity_difference', data=df, fit_reg=False)
#g.set_axis_labels("Location in Chromosome (kB)", "Inter-conditional difference in mean signal intensity")
#g.set(title = "Mean difference in signal intensity throughout the genome between conditions")
#g.savefig(output_graph_path)
print("done")
