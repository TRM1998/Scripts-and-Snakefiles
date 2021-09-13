import sys
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez
import seaborn as sns

pd.options.mode.chained_assignment = None  # default='warn'

    
def match_series(key, series):
    """
    Returns the value in a supplied series for a particular key
    Args:
        key (string / int): Key
        Series: Series
    Returns:
       matched (string / int): The matched value for the key from the series
    """
    matched = series.get(key)
    return matched

def calculate_condition_difference(WT1, WT2, SSS1, SSS2):
    average_WT = (WT1 + WT2) /2
    average_SSS = (SSS1 + SSS2) /2
    difference = abs(average_WT - average_SSS)
    return difference

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

def  get_transcript_length(transcript_id):
    #transcript_id = "NM_001178202.2"
    Entrez.email = "zx20935@bristol.ac.uk"
    handle = Entrez.efetch(db="nucleotide", id=transcript_id, rettype="gb", retmode="text")
    read = handle.read()
    matched_line = ([line for line in read.split('\n') if "LOCUS" in line])[0]
    transcript_length = (matched_line.split())[2]
    return transcript_length

def calculate_relative_position(pos, length):
    relative_position = (int(pos) / int(length)) * 100
    return relative_position

input_file_path = "../../data2/tom/RNA/all_RNA_mod_rate_changes.tsv"
#output_file_path = "../../../../data2/tom/RNA/outnanocompore_density_transcript_position.tsv"
output_graph_path = "../../../../data2/tom/RNA/modification_density_transcript_position_graph.png"
df_original = pd.read_csv(input_file_path, sep='\t', header=0)
df2 = df_original[df_original['GMM_logit_pvalue'].astype(str) != "nan"]
df3 = df2.loc[df2['ref_id'].str.contains("mrna")]
df = df3.loc[(df3['GMM_logit_pvalue'].astype(float) <= 0.01)]
df_filtered = df #df.loc[(df['relative_average_modification_rate'].astype(float) >= 20)]
transcript_lengths = {}
for item in df_filtered['transcript_reference'].unique():
    transcript_lengths[item] = get_transcript_length(item)
transcript_modification_counts = df_filtered['transcript_reference'].value_counts()
#df_filtered['WT_SSS_difference'] = df_filtered.apply(lambda row : calculate_condition_difference(row['WT1'], row['WT2'], row['SSS1'], row['SSS2']), axis = 1)
df_filtered['number_of_differential_modifications'] = df_filtered.apply(lambda row : match_series(row['transcript_reference'], transcript_modification_counts), axis = 1)
df_filtered['transcript_length'] = df_filtered.apply(lambda row : match(row['transcript_reference'], transcript_lengths), axis = 1)
df_filtered['relative_position_in_transcript'] = df_filtered.apply(lambda row : calculate_relative_position(row["pos"], row['transcript_length']), axis = 1)

#Introduce a cut-off for % difference
#df_filtered = df_filtered.loc[(df_filtered['relative_average_modification_rate'].astype(float) <= 1000)]

plot = sns.kdeplot(df_filtered['relative_position_in_transcript'])
#fig = plt.gcf()
#fig.set_size_inches( 16, 10) #16, 10, 40, 
plot.set(xlabel = 'Relative position in transcript (%)', ylabel = 'Density of modification sites')
plot.set_xlim([0, 100])
plot.set(title = "Relative position in transcript against differential modification site density")
plot = plot.get_figure()
#plot.legend.set_title("Transcript id")
plot.savefig(output_graph_path)

