import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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

def get_end(reference):
    splitted = reference.split("_")
    end = splitted[len(splitted) - 1]
    return end

def calculate_diff_mod_level(WT_SSS_difference):
    level = ""
    if WT_SSS_difference >= 100:
        level = "Very High increase"
    elif WT_SSS_difference >= 50 and WT_SSS_difference < 100:
        level = "High increase"
    elif WT_SSS_difference >= 20 and WT_SSS_difference < 50:
        level = "Moderate increase"
    elif WT_SSS_difference <= -100:
       level = "Very High decrease"
    elif WT_SSS_difference <= -50 and WT_SSS_difference > -100:
       level = "High decrease"
    elif WT_SSS_difference <= -20 and WT_SSS_difference > -50:
       level = "Moderate decrease"
    elif WT_SSS_difference > -20 and WT_SSS_difference < 0:
       level = "Low decrease"
    else:
        level = "Low increase"
    return level

#input_tail_file_path = "../../data2/tom/RNA/poly-A-significant-changes.tsv"
input_DESeq_file_path = "../../data2/tom/RNA/Yeast_dRNA_WT_vs_dRNA_SSS_results.csv"
input_matrix_file_path = "../../data2/tom/RNA/count_matrix.tsv"
input_RNA_tail_mod_file_path =  "../../data2/tom/RNA/10%_diff-poly-A-RNA-mod_change_level.tsv" #"../../data2/tom/RNA/diff-poly-A-RNA-mod_level.tsv" #poly-A-RNA-mod.tsv"
output_file_path = output_file_path = "../../data2/tom/RNA/expression+RNA_change_10%.tsv" #"../../data2/tom/RNA/expression+poly-A+RNA.tsv"
output_graph_path = "../../data2/tom/RNA/RNA_mod_change_10%_+expression.png" #"../../data2/tom/RNA/cluster+expression+poly-A.png"

df_DESeq = pd.read_csv(input_DESeq_file_path, sep=',', header=0)
df_matrix = pd.read_csv(input_matrix_file_path, sep='\t', header=0)
df_mod_tail = pd.read_csv(input_RNA_tail_mod_file_path, sep='\t', header=0)
#df_DESeq = df_DESeq.dropna()#df_DESeq[df_DESeq['log2FoldChange'].notna()]#df_DESeq.dropna()
#df_matrix = df_matrix.dropna()
df_matrix = df_matrix.loc[df_matrix['transcript_name'].str.contains("mrna")]
df_matrix["transcript_ID_end"] = df_matrix.apply(lambda row : get_end(row["transcript_name"]), axis = 1)
id_dictionary = {}
for ID in df_DESeq['id'].unique():
    end = df_matrix.loc[df_matrix['transcript_ID_end'] == str(ID)]
    array = end['transcript_name'].values[0:]
    if len(array) == 1:
        transcript = array[0]
    id_dictionary[ID] = transcript
#print(id_dictionary)
df_DESeq["transcript_reference"] = df_DESeq.apply(lambda row : match(row["id"], id_dictionary), axis = 1)
df_DESeq["transcript_id"] = df_DESeq.apply(lambda row : extract_transcript_reference(row["transcript_reference"]), axis = 1)
df_DESeq = df_DESeq.loc[(df_DESeq['transcript_id'] != "Failed")]
#print(df_DESeq.head())
#print(df_DESeq.head())
#print(df_matrix.head())
expression_dictionary = {}
for transcript in df_DESeq['transcript_id'].unique():
    data = df_DESeq.loc[df_DESeq['transcript_id'] == str(transcript)]
    #log_change = df_DESeq.loc[df_DESeq['transcript_reference'] == transcript]['log2FoldChange'].values[0]
    log_change = data['log2FoldChange'].values[0:]
    expression_dictionary[transcript] = log_change[0]
#print(expression_dictionary)
#print(df_mod_tail)
df_mod_tail = df_mod_tail.loc[(df_mod_tail["mod_rate_change"].astype(float) <= 400)]
df_mod_tail = df_mod_tail.loc[(df_mod_tail["mod_rate_change"].astype(float) >= -400)]
df_mod_tail["log2FoldChange"] = df_mod_tail.apply(lambda row : match(row["transcript_reference"], expression_dictionary), axis = 1)
#df_mod_tail['diff_mod_level'] = df_mod_tail.apply(lambda row : calculate_diff_mod_level(row['WT_SSS_difference']), axis = 1)
df_mod_tail['diff_mod_level'] = df_mod_tail.apply(lambda row : calculate_diff_mod_level(row["mod_rate_change"]), axis = 1) #relative_average_modification_rate']), axis = 1)
#print(df_mod_tail.head())
df_mod_tail.to_csv(output_file_path, sep='\t', encoding='utf-8')
#want to plot diff_mod_level #WT_SSS_difference
plot = sns.lmplot(y='log2FoldChange', x='mod_rate_change', data=df_mod_tail, hue = "diff_mod_level", fit_reg=False)
fig = plt.gcf()
fig.set_size_inches( 16, 10) #16, 10, 40, 28
plot.set(xlabel = 'Change in differential modification rate', ylabel = 'Log 2 fold change in transcript expression')
plot.set(title = "Change in differential modification rate against change in transcript expression")
plot._legend.set_title("Change in average differential modification rate")#"Number of differential modifications")
plot.savefig(output_graph_path)
