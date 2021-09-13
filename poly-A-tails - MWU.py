import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
pd.options.mode.chained_assignment = None  # default='warn'

def extract_transcript_reference(ref_id):
    """
    Fetches the transcript RefSeq id the ref_id string
    Args:
        ref_id (string): Contains the RefSeq transcript and chromosome ids
    Returns:
        transcript_reference (string): Transcript id
    """
#    print(ref_id)
    transcript_reference = ref_id
    #ref_id = ref_id.split("_")
    #transcript_reference = ref_id[3] + "_" + ref_id[4]
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

input_file_path_WT_R1 = "../../data2/tom/RNA/2020-12_KK_11_Yeast_dRNA_WT_6h_R1/poly-A-nanopolish/polya_tails.tsv" #sys.argv[0]
input_file_path_WT_R2 = "../../data2/tom/RNA/2020-12_KK_12_Yeast_dRNA_WT_6h_R2/poly-A-nanopolish/polya_tails.tsv" #sys.argv[0]
input_file_path_SSS_R1 = "../../data2/tom/RNA/2020-12_KK_09_Yeast_dRNA_SSS_6h_R1/poly-A-nanopolish/polya_tails.tsv" #sys.argv[0]
input_file_path_SSS_R2 = "../../data2/tom/RNA/2020-12_KK_10_Yeast_dRNA_SSS_6h_R2/poly-A-nanopolish/polya_tails.tsv" #sys.argv[0]
output_file_path_WT = "../../data2/tom/RNA/poly-A-WT.tsv"
output_file_path_SSS = "../../data2/tom/RNA/poly-A-SSS.tsv"
output_file_path = "../../data2/tom/RNA/poly-A-MWU.tsv"
output_file_path_filtered = "../../data2/tom/RNA/poly-A-MWU-filtered.tsv"

df_WT_R1 = pd.read_csv(input_file_path_WT_R1, sep='\t', header=0)
df_WT_R2 = pd.read_csv(input_file_path_WT_R2, sep='\t', header=0)
df_WT = pd.concat([df_WT_R1, df_WT_R2])
df_WT = df_WT[df_WT['qc_tag'].astype(str) == "PASS"]
df_WT = df_WT.loc[df_WT["contig"].str.contains("mrna")]
df_WT['transcript_reference'] = df_WT.apply(lambda row : extract_transcript_reference(row['contig']), axis = 1) 

df_SSS_R1 = pd.read_csv(input_file_path_SSS_R1, sep='\t', header=0)
df_SSS_R2 = pd.read_csv(input_file_path_SSS_R2, sep='\t', header=0)
df_SSS = pd.concat([df_SSS_R1, df_SSS_R2])
df_SSS = df_SSS[df_SSS['qc_tag'].astype(str) == "PASS"]
df_SSS = df_SSS.loc[df_SSS['contig'].str.contains("mrna")]
df_SSS['transcript_reference'] = df_SSS.apply(lambda row : extract_transcript_reference(row['contig']), axis = 1) 

#Output the data frame to a .tsv file.
#df_WT.to_csv(output_file_path_WT, sep='\t', encoding='utf-8')
#df_SSS.to_csv(output_file_path_SSS, sep='\t', encoding='utf-8')

df_outputs = pd.DataFrame(columns=('transcript_reference', 'statistic', 'pvalue'))
for transcript in df_WT['transcript_reference'].unique():
    df_1 = df_WT.loc[df_WT['transcript_reference'] == transcript]
    df_2 = df_SSS.loc[df_SSS['transcript_reference'] == transcript]
    statistic, pvalue = mannwhitneyu(df_1['polya_length'], df_2['polya_length']) # Carrying out the Wilcoxon–Mann–Whitney test
#    print(p)
#    print(statistic)
#    print(type(results))
#    statistic = results.query("statistic")
#    pvalue = results.query("pvalue")
#    df_ouputs = df_outputs.append({"transcript_reference": transcript, 'statistic':statistic, 'pvalue': pvalue}, ignore_index=True)
    new_row = pd.Series({"transcript_reference": transcript, "statistic": statistic, "pvalue": pvalue})
    df_outputs = df_outputs.append(new_row, ignore_index = True)
df_outputs.to_csv(output_file_path, sep='\t', encoding='utf-8')
df_outputs_filtered = df_outputs.loc[df_outputs['pvalue'] < 0.05]
df_outputs_filtered.to_csv(output_file_path_filtered, sep='\t', encoding='utf-8')
