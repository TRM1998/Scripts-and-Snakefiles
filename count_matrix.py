import sys
import pandas as pd
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None  # default='warn'

def extract_transcript_reference(ref_id):
    """
    Fetches the transcript RefSeq id the ref_id string
    Args:
        ref_id (string): Contains the RefSeq transcript and chromosome ids
    Returns:
        transcript_reference (string): Transcript id
    """
    ref_id = ref_id.split("_")
    transcript_reference = ref_id[3] + "_" + ref_id[4]
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
    if key in dictionary:# and key != None:
        matched = dictionary.get(key)
    else:
        matched = 0
    return matched


def get_counts(df):
    dictionary = {}
    for transcript in df['transcript_name']:
        data =  df.loc[df['transcript_name'] == transcript, ['est_count']].iloc[0]
        data = data.get("est_count")
        dictionary[transcript] = data
#    dictionary = {k: v for k, v in dictionary.items() if v is not None}
#    dictionary = {k: v for k, v in dictionary.items() if k is not None}
    return dictionary

input_file_path_WT_R1 = "../../data2/tom/RNA/2020-12_KK_11_Yeast_dRNA_WT_6h_R1/NanoCount/transcript_counts.tsv" #sys.argv[0]
input_file_path_WT_R2 = "../../data2/tom/RNA/2020-12_KK_12_Yeast_dRNA_WT_6h_R2/NanoCount/transcript_counts.tsv" #sys.argv[0]
input_file_path_SSS_R1 = "../../data2/tom/RNA/2020-12_KK_09_Yeast_dRNA_SSS_6h_R1/NanoCount/transcript_counts.tsv" #sys.argv[0]
input_file_path_SSS_R2 = "../../data2/tom/RNA/2020-12_KK_10_Yeast_dRNA_SSS_6h_R2/NanoCount/transcript_counts.tsv" #sys.argv[0]
output_file_path = "../../data2/tom/RNA/count_matrix.tsv"

df_WT_R1 = pd.read_csv(input_file_path_WT_R1, sep='\t', header=0)
df_WT_R2 = pd.read_csv(input_file_path_WT_R2, sep='\t', header=0)
df_SSS_R1 = pd.read_csv(input_file_path_SSS_R1, sep='\t', header=0)
df_SSS_R2 = pd.read_csv(input_file_path_SSS_R2, sep='\t', header=0)
df = pd.concat([df_WT_R1, df_WT_R2, df_SSS_R1, df_SSS_R2])
#df = df.reindex(["transcript_name"], axis = 1)
#df = df.loc[df['transcript_name'].unique(), axis = 1]
df = df.drop_duplicates(subset= ["transcript_name"])
df = df.drop(columns=['raw', 'est_count', 'tpm'], axis=1)
#print(df.head())

dictionary_WT_R1 = get_counts(df_WT_R1)
#print(dictionary_WT_R1)
df['WT_R1_count'] = df.apply(lambda row : match(row['transcript_name'], dictionary_WT_R1), axis = 1) 
dictionary_WT_R2 = get_counts(df_WT_R2)
print(dictionary_WT_R2.get("lcl|NC_001139.9_mrna_NM_001181423.1_2557"))
df['WT_R2_count'] = df.apply(lambda row : match(row['transcript_name'], dictionary_WT_R2), axis = 1) 
dictionary_SSS_R1 = get_counts(df_SSS_R1)
df['SSS_R1_count'] = df.apply(lambda row : match(row['transcript_name'], dictionary_SSS_R1), axis = 1) 
dictionary_SSS_R2 = get_counts(df_SSS_R2)
df['SSS_R2_count'] = df.apply(lambda row : match(row['transcript_name'], dictionary_SSS_R2), axis = 1) 

#Output the data frame to a .tsv file.
df.to_csv(output_file_path, sep='\t', encoding='utf-8')
