import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Entrez
pd.options.mode.chained_assignment = None  # default='warn'

def get_gene_name(ref_id):
    Entrez.email = "zx20935@bristol.ac.uk"
    handle = Entrez.efetch(db="gene", id=ref_id, rettype="gb", retmode="text")
    read = handle.read()
    lines = read.splitlines()
    cut_down = lines[2].split("[")
    name = cut_down[0]
    return name

def get_actual_gene_id(ref_id):
    Entrez.email = "zx20935@bristol.ac.uk"
    handle = Entrez.efetch(db="gene", id=ref_id, rettype="gb", retmode="text")
    read = handle.read()
    lines = read.splitlines()
    cut_down = lines[1].split(" ")
    ID = cut_down[1]
    return ID 

def extract_transcript_reference(reference_id):
    if reference_id != None:
        ref_id = reference_id.split("_")
        if len(ref_id) == 6:
            transcript_reference = ref_id[3] +  "_" + ref_id[4]
        else:
            transcript_reference = "Failed"
        return transcript_reference

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

input_DESeq_file_path = "../../data2/tom/RNA/Yeast_dRNA_WT_vs_dRNA_SSS_results.csv"
input_matrix_file_path = "../../data2/tom/RNA/count_matrix.tsv"
#output_file_path1 = "../../data2/tom/RNA/DESeq_upregulated_2.0.tsv"
#output_file_path2 = "../../data2/tom/RNA/DESeq_downregulated_2.0.tsv"
output_file_path3 = "../../data2/tom/RNA/DESeq_0.5.tsv"
output_graph_path = "../../data2/tom/RNA/DESeq_padj_data_graph.png"

df_DESeq = pd.read_csv(input_DESeq_file_path, sep=',', header=0)
#df_DESeq = df_DESeq.loc[(df_DESeq['padj'].astype(float) <= 0.05)]
#df_DESeq = df_DESeq.loc[(df_DESeq["pvalue"].astype(float) <= 0.05)]
df_DESeq = df_DESeq.loc[(df_DESeq['log2FoldChange'].astype(float) <= 1.0) | (df_DESeq['log2FoldChange'].astype(float) >= 0.5)]
df_matrix = pd.read_csv(input_matrix_file_path, sep='\t', header=0)
#df_matrix = df_matrix.loc[df_matrix["transcript_name"].str.contains("mrna")]
df_matrix["transcript_ID_end"] = df_matrix.apply(lambda row : get_end(row["transcript_name"]), axis = 1)
id_dictionary = {}
for ID in df_DESeq['id'].unique():
    end = df_matrix.loc[df_matrix['transcript_ID_end'] == str(ID)]
    array = end['transcript_name'].values[0:]
    if len(array) == 1:
        transcript = array[0]
    id_dictionary[ID] = transcript
df_DESeq["transcript_reference"] = df_DESeq.apply(lambda row : match(row["id"], id_dictionary), axis = 1)
df_DESeq = df_DESeq.loc[df_DESeq["transcript_reference"].str.contains("mrna")]
df_DESeq["transcript_id"] = df_DESeq.apply(lambda row : extract_transcript_reference(row["transcript_reference"]), axis = 1)
df_DESeq = df_DESeq.loc[(df_DESeq["transcript_id"] != "Failed")]
#upregulated_df_DESeq = df_DESeq.loc[(df_DESeq['log2FoldChange'].astype(float) >= 2.0)]
#upregulated_df_DESeq["gene_id"] = upregulated_df_DESeq.apply(lambda row : get_gene_id(row["transcript_id"]), axis = 1)
#upregulated_df_DESeq = upregulated_df_DESeq.drop_duplicates(subset='gene_id', keep="last")
#downregulated_df_DESeq = df_DESeq.loc[(df_DESeq['log2FoldChange'].astype(float) <= -2.0)]
#downregulated_df_DESeq["gene_id"] = downregulated_df_DESeq.apply(lambda row : get_gene_id(row["transcript_id"]), axis = 1)
#downregulated_df_DESeq = downregulated_df_DESeq.drop_duplicates(subset='gene_id', keep="last")

df_DESeq["gene_id"] = df_DESeq.apply(lambda row : get_gene_id(row["transcript_id"]), axis = 1)
df_DESeq = df_DESeq.drop_duplicates(subset='gene_id', keep="last")
df_DESeq["actual_gene_id"] = df_DESeq.apply(lambda row : get_actual_gene_id(row["gene_id"]), axis = 1)
#df_DESeq["gene_name"] = df_DESeq.apply(lambda row : get_gene_name(row["gene_id"]), axis = 1)
#upregulated_df_DESeq.to_csv(output_file_path1, sep='\t', encoding='utf-8')
#downregulated_df_DESeq.to_csv(output_file_path2, sep='\t', encoding='utf-8')
df_DESeq.to_csv(output_file_path3, sep='\t', encoding='utf-8')

#plot = sns.lmplot(x='actual_gene_id', y='log2FoldChange', data=df_DESeq, fit_reg=False)
#plot.set(xlabel = 'Gene ID', ylabel = 'Log2 fold change in expression')#Gaussian mixed model and logistic regression p-values
#plot.set(title = "Change in expression against gene id for significant adjusted p-values")
#fig = plt.gcf()
#fig.set_size_inches(12, 4)
#plot.set_xticklabels(horizontalalignment="right", rotation=90)
#plot.savefig(output_graph_path)
