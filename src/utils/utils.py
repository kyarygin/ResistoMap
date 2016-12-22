from Bio import SeqIO
import pandas as pd
import re
import os

def count_reads(readfile_path):
    n_reads = 0
    for record in SeqIO.parse(readfile_path, 'fastq'):
        n_reads += 1
    return n_reads

def short_gene_id(gene_id):
    return re.search('.*\|([^\|]+)$', gene_id).group(1)

def normalize_name(antibiotic):
    return antibiotic.replace(' antibiotic', '').capitalize()

def substance_dict_to_table(rpkm_substance, sample_name):
    data = [[antibiotic, abund] for antibiotic, abund in rpkm_substance.items()]
    df = pd.DataFrame(data, columns=['category', 'abund'])
    df['sample'] = sample_name
    df = df[['category', 'sample', 'abund']]
    return df

def gene_dict_to_table(rpkm_gene, sample_name):
    data = [[antibiotic, gene_id, abund] for antibiotic, d in rpkm_gene.items() for gene_id, abund in d.items()]
    df = pd.DataFrame(data, columns=['category', 'gene id', 'abund'])
    df['sample'] = sample_name
    df = df[['category', 'gene id', 'sample', 'abund']]
    return df

def delete_folder(top):
    if os.path.isdir(top):
        for root, dirs, files in os.walk(top, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(top)
