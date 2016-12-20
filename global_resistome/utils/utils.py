from Bio import SeqIO
import re

def count_reads(readfile_path):
    n_reads = 0
    for record in SeqIO.parse(readfile_path, 'fastq'):
        n_reads += 1
    return n_reads

def short_gene_id(gene_id):
    return re.search('.*\|([^\|]+)$', gene_id).group(1)