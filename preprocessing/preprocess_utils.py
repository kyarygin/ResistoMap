from Bio import SeqIO
import os

def count_reads(readfile_path):
    n_reads = 0
    for record in SeqIO.parse(readfile_path, 'fastq'):
        n_reads += 1
    print(n_reads)
