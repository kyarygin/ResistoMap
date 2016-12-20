from collections import namedtuple, defaultdict
from itertools import count
from Bio import SeqIO
import os
scriptdir = os.path.dirname(os.path.realpath(__file__))

def load_snp_table(snp_table_path):
    headers = ['gene_name', 'pos', 'ref', 'alt', 'antibiotic', 'cite_ref']
    coltypes = [str, int, str, str, str, str, str]
    snp_record = namedtuple('snp_record', headers)
    with open(snp_table_path) as f:
        raw_table = (line.strip().split(',') for line in f)
        raw_table = ([coltype(value) for coltype, value in zip(coltypes, rec)] for rec in raw_table)
        snp_table = (snp_record(*rec) for rec in raw_table)
        snp_table = [rec for rec in snp_table if rec.gene_name != 'embB'] # filter
    return snp_table

def process_aligned_faa(faa_folder, ecoli_genes_id, snp_table):
    gene_snp = {}
    for gene_name in ['folP', 'gyrA', 'gyrB', 'parC', 'parE', 'rpoB']:
        ecoli_snp = defaultdict(list)
        pos_to_antibiotic = defaultdict(list)
        for rec in snp_table:
            if rec.gene_name == gene_name:
                ecoli_snp[rec.pos].append((rec.ref, rec.alt))
                pos_to_antibiotic[rec.pos] = rec.antibiotic

        path = os.path.join(faa_folder, '%s.faa.aligned' % gene_name)
        stack = {record.id: record.seq._data for record in SeqIO.parse(path, 'fasta')}

        ecoli_seq = stack[ecoli_genes_id[gene_name]]

        for gene_id, gene_seq in stack.items():
            ecoli_counter = count(start=1)
            gene_counter = count(start=1)
            ecoli_enum = [next(ecoli_counter) if aa != '-' else 0 for aa in ecoli_seq]
            gene_enum = [next(gene_counter) if aa != '-' else 0 for aa in gene_seq]
            ecoli_pos_to_gene_pos = {x: y for x, y in zip(ecoli_enum, gene_enum) if x}

            gene_snp[gene_id] = {ecoli_pos_to_gene_pos[x]: {'altref_list': ecoli_snp[x][:],
                                                            'antibiotic': pos_to_antibiotic[x]} for x in ecoli_snp}

    return gene_snp


def get_gene_snp():
    ecoli_genes_path = os.path.join(scriptdir, './data/ecoli_genes.csv')
    with open(ecoli_genes_path) as f:
        ecoli_genes_id = dict(line.strip().split(',') for line in f)

    snp_table_path = os.path.join(scriptdir, './data/AR_SNP_table.csv')
    snp_table = load_snp_table(snp_table_path)
    aligned_faa_path = os.path.join(scriptdir, './data/aligned_faa/')
    gene_snp = process_aligned_faa(aligned_faa_path, ecoli_genes_id, snp_table)

    return gene_snp

