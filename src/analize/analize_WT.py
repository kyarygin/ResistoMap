from global_resistome.analize.snp_processing.snp_analizer import extract_snp_data
from collections import defaultdict
from Bio import SeqIO
import pysam
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))
GENE_NAMES = ['folP', 'gyrA', 'gyrB', 'parC', 'parE', 'rpoB']
ANTIBIOTICS = ['Sulfonamides', 'Fluoroquinolones', 'Coumarin', 'Rifampin']

def parse_WT_sam(sam_file_path, gene_info):
    wt_read_count = defaultdict(int)
    wt_snp_count = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {'ref': 0, 'alt': 0})))

    for sam_record in pysam.AlignmentFile(sam_file_path):
        gene_id = sam_record.reference_name
        if gene_info[gene_id]['length']:
            read_snp_count = extract_snp_data(sam_record)
            wt_read_count[gene_id] += 1
            for antibiotic in read_snp_count:
                for snp_pos in read_snp_count[antibiotic]:
                    wt_snp_count[gene_id][antibiotic][snp_pos]['ref'] += read_snp_count[antibiotic][snp_pos]['ref']
                    wt_snp_count[gene_id][antibiotic][snp_pos]['alt'] += read_snp_count[antibiotic][snp_pos]['alt']
    return wt_read_count, wt_snp_count


def analize_WT_sam(sample_name, n_reads):
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', 'WT.{}.sam'.format(sample_name))
    gene_info = defaultdict(lambda: defaultdict())

    for gene_name in GENE_NAMES:
        WT_fasta_path = os.path.join(SCRIPTDIR, '../fasta/WT/%s.ffn' % gene_name)
        for record in SeqIO.parse(WT_fasta_path, 'fasta'):
            gene_info[record.id]['length'] = len(record)
            gene_info[record.id]['gene_name'] = gene_name

    wt_read_count, wt_snp_count = parse_WT_sam(sam_file_path, gene_info)

    # gene_id_resist_percent = {gene_id: {antibiotic: resist_percent}}
    gene_id_resist_percent = defaultdict(lambda: defaultdict(float))
    for gene_id in wt_snp_count:
        for antibiotic in wt_snp_count[gene_id]:
            ref_percent = 1
            for snp_pos in wt_snp_count[gene_id][antibiotic]:
                ref_percent *= wt_snp_count[gene_id][antibiotic][snp_pos]['ref']*1. / (wt_snp_count[gene_id][antibiotic][snp_pos]['ref'] + wt_snp_count[gene_id][antibiotic][snp_pos]['alt'])
            alt_percent = 1 - ref_percent
            gene_id_resist_percent[gene_id][antibiotic] = alt_percent

    # gene_id_mutant_arai = {gene_id: {antibiotic: arai}}
    gene_id_mutant_arai = defaultdict(lambda: defaultdict(float))
    for gene_id in gene_id_resist_percent:
        length = gene_info[gene_id]['length']
        for antibiotic in gene_id_resist_percent[gene_id]:
            gene_id_mutant_arai[gene_id][antibiotic] = wt_read_count[gene_id]*1. / length / n_reads
            gene_id_mutant_arai[gene_id][antibiotic] *= gene_id_resist_percent[gene_id][antibiotic]

    # target_genes_mutant_arai = {gene_name: {antibiotic: arai}}
    target_genes_mutant_arai = defaultdict(lambda: defaultdict(float))
    for gene_id in gene_id_mutant_arai:
        gene_name = gene_info[gene_id]['gene_name']
        for antibiotic in gene_id_mutant_arai[gene_id]:
            target_genes_mutant_arai[gene_name][antibiotic] += gene_id_mutant_arai[gene_id][antibiotic]

    # target_genes_mutant_arai_T = {antibiotic: {gene_name: arai}}
    target_genes_mutant_arai_T = defaultdict(lambda: defaultdict(float))
    for gene_name in target_genes_mutant_arai:
        for antibiotic in target_genes_mutant_arai[gene_name]:
            target_genes_mutant_arai_T[antibiotic][gene_name] += target_genes_mutant_arai[gene_name][antibiotic]


    antibiotic_mutant_arai = {antibiotic: 0. for antibiotic in ANTIBIOTICS}
    for gene_name in GENE_NAMES:
        for antibiotic in ANTIBIOTICS:
            antibiotic_mutant_arai[antibiotic] += target_genes_mutant_arai[gene_name][antibiotic]


    return antibiotic_mutant_arai, target_genes_mutant_arai_T
