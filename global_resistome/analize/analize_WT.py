from snp_processing.snp_analizer import analyze_sam_record
from aro_ontology import AROntology
from collections import defaultdict
from Bio import SeqIO
import commands
import pysam
import os
import re

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

def parse_sam(samfile_path, gene_info):
    wt_read_count = defaultdict(int)
    wt_snp_count = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {'ref': 0, 'alt': 0})))

    for sam_record in pysam.AlignmentFile(samfile_path):
        gene_id = sam_record.reference_name
        if gene_info[gene_id]['length']:
            read_snp_count = analyze_sam_record(sam_record)
            wt_read_count[gene_id] += 1
            for antibiotic in read_snp_count:
                for snp_pos in read_snp_count[antibiotic]:
                    wt_snp_count[gene_id][antibiotic][snp_pos]['ref'] += read_snp_count[antibiotic][snp_pos]['ref']
                    wt_snp_count[gene_id][antibiotic][snp_pos]['alt'] += read_snp_count[antibiotic][snp_pos]['alt']
    return rd_read_count, wt_read_count, wt_snp_count


def parse_ar_sam(samfile_path, sample):
    samfile = os.path.basename(samfile_path)
    gene_info = defaultdict(lambda: defaultdict())

    RD_fasta_path = os.path.join(SCRIPTDIR, '../fasta/resistance_determinants_genes/resistance_determinants.fasta')
    for record in SeqIO.parse(RD_fasta_path, 'fasta'):
        gene_info[record.id]['set'] = 'RD'
        gene_info[record.id]['length'] = len(record)

    for gene_name in ['folP', 'gyrA', 'gyrB', 'parC', 'parE', 'rpoB']:
        WT_fasta_path = os.path.join(SCRIPTDIR, '../fasta/wildtype_genes/%s.ffn' % gene_name)
        for record in SeqIO.parse(WT_fasta_path, 'fasta'):
            gene_info[record.id]['set'] = 'WT'
            gene_info[record.id]['length'] = len(record)
            gene_info[record.id]['gene_name'] = gene_name

    with open(os.path.join(SCRIPTDIR, '../sample_read_count.tsv')) as f:
        sample_read_count = {sample: int(n) for sample, n in (line.strip().split('\t') for line in f)}
    n_reads = sample_read_count[sample]


    rd_read_count, wt_read_count, wt_snp_count = parse_sam(samfile_path, gene_info)

    aro = AROntology()
    aro.initialize()

    def gene_id_to_antibiotics(gene_id):
        aro_id = re.search('ARO:\d+', gene_id).group(0)
        antibiotics = aro.get_antibiotics_names(aro_id)
        return antibiotics

    rd_antibiotic_arai = {}
    rd_gene_arai = {}
    for gene_id in gene_info:
        if gene_info[gene_id]['set'] == 'RD':
            for antibiotic in gene_id_to_antibiotics(gene_id):
                rd_antibiotic_arai[antibiotic] = 0
            rd_gene_arai[gene_id] = 0

    for gene_id, count in rd_read_count.items():
        length = gene_info[gene_id]['length']
        for antibiotic in gene_id_to_antibiotics(gene_id):
            rd_antibiotic_arai[antibiotic] += 1.0*count/length/n_reads
        rd_gene_arai[gene_id] += 1.0*count/length/n_reads


    # gene_id_resist_percent = {gene_id: {antibiotic: resist_percent}}
    gene_id_resist_percent = defaultdict(lambda: defaultdict(float))
    for gene_id in wt_snp_count:
        for antibiotic in wt_snp_count[gene_id]:
            ref_percent = 1
            for snp_pos in wt_snp_count[gene_id][antibiotic]:
                ref_percent *= wt_snp_count[gene_id][antibiotic][snp_pos]['ref']*1. / (wt_snp_count[gene_id][antibiotic][snp_pos]['ref'] + wt_snp_count[gene_id][antibiotic][snp_pos]['alt'])
            alt_percent = 1 - ref_percent
            gene_id_resist_percent[gene_id][antibiotic] = alt_percent

    # gene_id_resist_percent = {gene_id: {antibiotic: arai}}
    gene_id_mutant_arai = defaultdict(lambda: defaultdict(float))
    for gene_id in gene_id_resist_percent:
        length = gene_info[gene_id]['length']
        for antibiotic in gene_id_resist_percent[gene_id]:
            gene_id_mutant_arai[gene_id][antibiotic] = wt_read_count[gene_id]*1. / length / n_reads
            gene_id_mutant_arai[gene_id][antibiotic] *= gene_id_resist_percent[gene_id][antibiotic]

    # gene_id_resist_percent = {gene_name: {antibiotic: arai}}
    target_genes_mutant_arai = defaultdict(lambda: defaultdict(float))
    for gene_id in gene_id_mutant_arai:
        gene_name = gene_info[gene_id]['gene_name']
        for antibiotic in gene_id_mutant_arai[gene_id]:
            target_genes_mutant_arai[gene_name][antibiotic] += gene_id_mutant_arai[gene_id][antibiotic]

    antibiotics = ['Sulfonamides', 'Fluoroquinolones', 'Coumarin', 'Rifampin']
    target_genes = ['folP', 'gyrA', 'gyrB', 'parC', 'parE', 'rpoB']
    antibiotic_mutant_arai = {antibiotic: 0. for antibiotic in antibiotics}
    for gene_name in target_genes:
        for antibiotic in antibiotics:
            antibiotic_mutant_arai[antibiotic] += target_genes_mutant_arai[gene_name][antibiotic]


    return rd_antibiotic_arai, rd_gene_arai, antibiotic_mutant_arai, target_genes_mutant_arai
