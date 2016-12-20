from snp_processing.snp_analizer import analyze_sam_record
from aro_ontology import AROntology
from collections import defaultdict
from Bio import SeqIO
import commands
import pysam
import os
import re

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))


def parse_RD_sam(sample_name, n_reads):
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', 'BacMet.{}.sam'.format(sample_name))

    gene_length = {}
    RD_fasta_path = os.path.join(SCRIPTDIR, '../fasta/RD/CARD_resistance_determinants.fasta')
    for record in SeqIO.parse(RD_fasta_path, 'fasta'):
        gene_length[record.id] = len(record)

    rd_read_count = defaultdict(int)
    for sam_record in pysam.AlignmentFile(samfile_path):
        gene_id = sam_record.reference_name
            rd_read_count[gene_id] += 1

    aro = AROntology()
    aro.initialize()

    def gene_id_to_antibiotics(gene_id):
        aro_id = re.search('ARO:\d+', gene_id).group(0)
        antibiotics = aro.get_antibiotics_names(aro_id)
        return antibiotics

    rd_antibiotic_rpkm = {}
    rd_gene_rpkm = {}
    for gene_id in gene_info:
        if gene_info[gene_id]['set'] == 'RD':
            for antibiotic in gene_id_to_antibiotics(gene_id):
                rd_antibiotic_rpkm[antibiotic] = 0
            rd_gene_rpkm[gene_id] = 0

    for gene_id, count in rd_read_count.items():
        length = gene_info[gene_id]['length']
        for antibiotic in gene_id_to_antibiotics(gene_id):
            rd_antibiotic_rpkm[antibiotic] += 1.0*count/length/n_reads
        rd_gene_rpkm[gene_id] += 1.0*count/length/n_reads

    return rd_antibiotic_rpkm, rd_gene_rpkm
