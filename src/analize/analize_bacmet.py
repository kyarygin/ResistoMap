from src.utils.utils import short_gene_id
from collections import defaultdict
from Bio import SeqIO
import os
import re

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

def load_bacmet_fasta():
    output = defaultdict(dict)

    bio_path = os.path.join(SCRIPTDIR, '../fasta/BacMet/BacMet_EXP.BIOCIDES.fasta')
    met_path = os.path.join(SCRIPTDIR, '../fasta/BacMet/BacMet_EXP.METALS.fasta')

    def process_fasta(path, type_, output):
        for record in SeqIO.parse(path, 'fasta'):
            output[record.id]['set'] = type_
            output[record.id]['length'] = len(record)

    process_fasta(bio_path, 'Biocides', output)
    process_fasta(met_path, 'Metals', output)

    return output

def analize_bacmet_sam(sample_name, n_reads):
    bacmet_data = load_bacmet_fasta()

    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', 'BacMet.{}.sam'.format(sample_name))

    counts = defaultdict(int)
    with open(sam_file_path) as f:
        reader = (line.strip().split('\t') for line in f if not line.startswith('@'))
        for record in reader:
            gene_id = record[2]
            counts[gene_id] += 1

    arais = {gene_id: 1. * count / n_reads / bacmet_data[gene_id]['length'] for gene_id, count in counts.items()}
    bm_arai_substance = {'Biocides': 0., 'Metals': 0.}
    bm_arai_gene = {'Biocides': defaultdict(float), 'Metals': defaultdict(float)}
    for gene_id, arai in arais.items():
        bm_arai_substance[bacmet_data[gene_id]['set']] += arai
        bm_arai_gene[bacmet_data[gene_id]['set']][short_gene_id(gene_id)] += arai

    return bm_arai_substance, bm_arai_gene
