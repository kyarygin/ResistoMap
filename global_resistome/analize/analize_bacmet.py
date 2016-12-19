from collections import defaultdict
from Bio import SeqIO
import pysam
import os
import re

scriptdir = os.path.dirname(os.path.realpath(__file__))

def load_bacmet_fasta(folder_path):
    output = defaultdict(dict)

    bio_chr_path = os.path.join(folder_path, 'BacMet_EXP.BAC0000.Chromosomal.BIOCIDES.288.fasta')
    met_chr_path = os.path.join(folder_path, 'BacMet_EXP.BAC0000.Chromosomal.METALS.283.fasta')
    bio_pla_path = os.path.join(folder_path, 'BacMet_EXP.BAC0000.Plasmid.BIOCIDES.37.fasta')
    met_pla_path = os.path.join(folder_path, 'BacMet_EXP.BAC0000.Plasmid.METALS.161.fasta')

    def process_fasta(path, type_, output):
        for record in SeqIO.parse(path, 'fasta'):
            output[record.id]['set'] = type_
            output[record.id]['length'] = len(record)

    process_fasta(bio_chr_path, 'Biocides', output)
    process_fasta(met_chr_path, 'Metals', output)
    process_fasta(bio_pla_path, 'Biocides', output)
    process_fasta(met_pla_path, 'Metals', output)

    return output

def parse_bacmet_sam(samfile_path, sample):
    bacmet_data = load_bacmet_fasta(os.path.join(scriptdir, '../borya_metals/'))
    counts = {gene_id: 0 for gene_id in bacmet_data}

    with open(os.path.join(scriptdir, '../sample_read_count.tsv')) as f:
        data = (line.strip().split('\t') for line in f)
        sample_sizes = {sample: int(size) for sample, size in data}
    sample_size = sample_sizes[sample]

    with open(samfile_path) as f:
        reader = (line.strip().split('\t') for line in f if not line.startswith('@'))
        for record in reader:
            gene_id = record[2]
            counts[gene_id] += 1

    rpkms = {gene_id: 1.*count/sample_size/bacmet_data[gene_id]['length'] for gene_id, count in counts.items()}
    bm_rpkm_substance = {'Biocides': 0., 'Metals': 0.}
    bm_rpkm_gene = {'Biocides': defaultdict(float), 'Metals': defaultdict(float)}
    for gene_id, rpkm in rpkms.items():
        bm_rpkm_substance[bacmet_data[gene_id]['set']] += rpkm
        bm_rpkm_gene[bacmet_data[gene_id]['set']][gene_id] += rpkm

    return bm_rpkm_substance, bm_rpkm_gene
