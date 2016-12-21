from global_resistome.CARD_ontology.CARD_ontology import Ontology
from global_resistome.utils.utils import short_gene_id, normalize_name
from collections import defaultdict
from Bio import SeqIO
import pysam
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))
TRASH_CATEGORIES = ['Miscellaneous', 'Mixture', 'Acridine dye']
aro = Ontology()
aro.initialize()

def analize_RD_sam(sample_name, n_reads):
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', 'RD.{}.sam'.format(sample_name))

    gene_length = {}
    RD_fasta_path = os.path.join(SCRIPTDIR, '../fasta/RD/CARD_resistance_determinants.fasta')
    for record in SeqIO.parse(RD_fasta_path, 'fasta'):
        gene_length[record.id] = len(record)

    rd_read_count = defaultdict(int)
    for sam_record in pysam.AlignmentFile(sam_file_path):
        gene_id = sam_record.reference_name
        rd_read_count[gene_id] += 1

    rd_arai_antibiotic = defaultdict(float)
    rd_arai_gene = defaultdict(lambda: defaultdict(float))

    for gene_id, count in rd_read_count.items():
        length = gene_length[gene_id]
        antibiotics = aro.gene_id_to_antibiotics(gene_id)
        antibiotics = [normalize_name(antibiotic) for antibiotic in antibiotics]
        for antibiotic in antibiotics:
            rd_arai_antibiotic[antibiotic] += 1.*count / length / n_reads
            rd_arai_gene[antibiotic][short_gene_id(gene_id)] += 1.*count / length / n_reads

    for antibiotic in TRASH_CATEGORIES:
        if antibiotic in rd_arai_gene:
            del rd_arai_gene[antibiotic]
        if antibiotic in rd_arai_antibiotic:
            del rd_arai_antibiotic[antibiotic]

    return rd_arai_antibiotic, rd_arai_gene
