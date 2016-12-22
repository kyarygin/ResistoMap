from global_resistome.analize.snp_processing.genes_snp import get_gene_snp
from collections import defaultdict
from Bio.Seq import Seq

gene_snp = get_gene_snp()


def extract_snp_data(sam_record):
    s_pos = sam_record.pos
    e_pos = s_pos + sam_record.query_length - 1
    s_pos_n = ((s_pos-1)//3 + 1) * 3
    e_pos_n = ((e_pos+1)//3 * 3) - 1

    query_n = sam_record.query[s_pos_n-s_pos:]
    if e_pos != e_pos_n:
        query_n = query_n[:-(e_pos-e_pos_n)]

    aa_query = Seq(query_n).translate()
    s_aa_pos = s_pos_n//3 + 1
    e_aa_pos = s_aa_pos + len(aa_query) - 1

    read_snp_count = defaultdict(lambda: defaultdict(lambda: {'ref': 0, 'alt': 0}))
    # read_snp_count = {antibiotic: {position: {'ref': 0, 'alt': 0}}}

    gene_id = sam_record.reference_name
    for snp_pos in gene_snp[gene_id]:
        if s_aa_pos <= snp_pos <= e_aa_pos:
            antibiotic = gene_snp[gene_id][snp_pos]['antibiotic']
            ref_aa = {aa for aa, _ in gene_snp[gene_id][snp_pos]['altref_list']}
            alt_aa = {aa for _, aa in gene_snp[gene_id][snp_pos]['altref_list']}
            if aa_query[snp_pos-s_aa_pos] in ref_aa:
                read_snp_count[antibiotic][snp_pos]['ref'] += 1
            if aa_query[snp_pos-s_aa_pos] in alt_aa:
                read_snp_count[antibiotic][snp_pos]['alt'] += 1
    return read_snp_count

