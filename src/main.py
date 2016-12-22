from global_resistome.mapping.map_utils import map_sample
from global_resistome.utils.utils import delete_folder, count_reads
from global_resistome.utils.utils import substance_dict_to_table, gene_dict_to_table
from global_resistome.analize.analize_bacmet import analize_bacmet_sam
from global_resistome.analize.analize_RD import analize_RD_sam
from global_resistome.analize.analize_WT import analize_WT_sam
import pandas as pd
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

def concat_tables():
    tables_names = os.listdir(os.path.join(SCRIPTDIR, 'tables'))
    gene_tables_names = [x for x in tables_names if x.startswith('gene_table')]
    ab_tables_names = [x for x in tables_names if x.startswith('ab_table')]

    gene_table_total = pd.DataFrame()
    for table_name in gene_tables_names:
        gene_table = pd.read_csv(os.path.join(SCRIPTDIR, 'tables', table_name), sep='\t')
        gene_table_total = pd.concat([gene_table_total, gene_table])

    ab_table_total = pd.DataFrame()
    for table_name in ab_tables_names:
        ab_table = pd.read_csv(os.path.join(SCRIPTDIR, 'tables', table_name), sep='\t')
        ab_table_total = pd.concat([ab_table_total, ab_table])

    with open(os.path.join(SCRIPTDIR, '..', 'gene_table_total.tsv'), 'w') as f:
        gene_table_total.to_csv(f, sep='\t', index=False, float_format='%.2e')
    with open(os.path.join(SCRIPTDIR, '..', 'ab_table_total.tsv'), 'w') as f:
        ab_table_total.to_csv(f, sep='\t', index=False, float_format='%.2e')

def process_readfile(readfile_path, n_threads):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    n_reads = count_reads(readfile_path)

    map_sample(readfile_path, n_threads)

    bm_rpkm_substance, bm_rpkm_gene = analize_bacmet_sam(sample_name, n_reads)
    rd_rpkm_antibiotic, rd_rpkm_gene = analize_RD_sam(sample_name, n_reads)
    wt_rpkm_antibiotic, wt_rpkm_gene = analize_WT_sam(sample_name, n_reads)

    bm_rpkm_substance = substance_dict_to_table(bm_rpkm_substance, sample_name)
    rd_rpkm_antibiotic = substance_dict_to_table(rd_rpkm_antibiotic, sample_name)
    wt_rpkm_antibiotic = substance_dict_to_table(wt_rpkm_antibiotic, sample_name)

    ab_table = pd.concat([bm_rpkm_substance, rd_rpkm_antibiotic, wt_rpkm_antibiotic])
    with open(os.path.join(SCRIPTDIR, 'tables', 'ab_table.{}.tsv'.format(sample_name)), 'w') as f:
        ab_table.to_csv(f, sep='\t', index=False, float_format='%.2e')

    bm_rpkm_gene = gene_dict_to_table(bm_rpkm_gene, sample_name)
    rd_rpkm_gene = gene_dict_to_table(rd_rpkm_gene, sample_name)
    wt_rpkm_gene = gene_dict_to_table(wt_rpkm_gene, sample_name)

    gene_table = pd.concat([bm_rpkm_gene, rd_rpkm_gene, wt_rpkm_gene])
    with open(os.path.join(SCRIPTDIR, 'tables', 'gene_table.{}.tsv'.format(sample_name)), 'w') as f:
        gene_table.to_csv(f, sep='\t', index=False, float_format='%.2e')

def main(readfile_paths, n_threads):
    os.mkdir(os.path.join(SCRIPTDIR, 'sam'))
    os.mkdir(os.path.join(SCRIPTDIR, 'logs'))
    os.mkdir(os.path.join(SCRIPTDIR, 'temp'))
    os.mkdir(os.path.join(SCRIPTDIR, 'tables'))

    for readfile_path in readfile_paths:
        print('Processing {} ...'.format(os.path.basename(readfile_path)))
        process_readfile(readfile_path, n_threads)

    concat_tables()

    delete_folder(os.path.join(SCRIPTDIR, 'sam'))
    delete_folder(os.path.join(SCRIPTDIR, 'logs'))
    delete_folder(os.path.join(SCRIPTDIR, 'temp'))
    delete_folder(os.path.join(SCRIPTDIR, 'tables'))
