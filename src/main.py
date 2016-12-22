from src.utils.utils import delete_folder, count_reads
from src.utils.utils import substance_dict_to_table, gene_dict_to_table
from src.mapping.map_utils import map_sample
from src.analize.analize_bacmet import analize_bacmet_sam
from src.analize.analize_RD import analize_RD_sam
from src.analize.analize_WT import analize_WT_sam
import pandas as pd
import os

def concat_tables(root_path):
    tables_names = os.listdir(os.path.join(root_path, 'src', 'output'))
    gene_tables_names = [x for x in tables_names if x.startswith('gene_table')]
    ab_tables_names = [x for x in tables_names if x.startswith('ab_table')]

    gene_table_total = pd.DataFrame()
    for table_name in gene_tables_names:
        gene_table = pd.read_csv(os.path.join(root_path, 'src', 'output', table_name), sep='\t')
        gene_table_total = pd.concat([gene_table_total, gene_table])

    ab_table_total = pd.DataFrame()
    for table_name in ab_tables_names:
        ab_table = pd.read_csv(os.path.join(root_path, 'src', 'output', table_name), sep='\t')
        ab_table_total = pd.concat([ab_table_total, ab_table])

    with open(os.path.join(root_path, 'gene_table_total.tsv'), 'w') as f:
        gene_table_total.to_csv(f, sep='\t', index=False, float_format='%.2e')
    with open(os.path.join(root_path, 'ab_table_total.tsv'), 'w') as f:
        ab_table_total.to_csv(f, sep='\t', index=False, float_format='%.2e')

def process_readfile(readfile_path, n_threads, root_path):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    n_reads = count_reads(readfile_path)

    map_sample(readfile_path, n_threads, root_path)

    bm_rpkm_substance, bm_rpkm_gene = analize_bacmet_sam(sample_name, n_reads, root_path)
    rd_rpkm_substance, rd_rpkm_gene = analize_RD_sam(sample_name, n_reads, root_path)
    wt_rpkm_substance, wt_rpkm_gene = analize_WT_sam(sample_name, n_reads, root_path)

    bm_rpkm_substance = substance_dict_to_table(bm_rpkm_substance, sample_name)
    rd_rpkm_substance = substance_dict_to_table(rd_rpkm_substance, sample_name)
    wt_rpkm_substance = substance_dict_to_table(wt_rpkm_substance, sample_name)

    ab_table = pd.concat([bm_rpkm_substance, rd_rpkm_substance, wt_rpkm_substance])
    with open(os.path.join(root_path, 'src', 'output', 'ab_table.{}.tsv'.format(sample_name)), 'w') as f:
        ab_table.to_csv(f, sep='\t', index=False, float_format='%.2e')

    bm_rpkm_gene = gene_dict_to_table(bm_rpkm_gene, sample_name)
    rd_rpkm_gene = gene_dict_to_table(rd_rpkm_gene, sample_name)
    wt_rpkm_gene = gene_dict_to_table(wt_rpkm_gene, sample_name)

    gene_table = pd.concat([bm_rpkm_gene, rd_rpkm_gene, wt_rpkm_gene])
    with open(os.path.join(root_path, 'src', 'output', 'gene_table.{}.tsv'.format(sample_name)), 'w') as f:
        gene_table.to_csv(f, sep='\t', index=False, float_format='%.2e')

def main(readfile_paths, n_threads, root_path):
    temp_folders = ['sam', 'logs', 'daa', 'output']
    for folder_name in temp_folders:
        delete_folder(os.path.join(root_path, 'src', folder_name))
        os.mkdir(os.path.join(root_path, 'src', folder_name))

    for readfile_path in readfile_paths:
        print('Processing {} ...'.format(os.path.basename(readfile_path)))
        process_readfile(readfile_path, n_threads, root_path)

    concat_tables(root_path)

    for folder_name in temp_folders:
        delete_folder(os.path.join(root_path, 'src', folder_name))
