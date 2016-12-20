from global_resistome.mapping.map_utils import map_sample
from global_resistome.utils.utils import count_reads
from global_resistome.analize.analize_bacmet import analize_bacmet_sam
from global_resistome.analize.analize_RD import analize_RD_sam
from global_resistome.analize.analize_WT import analize_WT_sam
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))


def process_readfile(readfile_path, n_threads):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    n_reads = count_reads(readfile_path)

    map_sample(readfile_path, n_threads)

    bm_rpkm_substance, bm_rpkm_gene = analize_bacmet_sam(sample_name, n_reads)
    rd_rpkm_antibiotic, rd_rpkm_gene = analize_RD_sam(sample_name, n_reads)
    rd_rpkm_antibiotic, rd_rpkm_gene = analize_WT_sam(sample_name, n_reads)


def main(readfile_paths, n_threads):
    # os.mkdir(os.path.join(SCRIPTDIR, 'sam'))
    # os.mkdir(os.path.join(SCRIPTDIR, 'logs'))
    # os.mkdir(os.path.join(SCRIPTDIR, 'temp'))

    for readfile_path in readfile_paths:
        print(readfile_path)
        process_readfile(readfile_path, n_threads)
