from global_resistome.mapping.map_utils import map_sample
from global_resistome.utils.utils import count_reads
from global_resistome.analize.analize_bacmet import parse_bacmet_sam
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))


def process_readfile(readfile_path, n_threads):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)

    map_sample(readfile_path, n_threads)

    # n_reads = count_reads(readfile_path)
    # bm_rpkm_substance, bm_rpkm_gene = parse_bacmet_sam(sample_name, n_reads)
    # print(bm_rpkm_gene)

def main(readfile_paths, n_threads):
    os.mkdir(os.path.join(SCRIPTDIR, 'sam'))
    os.mkdir(os.path.join(SCRIPTDIR, 'logs'))
    os.mkdir(os.path.join(SCRIPTDIR, 'temp'))

    for readfile_path in readfile_paths:
        print(readfile_path)
        process_readfile(readfile_path, n_threads)
