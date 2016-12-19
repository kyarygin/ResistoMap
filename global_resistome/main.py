from global_resistome.mapping.map_utils import map_sample
from global_resistome.utils.utils import count_reads

def process_readfile(readfile_path, n_threads):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    # map_sample(readfile_path, n_threads)
    n_reads = count_reads(readfile_path)
    print(n_reads)

def main(readfile_paths, n_threads):
    for readfile_path in readfile_paths:
        print(readfile_path)
        process_readfile(readfile_path, n_threads)