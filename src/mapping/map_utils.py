import json
import os

def map_bacmet(readfile_path, exec_path_diamond, n_threads, root_path):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)

    bacmet_index_path = os.path.join(root_path, 'src', 'references', 'BacMet', 'BacMet_EXP.dmnd')

    sam_file_path = os.path.join(root_path, 'src', 'sam', 'BacMet.{}.sam'.format(sample_name))
    log_file_path = os.path.join(root_path, 'src', 'logs', 'BacMet.mapping.{}.log'.format(sample_name))

    blastx_cmd = '{diamond_path} blastx -d {index} -q {reads} -o {sam} -f 101 -e 1e-5 -k 1 -p {n_threads} > {log}'
    blastx_cmd = blastx_cmd.format(
        diamond_path=exec_path_diamond,
        index=bacmet_index_path,
        reads=readfile_path,
        sam=sam_file_path,
        n_threads=n_threads,
        log=log_file_path
    )

    os.system(blastx_cmd)

def map_RD(readfile_path, exec_path_bowtie2, n_threads, root_path):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)

    rd_index_path = os.path.join(root_path, 'src', 'references', 'RD', 'RD')

    sam_file_path = os.path.join(root_path, 'src', 'sam', 'RD.{}.sam'.format(sample_name))
    log_file_path = os.path.join(root_path, 'src', 'logs', 'RD.mapping.{}.log'.format(sample_name))

    bowtie_cmd = '{bowtie_path} -x {index} -U {reads} -S {sam} -k 1 -p {n_threads} --no-unal > {log} 2>{log}'
    bowtie_cmd = bowtie_cmd.format(
        bowtie_path=exec_path_bowtie2,
        index=rd_index_path,
        reads=readfile_path,
        sam=sam_file_path,
        n_threads=n_threads,
        log=log_file_path
    )

    os.system(bowtie_cmd)

def map_WT(readfile_path, exec_path_bowtie2, n_threads, root_path):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)

    wt_index_path = os.path.join(root_path, 'src', 'references', 'WT', 'WT')

    sam_file_path = os.path.join(root_path, 'src', 'sam', 'WT.{}.sam'.format(sample_name))
    log_file_path = os.path.join(root_path, 'src', 'logs', 'WT.mapping.{}.log'.format(sample_name))

    bowtie_cmd = '{bowtie_path} -x {index} -U {reads} -S {sam} -k 1 -p {n_threads} --no-unal > {log} 2>{log}'
    bowtie_cmd = bowtie_cmd.format(
        bowtie_path=exec_path_bowtie2,
        index=wt_index_path,
        reads=readfile_path,
        sam=sam_file_path,
        n_threads=n_threads,
        log=log_file_path
    )

    os.system(bowtie_cmd)

def map_sample(readfile_path, n_threads, root_path):
    with open(os.path.join(root_path, 'config.json')) as f:
        exec_pathes = json.load(f)
    exec_path_bowtie2 = exec_pathes['bowtie2']
    exec_path_diamond = exec_pathes['diamond']

    map_bacmet(readfile_path, exec_path_diamond, n_threads, root_path)
    map_RD(readfile_path, exec_path_bowtie2, n_threads, root_path)
    map_WT(readfile_path, exec_path_bowtie2, n_threads, root_path)
