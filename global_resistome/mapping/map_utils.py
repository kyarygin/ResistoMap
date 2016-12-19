import json
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))
BACMET_INDEX_PATH = os.path.join(SCRIPTDIR, '../references/BacMet/BacMet_EXP.dmnd')
RD_INDEX_PATH = os.path.join(SCRIPTDIR, '../references/resistance_determinants/RD')
WT_INDEX_PATH = os.path.join(SCRIPTDIR, '../references/resistance_mutations/WT')

with open(os.path.join(SCRIPTDIR, '../../config.json')) as f:
    EXEC_PATHES = json.load(f)


def map_bacmet(readfile_path, n_threads):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    daa_file_path = os.path.join(SCRIPTDIR, '../temp/', '{}.daa'.format(sample_name))
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', 'BacMet.{}.sam'.format(sample_name))
    log_file_path = os.path.join(SCRIPTDIR, '../logs/', 'BacMet.mapping.{}.log'.format(sample_name))

    blastx_cmd = '{diamond_path} blastx -d {index} -q {reads} -a {daa} -e 1e-5 -k 1 -p {n_threads} > {log}'
    blastx_cmd = blastx_cmd.format(
        diamond_path=EXEC_PATHES['diamond'],
        index=BACMET_INDEX_PATH,
        reads=readfile_path,
        daa=daa_file_path,
        n_threads=n_threads,
        log=log_file_path
    )

    log_file_path = os.path.join(SCRIPTDIR, '../logs/', 'BacMet.view.{}.log'.format(sample_name))
    view_cmd = '{diamond_path} view -a {daa} -o {sam} -f sam > {log}'
    view_cmd = view_cmd.format(
        diamond_path=EXEC_PATHES['diamond'],
        daa=daa_file_path,
        sam=sam_file_path,
        log=log_file_path
    )

    os.system(blastx_cmd)
    os.system(view_cmd)
    os.remove(daa_file_path)

def map_RD(readfile_path, n_threads):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', 'RD.{}.sam'.format(sample_name))
    log_file_path = os.path.join(SCRIPTDIR, '../logs/', 'RD.mapping.{}.log'.format(sample_name))

    bowtie_cmd = '{bowtie_path} -x {index} -U {reads} -S {sam} -k 1 -p {n_threads} --no-unal > {log} 2>{log}'
    bowtie_cmd = bowtie_cmd.format(
        bowtie_path=EXEC_PATHES['bowtie2'],
        index=RD_INDEX_PATH,
        reads=readfile_path,
        sam=sam_file_path,
        n_threads=n_threads,
        log=log_file_path
    )

    os.system(bowtie_cmd)

def map_WT(readfile_path, n_threads):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', 'WT.{}.sam'.format(sample_name))
    log_file_path = os.path.join(SCRIPTDIR, '../logs/', 'WT.mapping.{}.log'.format(sample_name))

    bowtie_cmd = '{bowtie_path} -x {index} -U {reads} -S {sam} -k 1 -p {n_threads} --no-unal > {log} 2>{log}'
    bowtie_cmd = bowtie_cmd.format(
        bowtie_path=EXEC_PATHES['bowtie2'],
        index=WT_INDEX_PATH,
        reads=readfile_path,
        sam=sam_file_path,
        n_threads=n_threads,
        log=log_file_path
    )

    os.system(bowtie_cmd)

def map_sample(readfile_path, n_threads):
    map_bacmet(readfile_path, n_threads)
    map_RD(readfile_path, n_threads)
    map_WT(readfile_path, n_threads)
