import json
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))
BACMET_INDEX_PATH = os.path.join(SCRIPTDIR, '../references/BacMet/BacMet_EXP.dmnd')
RD_INDEX_PATH = os.path.join(SCRIPTDIR, '../references/resistance_determinants/RD')
WT_INDEX_PATH = os.path.join(SCRIPTDIR, '../references/resistance_mutations/WT')

with open(os.path.join(SCRIPTDIR, '../config.json')) as f:
    EXEC_PATHES = json.load(f)


def map_bacmet(group_name, readfile_path, n_threads, temp_folder):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    daa_file_path = os.path.join(temp_folder, '{}.daa'.format(sample_name))
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', group_name, 'BacMet.{}.sam'.format(sample_name))

    blastx_cmd = '{diamond_path} blastx -d {index} -q {reads} -a {daa} -e 1e-5 -k 1 -p {n_threads}'
    blastx_cmd = blastx_cmd.format(
        diamond_path=EXEC_PATHES['diamond'],
        index=BACMET_INDEX_PATH,
        reads=readfile_path,
        daa=daa_file_path,
        n_threads=n_threads
    )

    view_cmd = '{diamond_path} view -a {daa} -o {sam} -f sam'
    view_cmd = view_cmd.format(
        diamond_path=EXEC_PATHES['diamond'],
        daa=daa_file_path,
        sam=sam_file_path
    )

    os.system(blastx_cmd)
    os.system(view_cmd)
    os.remove(daa_file_path)

def map_RD(group_name, readfile_path, n_threads, temp_folder):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', group_name, 'RD.{}.sam'.format(sample_name))

    bowtie_cmd = '{bowtie_path} -x {index} -U {reads} -S {sam} -k 1 -p {n_threads} --no-unal'
    bowtie_cmd = bowtie_cmd.format(
        bowtie_path=EXEC_PATHES['bowtie2'],
        index=RD_INDEX_PATH,
        reads=readfile_path,
        sam=sam_file_path,
        n_threads=n_threads
    )

    os.system(bowtie_cmd)

def map_WT(group_name, readfile_path, n_threads, temp_folder):
    readfile_name = os.path.basename(readfile_path)
    sample_name, ext = os.path.splitext(readfile_name)
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', group_name, 'WT.{}.sam'.format(sample_name))

    bowtie_cmd = '{bowtie_path} -x {index} -U {reads} -S {sam} -k 1 -p {n_threads} --no-unal'
    bowtie_cmd = bowtie_cmd.format(
        bowtie_path=EXEC_PATHES['bowtie2'],
        index=WT_INDEX_PATH,
        reads=readfile_path,
        sam=sam_file_path,
        n_threads=n_threads
    )

    os.system(bowtie_cmd)

def map_sample(group_name, read_path, n_threads, temp_folder):
    map_bacmet(group_name, read_path, n_threads, temp_folder)
    map_RD(group_name, read_path, n_threads, temp_folder)
    map_WT(group_name, read_path, n_threads, temp_folder)


def map_group(group_name, read_paths_list, n_threads):
    temp_folder = os.path.join(SCRIPTDIR, '../temp')
    os.mkdir(temp_folder)

    group_sam_folder = os.path.join(SCRIPTDIR, '../sam/', group_name)
    os.mkdir(group_sam_folder)

    for read_path in read_paths_list:
        map_sample(group_name, read_path, n_threads, temp_folder)

    os.rmdir(temp_folder)
