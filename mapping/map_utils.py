import json
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))
BACMET_INDEX_PATH = os.path.join(SCRIPTDIR, '../references/BacMet/BacMet_EXP.dmnd')
RD_INDEX_PATH = os.path.join(SCRIPTDIR, '../references/resistance_determinants/RD')
WT_INDEX_PATH = os.path.join(SCRIPTDIR, '../references/resistance_mutations/WT')

with open(os.path.join(SCRIPTDIR, '../config.json')) as f:
    EXEC_PATHES = json.load(f)


def map_bacmet(group_name, readfile_path, n_threads, temp_folder):
    readfile_name = os.path.basename(reads_path)
    sample_name, ext = os.path.splitext(readfile_name)
    daa_file_path = os.path.join(temp_folder, '{}.daa'.format(sample_name))
    sam_file_path = os.path.join(SCRIPTDIR, '../sam/', '{}.BacMet.sam'.format(sample_name))

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




# def map(group_name, read_path, n_threads):
#     temp_folder = os.path.join(SCRIPTDIR, '../temp')
#     os.mkdir(temp_folder)
#     map_bacmet(group_name, read_path, n_threads, temp_folder)


# card_index_path="$1"
# sample="$3"
# merged_sample_path="$4"
# group_name="$5"
# map_card="$6"
# map_bacmet="$7"

# card_sam_path="./sam/$group_name/RD_WT__$sample.sam"
# bacmet_sam_path="./sam/$group_name/BacMet__$sample.sam"

# if [ "$map_card" == "True" ]; then
#     bowtie2 -x "$card_index_path" -U "$merged_sample_path" -S "$card_sam_path" -k 1 -p 20 --no-unal;
# fi
# if [ "$map_bacmet" == "True" ]; then
#     ./diamond blastx -d "$bacmet_index_path" -q "$merged_sample_path" -a "./temp/$sample" -e 1e-5 -k 1 -p 20;
#     ./diamond view -a "./temp/$sample.daa" -o "$bacmet_sam_path" -f sam;
#     rm "./temp/$sample.daa";
# fi

