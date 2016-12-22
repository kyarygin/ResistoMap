from src.main import main
import argparse
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('read_file', metavar='read_file', type=str, nargs='+',
                        help='input read files')
    parser.add_argument('-n', '--n_threads', type=int,
                        help='number of bowtie/diamond threads (default: 1)', default=1)
    parser.add_argument('-o', '--output_folder', type=str,
                        help='output folder path (default: current dir)', default=os.getcwd())

    args = vars(parser.parse_args())

    read_files_pathes = [os.path.abspath(read_file) for read_file in args['read_file']]
    n_threads = args['n_threads']
    output_folder = os.path.abspath(args['output_folder'])

    main(read_files_pathes, n_threads, output_folder, SCRIPTDIR)
