from global_resistome.main import main
import argparse

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('read_file', metavar='read_file', type=str, nargs='+',
                        help='read files')
    parser.add_argument('-n', '--n_threads', type=int,
                        help='number of bowtie/diamond threads (default: 20)', default=20)
    args = vars(parser.parse_args())

    read_files = args['read_file']
    n_threads = args['n_threads']

    main(read_files, n_threads)