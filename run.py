from global_resistome.main import main
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('read_file', metavar='read_file', type=str, nargs='+',
                        help='read files')
    parser.add_argument('-n', '--n_threads', type=int,
                        help='number of bowtie/diamond threads (default: 1)', default=1)
    args = vars(parser.parse_args())

    read_files = args['read_file']
    n_threads = args['n_threads']

    main(read_files, n_threads)