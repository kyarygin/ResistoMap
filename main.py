import argparse
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('group_name', type=str,
                        help='matagenomes group name')
    parser.add_argument('read_file', metavar='read_file', type=str, nargs='+',
                        help='read files')
    parser.add_argument('-n', '--n_threads', type=int,
                        help='number of bowtie/diamond threads (default: 20)', default=20)
    args = vars(parser.parse_args())
    print(args)