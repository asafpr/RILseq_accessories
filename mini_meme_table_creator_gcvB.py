import sys, os
from collections import Counter
import argparse

def process_command_line(argv):
    """
    Return an args list
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='PolyASiteAnalyzer main app')

    # define options here:
    parser.add_argument(
        '-s', '--sig_file',
        help='A RIL-seq experiment file')
    parser.add_argument(
        '-o', '--output', default = "/home/hosts/disk16/niv/RILseq-files/meme_table_gcvB.txt",
        help='The mast output_dir to write to')

    args = parser.parse_args(argv)
    return args


def main(argv=None, statsTupleList=[]):
    args = process_command_line(argv)

    known_srnas = ["dsrA","chiX","mcaS","rprA","gcvB","ryjA","micL","omrB","omrA","rydC","oxyS","micA","sgrS"
            ,"fnrS","rybB","ryhB","micF","cyaR","sdsR","glmZ","spf","arcZ","micC", "mgrR", "dicF", "rseX"]

    with open(args.output, 'w') as meme_table:
        for known in known_srnas:
            meme_table.write(known+'\t'+args.sig_file+'\t'+known+'_pgcvb'+'\n')


if __name__ == '__main__':
    status = main()
    sys.exit(status)
