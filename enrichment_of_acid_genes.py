#!/usr/bin/env python

"""
Test if targets of a sRNA are enriched with a list of genes
"""
from scipy.stats import fisher_exact
import sys
import argparse
import csv
from collections import defaultdict

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object, replace the description
    parser = argparse.ArgumentParser(
        description='Script description here.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-l', '--geneslist', default='/home/users/yael/sahar-clash/acid-resistance/symbol-of-acid-related-genes.ids-with-letter-code',
        help='A table with the acid response genes, gene names in last column.')
    parser.add_argument(
        'results_table',
        help='The results of RILseq.')
    parser.add_argument(
        '-m', '--mingroup', type=int, default=5,
        help='Minimal number of genes in intersect.')
    parser.add_argument(
        '-n', '--ngenes', type=int, default=4147,
        help='Number of genes in the genome.')

    settings = parser.parse_args(argv)

    return settings

def test_fisher(list_tars, list_genes, n):
    """
    Run Fisher's exact test to see if the list_tars and list_genes have larger
    than expected of intersect
    Arguments:
    - `list_tars`: The list of targets
    - `list_genes`: The list of acid genes
    - `n`: Total number of genes in the genome
    """
    ints = len(set(list_tars) & set(list_genes))
    return ints, fisher_exact(
        [[ints, len(list_tars)-ints], [len(list_genes)-ints,
                                       n-len(list_tars)-len(list_genes)+ints]],
        alternative='greater')[1]


def main(argv=None):
    settings = process_command_line(argv)
    with open(settings.geneslist) as sin:
        acid_genes = [l.strip().split()[-1] for l in sin.readlines()]
    sRNA_cols = (1,2)
    tar_cols = (2,1)

    stars = defaultdict(set)
    for row in csv.DictReader(open(settings.results_table), delimiter='\t'):
        for sRNA_col, tar_col in zip(sRNA_cols, tar_cols):
            tarname = row['RNA%d name'%tar_col]
            if ("IGR" in tarname) or ("TU" in tarname):
                glist = tarname.split(".")[:2]
            else:
                glist = [tarname.split(".")[0]]
            stars[row['RNA%d name'%sRNA_col]] |= set(glist)
    for k, v in stars.items():
        if len(v) < settings.mingroup:
            continue
        ints, pv = test_fisher(v, acid_genes, settings.ngenes)
        print "%s (%d): %d, %g"%(k, len(v), ints, pv)
        print v & set(acid_genes)
    # application code here, like:
    # run(settings, args)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
