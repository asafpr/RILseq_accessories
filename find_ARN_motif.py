#!/usr/bin/env python

"""
Find ARN(x) motif in sequence and compare to shuffled seqs
"""

import sys
import argparse
from collections import defaultdict
import re
from Bio import SeqIO
import random

from RNA.dinuc_shuffle import shuffle_difreq

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object, replace the description
    parser = argparse.ArgumentParser(
        description='Find ARN motif.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'fastafile',
        help='Input fasta file.')
    parser.add_argument(
        '-s', '--shuffles', type=int, default=999,
        help='Number of times to shuffle each sequence. Significance will be'
        ' measured according to the number of ARN repeats in all sequences.')
    parser.add_argument(
        '-r', '--region', default=False, action='store_true',
        help='Choose random region instead of shuffling.')
    parser.add_argument(
        '-g', '--genome', default=None,
        help='If random is chosen supply a genome.')
    settings = parser.parse_args(argv)

    return settings

def count_ARNs(seqs):
    """
    Count the number of times the ARN appears in the sequences. Return a
    dictionary from number of repeats to counts
    Arguments:
    - `seqs`: A dictionary of sequences
    """
    counts = defaultdict(int)
    motif = '((A[A][ACGT])+)'
    patt = re.compile(motif)

    for s in seqs.values():
        lengths = set([0])
        for mr in patt.findall(str(s)):
            lengths.add(len(mr[0]))
        for l in range(max(lengths)/3):
            counts[l+1] += 1
    return counts



def main(argv=None):
    settings = process_command_line(argv)
    if settings.region:
        genome = SeqIO.read(settings.genome, 'fasta').seq
    fas_in = {}
    for sr in SeqIO.parse(settings.fastafile, 'fasta'):
        fas_in[sr.id] = str(sr.seq)
    counts_table = count_ARNs(fas_in)
    shf_counts = []
    for sh in range(settings.shuffles):
        shf_seqs = {}
        for n, s in fas_in.items():
            if settings.region:
                rlen = len(s)
                reg_from = random.randint(0, len(genome)-rlen)
                shf_seqs[n] = genome[reg_from:reg_from+rlen]
                if random.random() > 0.5:
                    shf_seqs[n] = shf_seqs[n].reverse_complement()
            else:
                shf_seqs[n] = shuffle_difreq(s)
        shf_counts.append(count_ARNs(shf_seqs))
    # Compute significance for each number of counts
    for cn, counts in counts_table.items():
        above_count = len([sc[cn] for sc in shf_counts if sc[cn] >= counts])
        print "%s: %s times, p-value: %g"%(cn, counts, (above_count+1)/float(len(shf_counts) + 1))
        
    # application code here, like:
    # run(settings, args)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
