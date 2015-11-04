#!/usr/bin/env python

"""
Read a summary file made by pro_clash_significant_regions.py of pro_clash
package and for the given sRNA generate a fasta file with all interacting
regions
"""

import sys
import argparse
import csv
from Bio import SeqIO

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
        'summary_file',
        help='A pro_clash results file.')
    parser.add_argument(
        'sRNA',
        help='Name of sRNA that is found in the fourth column (common name).'
        'if "ALL_ROWS" is given return all sequences.' )
    parser.add_argument(
        'genome',
        help='A fasta file of the genome.')
    parser.add_argument(
        '--first', default=False, action='store_true',
        help='Assume the sRNA is in the first column, default is second.')
    parser.add_argument(
        '--both', default=False, action='store_true',
        help='Treat as sRNA in first and second columns.')
    parser.add_argument(
        '--pad', type=int, default=50,
        help='Pad the sequence before and after.')
    parser.add_argument(
        '--EC_chrlist', default='COLI-K12,chr',
        help='A comma separated dictionary of chromosome names from the '
        'EcoCyc name to the bam name. See the names of chromosomes in bam file'
        ' using samtools view -H foo.bam.')


    settings = parser.parse_args(argv)

    return settings

def main(argv=None):
    settings = process_command_line(argv)
    if settings.EC_chrlist:
        chr_dict = dict(zip(
                settings.EC_chrlist.split(',')[0::2],
                settings.EC_chrlist.split(',')[1::2]))
    else:
        chr_dict = {}
    genome = {}
    for seqrec in SeqIO.parse(open(settings.genome), 'fasta'):
        genome[seqrec.id] = seqrec.seq
    sRNA_cols = (2,)
    tar_cols = (1,)
    if settings.both:
        sRNA_cols = (1,2)
        tar_cols = (2,1)
    elif settings.first:
        sRNA_cols = (1,)
        tar_cols = (2,)
    written = set()
    for row in csv.DictReader(open(settings.summary_file), delimiter='\t'):
        for sRNA_col, tar_col in zip(sRNA_cols, tar_cols):
            if row['RNA%d name'%sRNA_col]==settings.sRNA or settings.sRNA=='ALL_ROWS':
                try:
                    chrname = chr_dict[row['RNA%d chromosome'%tar_col]]
                except KeyError:
                    chrname = row['RNA%d chromosome'%tar_col]
                tarname = row['RNA%d name'%tar_col]#.split(".")[0]
                if ("IGR" in tarname) or ("TU" in tarname):
                    glist = tarname.split(".")[:2]
                else:
                    glist = [tarname.split(".")[0]]
#                sys.stderr.write("%s\t%s\n"%(tarname,str(glist)))
                if set(glist) & written:
                    continue
                else:
                    written |= set(glist)
                try:
                    rseq = genome[chrname][max(int(row['RNA%d from'%tar_col])-1-settings.pad, 0):int(row['RNA%d to'%tar_col])+settings.pad]
                except KeyError:
                    continue
                if row['RNA%d strand'%tar_col] == '-':
                    rseq = rseq.reverse_complement()
    #            rseq = rseq.transcribe()
                print ">%s_%s_%s_%s\n%s"%(row['RNA%d name'%tar_col], row['RNA%d from'%tar_col], row['RNA%d to'%tar_col], row['RNA%d strand'%tar_col], rseq)
    # application code here, like:
    # run(settings, args)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
