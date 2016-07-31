#!/usr/bin/env python

"""
Read a summary file made by RILseq_significant_regions.py of RILseq
package and for the given sRNA generate a fasta file with all interacting
regions
"""

import sys
import argparse
import csv
import collections
from Bio import SeqIO
import warnings

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
        help='A RILseq results file.')
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
    plus_dict = {}
    minus_dict = {}
    for row in csv.DictReader(open(settings.summary_file), delimiter='\t'):
        for sRNA_col, tar_col in zip(sRNA_cols, tar_cols):
            if row['RNA%d name'%sRNA_col]==settings.sRNA or settings.sRNA=='ALL_ROWS':
                try:
                    chrname = chr_dict[row['RNA%d chromosome'%tar_col]]
                except KeyError:
                    chrname = row['RNA%d chromosome'%tar_col]
                ## commented in order to fix a problem regarding the target picking.
                ## joins the different elements from that overlap so we won't miss targets.
                tarname = row['RNA%d name'%tar_col]#.split(".")[0]
#                 if ("IGR" in tarname) or ("TU" in tarname):
#                     glist = tarname.split(".")[:2]
#                 else:
#                     glist = [tarname.split(".")[0]]
#                 #sys.stderr.write("%s\t%s\n"%(tarname,str(glist)))
#                 if set(glist) & written:
#                     continue
#                 else:
#                     written |= set(glist)

                ## build a genome hash table with elements coordinates.

                if row['RNA%d strand'%tar_col] == '+':
                    try:
                        plus_dict[(max(int(row['RNA%d from'%tar_col]), 0),int(row['RNA%d to'%tar_col]))] = tarname
                    except KeyError:
                        continue
                elif row['RNA%d strand'%tar_col] == '-':
                    try:
                        minus_dict[(max(int(row['RNA%d from'%tar_col]), 0),int(row['RNA%d to'%tar_col]))] = tarname
                    except KeyError:
                        continue
                else:
                    warnings.warn("Strand has to be + or -")

                #print ">%s_%s_%s_%s\n%s"%(row['RNA%d name'%tar_col], row['RNA%d from'%tar_col], row['RNA%d to'%tar_col], row['RNA%d strand'%tar_col], rseq)
    ## joining the elements that overlap in their coordinates
    ordered_plus_dict = collections.OrderedDict(sorted(plus_dict.items()))
    ordered_minus_dict = collections.OrderedDict(sorted(minus_dict.items()))

    plus_elements_dictionary = {}
    minus_elements_dictionary = {}

    #sys.stderr.write("len ordered_plus_dict: %s\n"%(str(len(ordered_plus_dict))))
    #sys.stderr.write("ordered_plus_dict: %s\n"%(str(ordered_plus_dict)))
    join_overlapping_elements(ordered_plus_dict, plus_elements_dictionary)
    #sys.stderr.write("len plus_elements_dictionary: %s\n" % (str(len(plus_elements_dictionary))))

    join_overlapping_elements(ordered_minus_dict, minus_elements_dictionary)

    for key, value in plus_elements_dictionary.iteritems():
        try:
            rseq = genome[chrname][max(key[0]-1-settings.pad,0):key[1]+settings.pad]
            print ">%s_%s_%s_%s\n%s" % (value, max(key[0]-settings.pad,0), key[1]+settings.pad, '+', rseq)
        except KeyError:
            pass

    for key, value in minus_elements_dictionary.iteritems():
        try:
            rseq = genome[chrname][key[0]-1-settings.pad:key[1]+settings.pad]
            rseq = rseq.reverse_complement()
            print ">%s_%s_%s_%s\n%s" % (value, max(key[0]-settings.pad,0), key[1]+settings.pad, '-', rseq)
        except KeyError:
            pass

    # application code here, like:
    # run(settings, args)
    return 0        # success

"""
The function receives a dictionary and fills is up with the overlapping elements from a different dictionary.
input:
{(20,200): "tar1", (100,250):"tar2"}
output:
{(20,250):"tar1-tar2"}
"""
def join_overlapping_elements(elements_dict, overlapping_dict):
    sum_overlaps = 0
    for key, value in elements_dict.iteritems():
            to_del = []
            overlap = False
            for overlap_key in overlapping_dict:
                #sys.stderr.write("key: %s\n" % (str(key)))
                #sys.stderr.write("overlap_key: %s\n" % (str(overlap_key)))
                if not (key[0] > overlap_key[1] or key[1] < overlap_key[0]):
                    overlap = True
                    sum_overlaps += 1
                    coords = overlap_key
                    name = overlapping_dict[overlap_key]
                    to_del.append(overlap_key)
                    key = (min(coords[0], key[0]), max(coords[1], key[1]))
                    value = name+":"+value
                    if abs(key[1]-key[0]) >= 400:
                        sys.stderr.write("The overlapping element exceeds 400 bases. num: %s, name: %s\n" %
                                         (str(abs(key[1]-key[0])), value))
                    continue

            if overlap:
                for delete_key in to_del:
                    del overlapping_dict[delete_key]
                overlapping_dict[key] = value
            else:
                overlapping_dict[key] = value

if __name__ == '__main__':
    status = main()
    sys.exit(status)
