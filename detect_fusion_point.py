#!/usr/bin/env python

"""
Given already mapped fusions using the reads file (format:
gene1 gene2 position1 strand1 position2 strand2 read_name)
# Changed to format:
chr1 position1 strand1 chr2 position2 strand2 read_name
Use the original BAM file to locate the fusion on the sequence if possible.
Since the position of the first and second mates is known try to extend this
match on the 1st read, if the match doesn't reach the end of the read try
to find overlap between the left over and the second read and then extend
the second read using global alignment (actually glocal) including the region
that was aligned to the first position to see if there is some overlap, i.e.
nts that can map both of the positions. If the first read is fully mapped
try to map the second read and do the process again.

Report the positions in the genome that were fused and a PSSM of the nucleotides
before and after the fusion point.
"""

import sys
import optparse
import csv
from collections import defaultdict

from Bio.Seq import Seq
import pysam
from Bio import SeqIO
from RNA.structure import computeACC

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = optparse.OptionParser(
        formatter=optparse.TitledHelpFormatter(width=78),
        add_help_option=None)

    parser.add_option(
        '-l', '--list_reads',
        help='File with list of reads and their fused positions.')
    parser.add_option(
        '-b', '--bamfile',
        help='The original bam file with the full reads.')
    parser.add_option(
        '-g', '--genome',
        help='genome fasta file.')
    parser.add_option(
        '-p', '--printto',
        help='Output file name for the individual reads.')
    parser.add_option(
        '-o', '--overlap', type='int', default=4,
        help='Minimal overlap between mates to search for fusion point in.')
    parser.add_option(
        '-w', '--width', type='int', default=25,
        help='Width of PSSM to compute.')
    parser.add_option(      # customized description; put --help last
        '-h', '--help', action='help',
        help='Show this help message and exit.')

    settings, args = parser.parse_args(argv)

    # check number of arguments, verify values, etc.:
    if args:
        parser.error('program takes no command-line arguments; '
                     '"%s" ignored.' % (args,))

    # further process settings & args if necessary

    return settings, args

def get_reads_seqs(bamfile, rnames):
    """
    Return the sequences of all the reads from the bam file
    Arguments:
    - `bamfile`: The pysam file
    - `rnames`: reads names
    """
    r1_seqs = {}
    r2_seqs = {}
    rqns = set()
    reads = defaultdict(list)
    for read in bamfile.fetch(until_eof=True):
        rqns.add(read.qname)
        reads[read.qname].append(read)
    for rn in set(rnames) & rqns:
        for read in reads[rn]:
            if read.is_read1:
                outseq = Seq(read.seq)
                if not read.is_reverse:
                    outseq = outseq.reverse_complement()
                r1_seqs[read.qname] = str(outseq)
            else:
                outseq = Seq(read.seq)
                if read.is_reverse:
                    outseq = outseq.reverse_complement()
                r2_seqs[read.qname] = str(outseq)
    # r1_seqs is the 3' end of the second fused RNA, r2_seqs is the 5' of the
    # first fused RNA
    return r1_seqs, r2_seqs

def extend_alignment(rseq, pos5p, pos3p, is_read1, strand, genome):
    """
    Align the rseq to the genome in the specified position. Return the last
    position of the read mapped to the genome.
    Use local alignment
    Arguments:
    - `rseq`: Read sequence
    - `pos5p`: the 5' position, exact if read 2 or as limit if read 1
    - `pos3p`: the 3' position, exact if read 1 or as limit if read 2
    - `is_read1`: This read is read 1
    - `strand`: mapping strand
    - `genome`: The genome Seq object
    """
    rcnt = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    glen = len(genome)
    if is_read1:
        # Start from the last position and move on to the 5' end
        if strand == '-':
            ipos = 0
            try:
                while rcnt[genome[(pos3p+ipos)%glen]] == rseq[-(ipos+1)]:
                    ipos += 1
            except IndexError:
                return ipos - 1
            return ipos
        else:
            ipos = 0
            try:
                while genome[(pos3p-ipos)%glen] == rseq[-(ipos+1)]:
                    ipos += 1
            except IndexError:
                return ipos-1
            return ipos
    else:
        if strand == '-':
            ipos = 0
            try:
                while rcnt[genome[(pos5p-ipos)%glen]] == rseq[ipos]:
                    ipos += 1
            except IndexError:
                return ipos -1
            return ipos 
        else:
            ipos = 0
            try:
                while genome[(pos5p+ipos)%glen] == rseq[ipos]:
                    ipos += 1
            except IndexError:
                return ipos - 1
            return ipos 
        
                

        
def find_overlap(s1, s2):
    """
    Find overlaps between two reads. Assume they are both in the same
    orientation (r1 is revcomp)
    Return 3 seuqnces: s1, overlap, s2
    Arguments:
    - `s1`: first sequence, this is mate 2 actually in our experiments
    - `s2`: last sequence, mate 1
    """
    for i in range(min(len(s1), len(s2)))[::-1]:
        if s1[-i:]==s2[:i]:
            return s1[:-i], s1[-i:], s2[i:]
    return s1, '', s2


def main(argv=None):
    settings, args = process_command_line(argv)
    # Read the read names and positions
    read_5ps = {}
    read_3ps = {}
    read_genes = {}
    genome = SeqIO.read(settings.genome, 'fasta').seq
    for line in csv.reader(open(settings.list_reads), delimiter='\t'):
        if line[7] == 'single':
            continue
        read_5ps[line[6]] = [int(line[1])-1, line[2]]
        read_3ps[line[6]] = [int(line[4])-1, line[5]]
        read_genes[line[6]] = []#[line[1], line[4]]
    # Read the bam file and return the long sequences
    r1_seqs, r2_seqs = get_reads_seqs(
        pysam.Samfile(settings.bamfile), read_genes.keys())
    outer = csv.writer(open(settings.printto, 'w'), delimiter='\t')
    pssm_before = defaultdict(lambda: defaultdict(int))
    pssm_after = defaultdict(lambda: defaultdict(int))
    # For each read find the overlap, if exists and find the fusion point
    for rname in set(r1_seqs.keys())&set(r2_seqs.keys()):
        s1, overlap, s2 = find_overlap(r2_seqs[rname], r1_seqs[rname])
        if len(overlap) < settings.overlap:
            continue
#        print rname
#        print "%s %s %s"%(s1, overlap, s2)
        side_5p_len = extend_alignment(
            s1+overlap+s2, read_5ps[rname][0], 0, False, read_5ps[rname][1],genome)
        side_3p_len = extend_alignment(
            s1+overlap+s2, 0, read_3ps[rname][0], True, read_3ps[rname][1], genome)
#        print "%d %d %d %d %d"%(side_5p_len, side_3p_len, len(s1), len(overlap), len(s2))
        if side_5p_len + side_3p_len == len(s1) + len(overlap) + len(s2):
            # Report this as a fusion point
            f1_point = read_5ps[rname][0] + side_5p_len - 1
            f1_seq = genome[read_5ps[rname][0]:f1_point+1]
#            f1_seq = genome[f1_point-settings.width: f1_point+1]
            if read_5ps[rname][1] == '-':
                f1_point = read_5ps[rname][0] - side_5p_len + 1
                f1_seq = genome[f1_point:read_5ps[rname][0]+1].reverse_complement()
#                f1_seq = genome[f1_point:f1_point+settings.width+1].reverse_complement()
            f2_point = read_3ps[rname][0] - side_3p_len + 1
            f2_seq = genome[f2_point:read_3ps[rname][0]+1]
#            f2_seq = genome[f2_point:f2_point+settings.width+1]
            if read_3ps[rname][1] == '-':
                f2_point = read_3ps[rname][0] + side_3p_len - 1
                f2_seq = genome[read_3ps[rname][0]:f2_point+1].reverse_complement()
#                f2_seq = genome[f2_point-settings.width:f2_point+1].reverse_complement()
            acc1 = computeACC((s1+overlap)[:side_5p_len], winlen=1)
            acc2 = computeACC((overlap+s2)[-side_3p_len:], winlen=1)
            outer.writerow([rname] + read_genes[rname] + read_5ps[rname] + read_3ps[rname] + [f1_point, f1_seq, f2_point, f2_seq]+acc1[-settings.width:]+acc2[:settings.width])
            # Add the sequences to the pssm
            for i, nt in enumerate(f1_seq):
                pssm_before[len(f1_seq)-i][nt] += 1
            for i, nt in enumerate(f2_seq):
                pssm_after[i][nt] += 1
            print "%sNN%s"%(f1_seq, f2_seq)
#    print pssm_before
#    print pssm_after
    for nt in ('A', 'C', 'G', 'T'):
        for i in range(settings.width+1,0,-1):
            sys.stdout.write("%s\t"%pssm_before[i][nt])
        sys.stdout.write("\n")
    for nt in ('A', 'C', 'G', 'T'):
        for i in range(settings.width+1):
            sys.stdout.write("%s\t"%pssm_after[i][nt])
        sys.stdout.write("\n")
        
            
            
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
