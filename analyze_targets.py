import sys, os
import argparse
import warnings
import re

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
        '-g', '--gff', default='/home/users/niv/chimeric_help/data/refseq_ver2_genes_chr.gff',
        help='genome gff file.')
    parser.add_argument(
        '-t', '--tarsfile', default='/home/users/niv/chimeric_help/results/results-new-fix-known/detailed_mrna_all.txt',
        help='targets file')

    settings = parser.parse_args(argv)

    return settings

def main(argv=None):
    settings = process_command_line(argv)
    header = ['sRNA', 'Target', 'Position relative to AUG', 'Absolute Position start',
              'Width', 'Strand', 'Sequence', 'p_value', 'sRNA_from', 'sRNA_to', 'sRNA_strand', 'sRNA_length']
    print "\t".join(header)

    targets_file = open(settings.tarsfile, 'r')
    for line in targets_file.readlines():
        spl = line.strip().split('\t')
        outvec = []
        outvec.append(spl[0].split('_')[0]) # sRNA
        outvec.append(spl[1]) # Target

        pra = combine_genome_and_tars(settings.gff, line)

        outvec.append(pra) # Position relative to AUG
        outvec.append(spl[2]) # Absolute Start
        outvec.append(abs(int(spl[3])-int(spl[2]))+1) # Width
        outvec.append(spl[4]) # Strand
        outvec.append(spl[5]) # Sequence
        outvec.append(spl[6]) # p_value
        outvec.append(spl[7]) # sRNA_from
        outvec.append(spl[8]) # sRNA_to
        outvec.append(spl[9]) # sRNA_strand
        outvec.append(spl[10]) # sRNA_length

        print "\t".join([str(o) for o in outvec])

    targets_file.close()

    return 0        # success

def combine_genome_and_tars(gff_filename, tar_line):
    """
    The function returns the position relative to the AUG of the target based on a gff file.
    :param gff_filename:
    :param tar_line:
    :return: position relative to the AUG of the target based on a gff file.
    """
    with open(gff_filename, 'r') as gff_file:
        pra_list = []
        for gff_line in gff_file.readlines():
            gff_start = int(gff_line.strip().split()[3])
            gff_end = int(gff_line.strip().split()[4])
            gff_name = gff_line.strip().split()[9][1:-2]
            tar_names = re.split('\.|\:',tar_line.strip().split('\t')[1])
            tar_start = int(tar_line.strip().split('\t')[2])
            tar_end = int(tar_line.strip().split('\t')[3])
            tar_strand = tar_line.strip().split('\t')[4]

            if gff_name in tar_names:  # and gff_end >= tar_end and gff_start <= tar_start:
                #  There is a match
                if tar_strand == '+':
                    pra_list.append(tar_start - gff_start)
                elif tar_strand == '-':
                    pra_list.append(gff_end - tar_end)
                else:
                    warnings.warn("Strand problem on Target line")
                    pra_list.append('NA')

        if len(pra_list) > 1:
            warnings.warn("More than 1 item in AUG list")
        elif len(pra_list) == 0:
            return 'NA'
        return pra_list[0]

def strand_processing(strand):
    """
    The function gets a strand string or integer and returns
    the strand in a integer 1/-1 format
    """
    if not (strand=="+" or strand=="-" or strand==1 or strand==-1):
        warnings.warn("Strand must be 1 or -1!")
    try:
        return int(strand+"1")
    except TypeError:
        return strand
if __name__ == '__main__':
    status = main()
    sys.exit(status)
