#!/usr/bin/env python

"""
Summarizes all the data to generate table S4
"""

import sys, os
import argparse
import glob
from Bio import SeqIO
from os import path
import warnings
import imp
foo = imp.load_source('MastHitListParser', '/home/users/niv/chimeric_help/bin/MastHitListParser.py')
import MastHitListParser

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
        '-n', '--niv_res', default='/home/users/niv/chimeric_help/results/results-new-fix',
        help='Dir of Niv results.')
    parser.add_argument(
        '-m', '--meme_runs', default='/home/hosts/disk16/niv/RILseq-meme-runs-new-fix',
        help='The meme output dir.')
    parser.add_argument(
        '-t', '--tarsfile', default='/home/users/niv/chimeric_help/data/updated-temp-with-sdsr-muts-sep6-2015-combined-genomic-binding-coords-with-seq',
        help='targets file from Yael.')
    

    settings = parser.parse_args(argv)

    return settings

def count_targets(dirname, tarsfile, inputfile, total_input, mrna_res):
    """
    Count the number of targets of a sRNA, recalled by RILseq and found binding
    Arguments:
    - `dirname`: name of directory 
    - `tarsfile`: file of targets binding sites
    - `inputfile`: targets sequences
    - `mrna_res`: results of mRNA binding sites in MEME
    """
    tarslist = set()
    sdname = dirname.split('/')[-1]
    srna = sdname.split('_')[0]
    # List of targets
    with open(tarsfile) as tin:
        for line in tin:
            spl = line.strip().split()
            sn = spl[10]
            if sn[0].lower()+sn[1:] == srna:
                tarslist.add(spl[2])
    foundlist = set()
    # Targets in input file
    input_len = 0
    for sseq in SeqIO.parse(inputfile, 'fasta'):
        input_len += 1
        for tar in tarslist:
            overlap = sseq.id.split('_')[0].split(':')
            overlap_trim = [i for i in overlap if not "IGR" in i]
            if tar in ':'.join(overlap_trim):
                foundlist.add(tar)
    # Targets with MEME match
    true_bs = set()
    with open(mrna_res) as fin:
        for line in fin:
            if sdname in line:
                for tar in foundlist:
                    if tar in line and line.strip().endswith('True'):
                        true_bs.add(tar)
    totalinput_len = 0
    for sseq in SeqIO.parse(total_input, 'fasta'):
        totalinput_len += 1

    return max(totalinput_len, input_len), len(tarslist), len(foundlist), len(true_bs)

def match_known_bs(dirname, srfile):
    """
    Return True or False if there is a match between the motif and known BS
    Arguments:
    - `dirname`: The meme dir
    - `srfile`: detailed sRNA coverage
    """
    dirname = dirname.split('/')[-1]
    sname = dirname.split('_')[0]
    with open(srfile) as sin:
        for line in sin:
            if line.startswith(sname):
                cl = sin.next()
                bpos = set([i for i,d in enumerate(cl.strip().split()[-1]) if d != '0'])
                mc = sin.next()
                if not mc.startswith('m_cov'):
                    return False
                for l2 in sin:
                    if dirname in l2:
                        xpos = set([i for i, d in enumerate(l2.strip().split()[1]) if d != '-'])
                        return bool(xpos & bpos)
    return False

def read_meme(fname):
    """
    Read the meme output file and return the p-values and consensus sequences
    of each motif. turn T into U
    Arguments:
    - `fname`: meme.txt file
    """
    pvals = []
    conss = []
    seqs = []
    frac = []
    n = 0
    try:
        with open(fname) as fin:
            for line in fin:
                if line.startswith("data"):
                    n = line.strip().split()[-1]
                if line.startswith("MOTIF"):
                    pvals.append("%.3g"%float(line.strip().split()[-1]))
                if "regular expression" in line:
                    _ = fin.next()
                    conss.append(fin.next().strip().replace('T', 'U'))
                if "seqs=" in line:
                    seqs.append(int(line.split("seqs=")[1].strip()))

    except IOError:
        return 0, [], [], 0

    # Taking the fraction of sequences that entered the motif.
    for seq in seqs:
        if n > 0:
            frac.append('%.2f' % (seq/float(n)))
        else:
            frac.append('0')
    return n, pvals, conss, frac, seqs

def read_mast(fname):
    """
    Read the mast file and return the e-value and pattern
    Arguments:
    - `fname`: mast.txt file
    """
    eval = None
    patt = None
    try:
        with open(fname) as fin:
            for line in fin:
                if '+ strand' in line: # Skip the + strand
                    for _ in range(3):
                        _ = fin.next()
                if line.startswith("  LENGTH"):
                    eval = "%.3g"%float(line.strip().split()[-1])
                    patt = fin.next().strip().split()[-1]
                    break
    except IOError:
        return 'NA', 'NA'
    return eval, patt

def get_known_sRNAs(knownbs):
    """
    Return a list of sRNAs with known binding site
    Arguments:
    - `knownbs`: The targets binding sites file
    """
    srs = set()
    with open(knownbs) as sin:
        for line in sin:
            sn = line.split()[10]
            srs.add(sn[0].lower()+sn[1:])
    return srs
    

def main(argv=None):
    settings = process_command_line(argv)
    known_bs = get_known_sRNAs(settings.tarsfile)
    header = ['Gene name', 'Condition', 'RNA location in chimera', '# of targets in analyzed set',
              'MEME e-value', 'MEME motif', '# of targets with motif', 'Total # of targets',
              'Has known binding site?', '# of known targets with binding site',
              '# of known targets re-discovered by RIL-seq', 'Binding site recognized in targets', "MAST p-value",
              'Matches known binding site']
    print "\t".join(header)

    for dname in glob.glob("%s/*"%settings.meme_runs):
        outvec = []
        if not path.exists("%s/meme.txt"%dname):
            continue
        if len(dname.split('/')[-1].split('_')) == 2:
            outvec.extend(dname.split('/')[-1].split('_'))
            outvec.extend(["second"])

        elif len(dname.split('/')[-1].split('_')) == 3:
            outvec.extend(dname.split('/')[-1].split('_'))
        else:
            warnings.warn("Problem parsing the folder name: %s" % dname)
        n, meme_evs, meme_mots, meme_frac, meme_targ = read_meme("%s/meme.txt"%dname)
        if int(n) < 4:
            continue
        if meme_evs:
            outvec.extend([n, meme_evs[0], meme_mots[0], meme_targ[0]])
        else:
            outvec.extend([n, 'NA', 'NA', 'NA'])

        total_input = settings.meme_runs+'/'+dname.split('/')[-1].split('_')[0]+"_all_both/inputfile.fa"
        ti, kt, rt, bt = count_targets(
            dname, settings.tarsfile, "%s/inputfile.fa"%dname, total_input, "%s/detailed_mrna.txt"%settings.niv_res)
        outvec.extend([ti])

        outvec.append(dname.split('/')[-1].split('_')[0] in known_bs)

        outvec.extend([kt, rt, bt])
        #motif_list = []
        min_pval = "NMH"
        if os.path.isfile("%s/mast.hit_list_001"%dname) and os.stat("%s/mast.hit_list_001"%dname).st_size>0:
            mast_parser = MastHitListParser.MastHitListParser(open("%s/mast.hit_list_001"%dname).read(), is_file=False)# args.mast_output_dir+'/'+directory+'.mast.hit_list')
            mast_results = mast_parser.ParseMastHitList()

            mast_pvals = [i[4] for i in mast_results[dname.split('/')[-1].split('_')[0]] if i[0] == -1]  # only motif 1
            if mast_pvals:
                min_pval = str(min(mast_pvals))

        #mastev, _ = read_mast("%s/mastout_m1/mast.txt"%dname)
        outvec.append(min_pval)
        match = match_known_bs(dname, "%s/detailed_srna_coverage.txt"%settings.niv_res)
        outvec.append(match)
        print "\t".join([str(o) for o in outvec]) 
        
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
