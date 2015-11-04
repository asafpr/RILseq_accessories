#!/usr/bin/env python
from subprocess import call, Popen, PIPE
import fileinput
import os
import os.path
import shutil
from Bio import SeqIO

def run_line(line, genome, process, meme, seqer, mast, bnumer, cog_find):
    """
    Run one line of the table
    """
    srname, sfile, oroutdir = line.strip().split()[:3]
#    outdir = "%s_%s_meme_run"%(srname, sfile.split("/")[-1].split("_")[0])
    for outdir in (oroutdir, oroutdir+"_first", oroutdir+"_both"):
        try:
            os.makedirs(outdir)
        except OSError:
            pass
        with open("%s/inputfile.fa"%outdir, 'w') as outf:
            runvec = [process, sfile, srname, genome]
            if "first" in outdir:
                runvec = [process, '--first', sfile, srname, genome]
            elif "both" in outdir:
                runvec = [process, '--both', sfile, srname, genome]
            call(runvec, stdout=outf)
        # Run COG analyzer
        ids = []
        for sr in SeqIO.parse("%s/inputfile.fa"%outdir, 'fasta'):
            ids.append(sr.id.split('_')[0].split('.')[0])
        print ids
        with open("%s/cog_analysis_results.txt"%outdir, 'w') as cout:
            bber = Popen(bnumer, stdin=PIPE, stdout=PIPE)
            cfr = Popen(cog_find, stdin=bber.stdout, stdout=cout)
            for idd in ids:
                bber.stdin.write("%s\n"%idd)
            bber.stdin.close()
        # Run meme
        call(meme.split() + ['-oc', outdir, "%s/inputfile.fa"%outdir])
        if os.path.exists("manual_sRNA_seqs/%s_sRNA_seq.fa"%srname):
            shutil.copyfile(
                "manual_sRNA_seqs/%s_sRNA_seq.fa"%srname,
                "%s/srna_seq.fa"%outdir)
        else:
            with open("%s/srna_seq.fa"%outdir, 'w')as sf:
                sf.write(">%s\n"%srname)
                sf.flush()
                sq = Popen(seqer, stdout=sf, stdin=PIPE)
                sqline = srname
                if '3UTR' in srname:
                    sqline += " -50 100"
                if '5UTR' in srname:
                    sqline += " 50 -100"
                sq.communicate(input=sqline)
        
        call(mast.split() + ['-m', '1', '-oc', "%s/mastout_m1"%outdir, "%s/meme.txt"%outdir, "%s/srna_seq.fa"%outdir])
        call(mast.split() + ['-m', '2', '-oc', "%s/mastout_m2"%outdir, "%s/meme.txt"%outdir, "%s/srna_seq.fa"%outdir])
        call(mast.split() + ['-m', '3', '-oc', "%s/mastout_m3"%outdir, "%s/meme.txt"%outdir, "%s/srna_seq.fa"%outdir])

meme = "meme -dna -mod zoops -maxsize 300000 -minw 6 -maxw 15 -nmotifs 3"
mast = "mast -mt 0.01 -sep"
seqer = "/home/users/assafp/lib/sequence.pl"
cog_find = "/home/users/assafp/lib/cog_finder.pl"
bnumer = "/home/users/assafp/lib/gene2bnum.pl"
process = "/home/hosts/disk19/E_coli/Sahar/bin/prepare_meme_from_summary.py"
genome = "/home/hosts/disk19/reference_genomes/E_coli/genome.fa"
for line in fileinput.input():
    run_line(line, genome, process, meme, seqer, mast, bnumer, cog_find)


