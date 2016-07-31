import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
import argparse

def process_command_line(argv):
    """
    Return an args list
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='plot_srna_meme_figure main app')

    # define options here:
    parser.add_argument(
        '-f', '--srna_hits', default = "/home/users/niv/chimeric_help/results/results-new-fix-known",
        help='A known_srna file for parsing')
    parser.add_argument(
        '-f2', '--srna_hits_0001', default = "/home/users/niv/chimeric_help/results/results-new-fix-known-0.001",
        help='The mast p_value threshold')
    parser.add_argument(
        '-d', '--meme_parent_dir', default = "/home/hosts/disk16/niv/RILseq-meme-runs-new-fix",
        help='The meme parent directory to use')
    parser.add_argument(
        '-o', '--figure_path', default = "/home/users/niv/chimeric_help/results/figures/Fig4a.png",
        help='The figure path for saving the file')

    args = parser.parse_args(argv)

    return args

def main(argv=None, statsTupleList=[]):

    args = process_command_line(argv)
    srna_hits_table = args.srna_hits + "/table_srna.csv"
    srna_hits_table_0001 = args.srna_hits_0001 + "/table_srna.csv"
    homedir = args.meme_parent_dir

    #known_srnas = ["dsrA","chiX","mcaS","rprA","gcvB","ryjA","micL","omrB","omrA","rydC","oxyS","micA","sgrS"
    #    ,"fnrS","rybB","ryhB","micF","cyaR","sdsR","glmZ","spf","arcZ","micC","mgrR","dicF"]

    known_srnas = ['ArcZ', 'ChiX', 'CyaR', 'DsrA', 'FnrS', 'GcvB', 'McaS', 'MicA',
                   'OmrB', 'RprA', 'RybB', 'RydC', 'RyhB', 'SdsR', 'Spf', 'MgrR']

    manual_motifs = ["arcZ_all:1", "chiX_iron:1", "cyaR_all_both:1", "dsrA_all:2", "fnrS_stat:1", "gcvB_iron:1",
                     "mcaS_iron_both:1", "micA_stat_both:1", "omrB_all_first:1", "rprA_all_both:1",
                     "rybB_iron_both:1", "rydC_all_both:1", "ryhB_all_both:1", "sdsR_all:1", "spf_stat_both:1",
                     "mgrR_iron:1"]

    sequence_mot = []
    sequence_all = []
    srna_hits = []
    srna_hits_0001 = []

    directory_list = next(os.walk(homedir))[1]
    for motif in manual_motifs:
        seqs_mot = []
        pvals = []
        seqs_all = []
        for directory in directory_list:
            if motif.split(':')[0] == directory:
                meme_filename = homedir+'/'+directory+'/meme.txt'
                if os.path.isfile(meme_filename) and os.stat(meme_filename).st_size>0:
                    #print directory
                    with open(meme_filename) as meme_file:
                        lines = meme_file.readlines()
                        for line in lines:
                            if line.startswith("data"):
                                seqs_all.extend([line.strip().split()[-1], line.strip().split()[-1], line.strip().split()[-1]])
                            if line.startswith("MOTIF"):
                               pvals.append(float("%.3g"%float(line.strip().split()[-1])))
                            if "seqs" in line:
                                seqs_mot.append(int(line.split("seqs=")[1].strip()))

                if pvals:
                    motif_index = int(motif.split(':')[1])-1
                    sequence_mot.append(int(seqs_mot[motif_index]))
                    sequence_all.append(int(seqs_all[motif_index]))

                else:
                    sequence_mot.append(0)
                    sequence_all.append(0)

        directory_list = next(os.walk(homedir))[1]
        sequence_all_both = []
        for srna in known_srnas:
            srna = srna[0].lower()+srna[1:]
            all_both_seqs = 0
            for directory in directory_list:
                if srna+"_all_both" == directory:
                    input_filename = homedir+'/'+directory+'/inputfile.fa'
                    if os.path.isfile(input_filename) and os.stat(input_filename).st_size>0:
                        all_both_seqs = len(open(input_filename).readlines())/2
            sequence_all_both.append(int(all_both_seqs))

        for line in open(srna_hits_table, 'r'):
            if motif.split('_')[0] in line:
               srna_hits.append(line.strip().split(',')[-1])
        for line in open(srna_hits_table_0001, 'r'):
            if motif.split('_')[0] in line:
               srna_hits_0001.append(line.strip().split(',')[-1])

    if not os.path.exists('/'.join(args.figure_path.split('/')[:-1])):
        os.makedirs('/'.join(args.figure_path.split('/')[:-1]))

    plot_bar_chart(known_srnas, sequence_all, sequence_all_both, sequence_mot, srna_hits, srna_hits_0001,
                   args.figure_path)

def plot_bar_chart(known_srnas, sequence_all, sequence_all_both, sequence_mot, srna_hits, srna_hits_0001, path):

    # print len(sequence_all)
    # print len(sequence_mot)

    frac_list = []

    for x,i in enumerate(sequence_all):
        try:
            frac_list.append(sequence_mot[x]/float(sequence_all_both[x]))
        except ZeroDivisionError:
            frac_list.append(0)
    n_groups = len(known_srnas)

    print len(frac_list)

    fig, ax = plt.subplots()
    fig = plt.figure(figsize=(24.0, 12.0), dpi=100)
    fig.patch.set_visible(False)

    index = np.arange(n_groups)
    bar_width = 0.6
    opacity = 0.7

    rects1 = plt.bar(index, tuple(frac_list), bar_width,
                     alpha=opacity,
                     color='darkslategray')
                    #color= 'green')

    rcParams['font.family'] = 'sans-serif'
    # hfont_32 = {'fontname':'Helvetica', 'fontsize':'32'}
    # hfont_25_bold = {'fontname':'Helvetica', 'fontsize':'25', 'weight':'bold'}
    # hfont_25 = {'fontname':'Helvetica', 'fontsize':'25'}
    # hfont_28 = {'fontname':'Helvetica', 'fontsize':'38'}

    for x,rect in enumerate(rects1):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.01,
                '%d\n---\n%d'%(sequence_mot[x],sequence_all_both[x]), ha='center', va='bottom', weight='bold', fontsize='25')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    rects1[-1].set_alpha(0.35)

    plt.xlabel('Known sRNAs', fontsize='36')
    plt.ylabel('Fraction with motif', fontsize='36')
    # plt.title('Fraction of targets in MOTIF for known sRNAs', fontsize='24')
    plt.xticks(index + bar_width/2, [i for i in known_srnas], fontsize='28')
    plt.yticks([(i+1)*0.2 for i in range(5)], fontsize='32')
    plt.margins(0.01, 0)
    plt.tick_params(bottom='off', top='off', left='off', right='off')

    sig_srna_hits = [i for i in range(len(srna_hits)) if srna_hits[i]=='True']
    sig_srna_hits_0001 = [i for i in range(len(srna_hits_0001)) if srna_hits_0001[i]=='True']

    [i.set_color("purple") for x,i in enumerate(plt.gca().get_xticklabels()) if x in sig_srna_hits_0001]
    #[rects1[x].set_color("purple") for x,i in enumerate(plt.gca().get_xticklabels()) if x in sig_srna_hits_0001]
    [i.set_color("blue") for x,i in enumerate(plt.gca().get_xticklabels()) if x in sig_srna_hits]
    #[rects1[x].set_color("blue") for x,i in enumerate(plt.gca().get_xticklabels()) if x in sig_srna_hits]

    stringent = mpatches.Patch(color='blue', label='0.0001')
    medium = mpatches.Patch(color='purple', label='0.001')

    #plt.tight_layout()
    legend = plt.legend(frameon=False, title="Threshold", handles = [medium, stringent], loc='upper right', bbox_to_anchor=(1.05, 1), fontsize='28')
    legend.get_title().set_fontsize('32')
    legend.get_title().set_weight('bold')

    fig.savefig(path, dpi=fig.dpi, bbox_inches='tight')
    #plt.show()

if __name__ == '__main__':
    status = main()
    sys.exit(status)
