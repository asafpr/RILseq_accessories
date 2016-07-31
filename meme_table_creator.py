import sys, os
from collections import Counter
import argparse

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
        '-l', '--sigfile_log', default='home/users/yael/sahar-clash/paper-anal-jul-2015/'
                                       'paper-final-revision-anal-feb-2016/selected-tables-for-final-anal/'
                                       'with-count-log-signif-libs-code_assign-type-to-signif-chimeras-of-flipped-'
                                       'MG_hfq-FLAG101-104_108_109_unified_Log_bwa.bam_all_fragments_l25.txt_sig_'
                                       'interactions_with-total.txt.with-known',
        help='Significant interactions fle for log condition')
    parser.add_argument(
        '-s', '-sigfile_stat', default='home/users/yael/sahar-clash/paper-anal-jul-2015/'
                                       'paper-final-revision-anal-feb-2016/selected-tables-for-final-anal/'
                                       'with-count-stat-signif-libs-code_assign-type-to-signif-chimeras-of-flipped-'
                                       'MG_hfq-FLAG209_210_312_bwa.bam_unified_Stationary_all_fragments_l25.txt_sig_'
                                       'interactions_with-total.txt.with-known',
        help='Significant interactions fle for stat condition')
    parser.add_argument(
        '-i', '--sigfile_iron', default='home/users/yael/sahar-clash/paper-anal-jul-2015/'
                                        'paper-final-revision-anal-feb-2016/selected-tables-for-final-anal/'
                                        'with-count-iron-signif-libs-code_assign-type-to-signif-chimeras-of-flipped-'
                                        'MG_hfq-FLAG207_208_305_unified_Iron_bwa.bam_all_fragments_l25.txt_sig_'
                                        'interactions_with-total.txt.with-known',
        help='Significant interactions fle for iron condition')
    parser.add_argument(
        '-a', '--sigfile_all', default="/home/hosts/disk16/niv/RILseq-files/unified_all.new",
        help='Significant interactions fle for all condition')

    parser.add_argument(
        '-t', '--table_name', default="/home/hosts/disk16/niv/RILseq-files/"
                                 "meme_table_final_unique_at_least_3_with_all.txt.new",
        help='Significant interactions fle for all condition')

    settings = parser.parse_args(argv)
    return settings


def main(argv=None):

    settings = process_command_line(argv)

    sig_filenames = [settings.sigfile_log, settings.sigfile_iron, settings.sigfile_stat, settings.sigfile_all]
    sig_conditions = ['log', 'iron', 'stat', 'all']
    element_list = []

    with open(settings.table_name, "w") as meme_table:
        for i, sig_filename in enumerate(sig_filenames):
            sig_file = open(sig_filename)
            for line in sig_file.readlines()[1:]:
                element_list.append(str(line.strip().split('\t')[4]))
                element_list.append(str(line.strip().split('\t')[5]))

            counted_dict = Counter(element_list)
            for key in counted_dict:
                if counted_dict[key]>=3:
                    meme_table.write("%s\t%s\t%s\n" % (key,sig_folder+'/'+sig_filename, key+'_'+sig_conditions[i]))
            element_list = []
            sig_file.close()

if __name__ == '__main__':
    status = main()
    sys.exit(status)

