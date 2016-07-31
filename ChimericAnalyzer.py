"""
This script will analyze a RILseq run.
The script will get an element file from the experiment and a meme output parent directory,
and will return for each element, the coverage of the binding site according to the known targets DB.
It will also return the binding site from the experiment for graphical impression.
"""
import os, sys
import argparse
import KnownSrnaParser
import MastHitListParser
import mRNATargetsParser
import MemeFileParser
from collections import defaultdict
import subprocess
import csv
import warnings
import re
from Bio import SeqIO

def process_command_line(argv):
    """
    Return an args list
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]
        
    # initialize the parser object:
    parser = argparse.ArgumentParser(description='PolyASiteAnalyzer main app')

    # define options here:
    parser.add_argument(
        '-k', '--known_srna_file', default="/home/users/niv/chimeric_help/data/updated-temp-with-sdsr-muts-sep6-2015-srna-binding-coord-for-debug",
        help='A known_srna file for parsing')
    parser.add_argument(
        '-p', '--p_value_mast_threshold', default=0.0001, type=float,
        help='The mast p_value threshold')
    parser.add_argument(
        '--mev', '--e_value_mast_threshold',  default=10, type=float,
        help='The meme E-value threshold')
    parser.add_argument(
        '-e', '--mast_exe', default="mast",
        help='The mast executable to run')
    parser.add_argument(
        '-m', '--mrna_targets_file', default="/home/users/niv/chimeric_help/data/updated-temp-with-sdsr-muts-sep6-2015-combined-genomic-binding-coords-with-seq",
        help='The mast output_dir to write to')
    parser.add_argument(
        # '-d', '--meme_parent_dir', default = "/home/hosts/disk19/E_coli/Sahar/CLIP_Lig/CLASH_4th_run/cluster_results/meme_runs",
        '-d', '--meme_parent_dir', default="/home/hosts/disk16/niv/RILseq-meme-runs-new-fix-known",
        help='The meme parent directory to use')
    parser.add_argument(
        '-g', '--gff_file', default='/home/users/niv/chimeric_help/data/refseq_ver2_genes_chr.gff',
        help='genome gff file.')
    parser.add_argument(
        '--output_dir', default="../results/testdir/",
        help='Output csv filename for writing statistics')
    parser.add_argument(
        '--accessory_info', action='store_true', default=False,
        help='Whether to add accessory information')
    
    args = parser.parse_args(argv)
    
    return args

def main(argv=None, statsTupleList=[]):
    args = process_command_line(argv)

    try:
        os.makedirs(args.output_dir)
    except OSError:
        pass
    
    motif_dict = defaultdict(list)
    # known_srnas = ["dsrA","chiX","mcaS","rprA","gcvB","ryjA","micL","omrB","omrA","rydC","oxyS","micA","sgrS"
    #     ,"fnrS","rybB","ryhB","micF","cyaR","sdsR","glmZ","spf","arcZ","micC"]

    manual_motifs = ["arcZ_all:1", "chiX_iron:1", "cyaR_all_both:1", "dsrA_all:2", "fnrS_stat:1", "gcvB_iron:1",
                 "mcaS_iron_both:1", "micA_stat_both:1", "micF_all:1", "omrB_all_first:1", "rprA_all_both:1",
                 "rybB_iron_both:1", "rydC_all_both:1", "ryhB_all_both:1", "sdsR_all:1", "spf_stat_both:1",
                 "mgrR_iron:1", "sucD.EST3UTR_all_both:1", "uhpT.EST3UTR_all_both:1", "glnA.3UTR_iron_both:1",
                 "ykgH.EST3UTR_iron_both:1", "cpxP.3UTR_all_both:1", "fliC.fliA.IGR_iron_both:1",
                 "aceK_stat_both:2", "ymfH_all_both:1", "raiA_log_both:1", "gadE.mdtE.TU_stat_both:1"]
    maunal_motifs_input = [335,43,328,13,72,157,36,96,16,86,68,14,5,192,92,71,130,37,115,11,19,248,22,15,6,28,26]

    motif_hits_S5 = ["gcvB_iron:1", "chiX_iron:1", "fnrS_stat:1", "ryhB_all_both:1",
                     "micA_stat_both:1", "gadE.mdtE.TU_stat_both:1", "sucD.EST3UTR_all_both:1"]
    
    srna_parser = KnownSrnaParser.KnownSrnaParser(args.known_srna_file)
    known_srna_elements = srna_parser.ParseKnownSrna()
    known_srna_coverage_dict = srna_parser.GetCoverageDict()
    targetParser = mRNATargetsParser.mRNATargetsParser(args.mrna_targets_file)
    targets = targetParser.ParsemRNATargets()

    srna_covfile = open(args.output_dir+"/detailed_srna_coverage.txt", 'wb')
    all_srna_covfile = open(args.output_dir+"/detailed_all_srna_coverage.txt", 'wb')
    tablefile = open(args.output_dir+"/table_all_srna.csv", 'wb')
    csvwriter_table = csv.writer(tablefile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    csvfile = open(args.output_dir+"/table_srna.csv", 'wb')
    csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    mrna_covfile = open(args.output_dir+"/detailed_mrna.txt", 'wb')
    pra_file = open(args.output_dir+"/detailed_pra.txt", 'wb')
    figs_file = open(args.output_dir+"/data_for_figures_4bcd_S5.txt", 'wb')
    tableS5_file = open(args.output_dir+"/data_for_table_S5.txt", 'wb')

    directory_list = next(os.walk(args.meme_parent_dir))[1]
    for directory in directory_list:
        meme_filename = args.meme_parent_dir+'/'+directory+'/meme.txt'
        seq_filename = args.meme_parent_dir+'/'+directory+'/srna_seq.fa'
        input_filename = args.meme_parent_dir+'/'+directory+'/inputfile.fa'
        output = None
        if os.path.isfile(meme_filename) and os.path.isfile(input_filename) and \
                        os.path.isfile(seq_filename) and os.stat(meme_filename).st_size > 0:
            memeParser = MemeFileParser.MemeFileParser(meme_filename)
            inputParser = MemeFileParser.InputFileParser(input_filename)

            # Detailed mRNA analysis:
            # TODO: This function uses MOTIF 1 only.
            create_detailed_mRNA(mrna_covfile, targets, directory, memeParser, inputParser)

            if args.accessory_info:
                output = runMast(args.mast_exe, 0.01, args.mev,
                             args.meme_parent_dir+'/'+directory+'/'+'mast.hit_list_0.01', meme_filename, seq_filename)
            else:
                output = runMast(args.mast_exe, args.p_value_mast_threshold, args.mev,
                             args.meme_parent_dir+'/'+directory+'/'+'mast.hit_list_'+str(args.p_value_mast_threshold), meme_filename, seq_filename)
        
            if output:
                mast_parser = MastHitListParser.MastHitListParser(output, is_file=False)#args.mast_output_dir+'/'+directory+'.mast.hit_list')
                mast_elements = mast_parser.ParseMastHitList()
                mast_parser.FilterByElement(directory, directory.split('_')[0], motif_dict)

                # Detailed Coverage Motifs for all seqs
                create_detailed_seqs(all_srna_covfile, csvwriter_table, mast_elements, seq_filename)

                if mast_elements and args.accessory_info:
                    # Accessory files
                    create_accessory_files(pra_file, figs_file, tableS5_file, directory, seq_filename,
                                           memeParser, inputParser, mast_elements, manual_motifs, maunal_motifs_input,
                                           motif_hits_S5, args.meme_parent_dir, args.gff_file)

    # Detailed Coverage Motifs
        create_detailed_sRNA(srna_covfile, known_srna_coverage_dict, motif_dict)

    # Table sRNA analysis:
    create_table(csvwriter, known_srna_elements, known_srna_coverage_dict, motif_dict)

def runMast(mast_exe, p_value_mast_threshold, e_value_mast_threshold, output_filename, meme_txt, srna_fasta):
    """
    Runs the mast according to the params given.
    """
    # if not os.path.exists('/'.join(output_file.split('/')[:-1])):
    #     os.makedirs('/'.join(output_file.split('/')[:-1]))
    try:    
        output = subprocess.check_output([mast_exe, meme_txt, srna_fasta, '-mt',  str(p_value_mast_threshold), '-mev', str(e_value_mast_threshold), '-hit_list', '-sep'])
        open(output_filename, 'w').write(output)
    except subprocess.CalledProcessError as exc:
        return None
    return output


def create_detailed_sRNA(srna_covfile, known_srna_coverage_dict, motif_dict):
    for known_key in known_srna_coverage_dict:
        print_string = ""
        mot_print_string = ""
        srna_covfile.write( known_key[0] + '\t' + known_key[1] + '\n')
        for i in known_srna_coverage_dict[known_key]:
            if i >= 10:
                print_string += 'i'
            else:
                print_string += str(i)
        srna_covfile.write( 'cov\t' + print_string+'\n' )
        if motif_dict.has_key(known_key[0]):
            mot_coverage = {}
            mot_coverage[known_key[0]] = [0 for i in range(len(known_key[1]))] # creating a new zero list in the length of the srna.
            for motif_site_tup in motif_dict[known_key[0]]: #running through motif_sites tuples
                for i in range(motif_site_tup[2]-1, motif_site_tup[3]):
                    mot_coverage[known_key[0]][i] += 1
            for i in mot_coverage[known_key[0]]:
                if i >= 10:
                    mot_print_string += 'i'
                else:
                    mot_print_string += str(i)
            srna_covfile.write( 'm_cov:\t'+mot_print_string + '\n')
            for x, motif_site_tup in enumerate(motif_dict[known_key[0]]): #running through motif_sites tuples
                srna_covfile.write('mot\t'+'-'*(int(motif_site_tup[2])-1)+'x'*(int(motif_site_tup[3])-int(motif_site_tup[2])+1)+'-'*(len(known_key[1])-int(motif_site_tup[3])) + '\t' + motif_site_tup[0]+'_'+str(motif_site_tup[1])+'_'+str(motif_site_tup[4]) +' rank: '+ str(motif_site_tup[5])+"\n")
        srna_covfile.write('\n')


def create_detailed_seqs(all_srna_covfile, csvwriter_table, mast_elements, seq_filename):
    print_string = ""
    seq_coverage = {}
    fasta_sequences = SeqIO.parse(open(seq_filename), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
    for key in mast_elements:
        is_coverage = False
        seq_coverage[key] = [0 for i in range(len(sequence))]  # creating a new zero list in the length of the srna.
        for mast_tuple in mast_elements[key]:
            for i in range(mast_tuple[1]-1, mast_tuple[2]):
                seq_coverage[key][i] += 1
        all_srna_covfile.write("{0: <30}\t{1}\n".format(name, sequence))
        for i in seq_coverage[key]:
            if i > 0:
                is_coverage = True
            if i >= 10:
                print_string += 'i'
            else:
                print_string += str(i)
        all_srna_covfile.write('{0: <30}\t{1}\n'.format("cov:", print_string))
        all_srna_covfile.write('\n')
        csvwriter_table.writerow([name, is_coverage])


def create_table(csvwriter, known_srna_elements, known_srna_coverage_dict, motif_dict):
    for known_key in known_srna_elements:
        is_motif_in_binding_site = "False"
        num_targets = len(known_srna_elements[known_key])
        # TODO: ask yael, 2 different binding sites in the same mRNA-target is considered 2 or 1?
        for motif_site_tup in motif_dict[known_key[0]]:  # running through motif_sites tuples
            hit = 0
            for i in range(motif_site_tup[2]-1, motif_site_tup[3]):
                if known_srna_coverage_dict[known_key][i] > 0:
                    hit += 1
            if (motif_site_tup[3]-motif_site_tup[2]+1)/float(2) <= hit:
                is_motif_in_binding_site = "True"

        csvwriter.writerow([known_key[0], num_targets, is_motif_in_binding_site])


def create_detailed_mRNA(mrna_covfile, targets, directory, memeParser, inputParser):
    for key in targets:
        for target_tup in targets[key]:
            if directory.split('_')[0] == key:  # only directories with the correct srna name.
                for input_key, input_val in inputParser.ParseInputFile().iteritems():
                    input_string = input_key[0]+'_'+str(input_key[1])+'_'+str(input_key[2])+'_'+str(input_key[3])
                    meme_dict, meme_list, meme_widths, meme_pvals, meme_regex  = memeParser.ParseMemeFile()
                    for meme_key in meme_list[0]:
                        is_motif_in_binding_site_mrna = "False"
                        mot_sequence = meme_list[0][meme_key][3]+meme_list[0][meme_key][4]
                        mot_index = meme_list[0][meme_key][0]-1
                        if meme_key.split('*')[0] in input_string and target_tup[0] in input_key[0]\
                                and input_val.find(mot_sequence) == mot_index:
                            # The motif in the meme.txt file contains the input_file string
                            # and the mRNA is in the targets file and the sequences match
                            if input_key[3] == -1:
                                motif_coords = (input_key[2]-meme_list[0][meme_key][0]+1-meme_widths[0]+1,input_key[2]-meme_list[0][meme_key][0]+1)
                            elif input_key[3] == 1:
                                motif_coords = (input_key[1]+meme_list[0][meme_key][0]-1,input_key[1]+meme_list[0][meme_key][0]-1+meme_widths[0]-1)
                            else:
                                warnings.warn("Not a correct strand!")
                            hit = 0
                            for i in range(motif_coords[0], motif_coords[1]+1):
                                if i >= target_tup[3] and i <= target_tup[4] > 0:
                                    hit += 1
                            if (target_tup[4]-target_tup[3]+1)/float(3) <= hit:
                                is_motif_in_binding_site_mrna = "True"
                            mrna_covfile.write( "HIT: exp_condition: {0}, target_element and start: {1}, {2}, calculated_coords: {3}, target_coords: {4}, is_binding: {5}\n".format(
                                directory, input_key[0]+'_'+str(input_key[1])+'_'+str(input_key[2])+'_'+str(input_key[3]), meme_list[0][meme_key][0], motif_coords, (target_tup[3],target_tup[4]) , is_motif_in_binding_site_mrna) )


def create_accessory_files(pra_file, figs_file, tableS5_file, directory, seq_filename, memeParser, inputParser,
                           mast_elements, manual_motifs, maunal_motifs_input, motif_hits_S5, meme_parent_dir, gff_file):

    for x, motif in enumerate(motif_hits_S5):
        motif_name = motif.split(':')[0]
        motif_index = int(motif.split(':')[1])-1
        if directory == motif_name:  # only directories with the correct srna name.
            srna_len = len("".join(open(seq_filename).readlines()[1:]).replace('\n','').replace(' ',''))
            for input_key, input_val in inputParser.ParseInputFile().iteritems():
                input_string = input_key[0]+'_'+str(input_key[1])+'_'+str(input_key[2])+'_'+str(input_key[3])
                meme_dict, meme_list, meme_widths, meme_evals, meme_regex = memeParser.ParseMemeFile()
                mast_pvals = [i[4] for i in mast_elements[motif_name.strip().split('_')[0]]]
                mast_index = mast_pvals.index(min(mast_pvals))
                mast_from = str(mast_elements[motif_name.strip().split('_')[0]][mast_index][1])
                mast_to = str(mast_elements[motif_name.strip().split('_')[0]][mast_index][2])
                mast_strand = "NA"
                if int(mast_elements[motif_name.strip().split('_')[0]][mast_index][0]) > 0:
                    mast_strand = '+'
                elif int(mast_elements[motif_name.strip().split('_')[0]][mast_index][0]) < 0:
                    mast_strand = '-'
                for meme_key in meme_list[motif_index]:
                    mot_sequence = meme_list[motif_index][meme_key][3]+meme_list[motif_index][meme_key][4]
                    mot_index = meme_list[motif_index][meme_key][0]-1
                    if meme_key.split('*')[0] in input_string and input_val.find(mot_sequence) == mot_index:
                        #  The motif in the meme.txt file contains the input_file string
                        #  and the mRNA is in the targets file.
                        if input_key[3] == -1:
                            motif_coords = (input_key[2]-meme_list[motif_index][meme_key][0]+1-meme_widths[motif_index]+1,input_key[2]-meme_list[motif_index][meme_key][0]+1)
                            input_strand = '-'
                        elif input_key[3] == 1:
                            motif_coords = (input_key[1]+meme_list[motif_index][meme_key][0]-1,input_key[1]+meme_list[motif_index][meme_key][0]-1+meme_widths[motif_index]-1)
                            input_strand = '+'
                        else:
                            warnings.warn("Not a correct strand!")
                            motif_coords = (0, 0)
                            input_strand = 'NA'
                        tableS5_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (directory, input_key[0],
                                                                             str(motif_coords[0]),
                                                                             str(motif_coords[1]),
                                                                             input_strand,
                                                                             str(meme_list[motif_index][meme_key][1]),
                                                                             meme_list[motif_index][meme_key][3],
                                                                             mast_from,
                                                                             mast_to,
                                                                             mast_strand,
                                                                             srna_len))
    for x, motif in enumerate(manual_motifs):
        motif_name = motif.split(':')[0]
        motif_index = int(motif.split(':')[1])-1
        if directory == motif_name:  # only directories with the correct srna name.
            srna_len = len("".join(open(seq_filename).readlines()[1:]).replace('\n','').replace(' ',''))

            seq_len, mot_pos_list, mlens, evals, conss = memeParser.ParseMemeFile()
            # output = runMast(args.mast_exe, 0.01, args.mev, args.meme_parent_dir+'/'+directory+'/'+'mast.hit_list', meme_filename, seq_filename)
            mast_pvals = [i[4] for i in mast_elements[motif_name.strip().split('_')[0]]]
            mast_index = mast_pvals.index(min(mast_pvals))
            mast_from = str(mast_elements[motif_name.strip().split('_')[0]][mast_index][1])
            mast_to = str(mast_elements[motif_name.strip().split('_')[0]][mast_index][2])
            mast_strand = "NA"
            if int(mast_elements[motif_name.strip().split('_')[0]][mast_index][0]) > 0:
                mast_strand = '+'
            elif int(mast_elements[motif_name.strip().split('_')[0]][mast_index][0]) < 0:
                mast_strand = '-'

            figs_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (motif_name, len(mot_pos_list[motif_index]),
                                                          str(maunal_motifs_input[x]),
                                                          evals[motif_index],
                                                          str(mast_elements[motif_name.strip().split('_')[0]][mast_index][4]),
                                                          meme_parent_dir+directory+"/logo"+str(motif_index+1)+".png",
                                                            mast_from,
                                                            mast_to,
                                                            mast_strand,
                                                            meme_parent_dir+directory+"/srna_seq.fa"))

            for input_key, input_val in inputParser.ParseInputFile().iteritems():
                input_string = input_key[0]+'_'+str(input_key[1])+'_'+str(input_key[2])+'_'+str(input_key[3])
                meme_dict, meme_list, meme_widths, meme_evals, meme_regex = memeParser.ParseMemeFile()
                for meme_key in meme_list[motif_index]:
                    mot_sequence = meme_list[motif_index][meme_key][3]+meme_list[motif_index][meme_key][4]
                    mot_index = meme_list[motif_index][meme_key][0]-1
                    if meme_key.split('*')[0] in input_string and input_val.find(mot_sequence) == mot_index:
                        #  The motif in the meme.txt file contains the input_file string
                        #  and the mRNA is in the targets file.
                        if input_key[3] == -1:
                            motif_coords = (input_key[2]-meme_list[motif_index][meme_key][0]+1-meme_widths[motif_index]+1,input_key[2]-meme_list[motif_index][meme_key][0]+1)
                            input_strand = '-'
                        elif input_key[3] == 1:
                            motif_coords = (input_key[1]+meme_list[motif_index][meme_key][0]-1,input_key[1]+meme_list[motif_index][meme_key][0]-1+meme_widths[motif_index]-1)
                            input_strand = '+'
                        else:
                            warnings.warn("Not a correct strand!")
                            motif_coords = (0, 0)
                            input_strand = 'NA'
                        pra_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                                (directory.split('_')[0],
                                                 input_key[0],
                                                 combine_genome_and_tars(gff_file,
                                                                         input_key[0], motif_coords[0], motif_coords[1],
                                                                         input_strand),
                                                 str(motif_coords[0]),
                                                 abs(int(motif_coords[1])-int(motif_coords[0]))+1,
                                                 input_strand,
                                                 str(meme_list[motif_index][meme_key][1]),
                                                 meme_list[motif_index][meme_key][3],
                                                 mast_from,
                                                 mast_to,
                                                 mast_strand,
                                                 srna_len))


def combine_genome_and_tars(gff_filename, tar_name, tar_start, tar_end, tar_strand):
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
            tar_names = re.split('\.|\:', tar_name) #tar_line.strip().split('\t')[1])
            # tar_start = int(tar_line.strip().split('\t')[2])
            # tar_end = int(tar_line.strip().split('\t')[3])
            # tar_strand = tar_line.strip().split('\t')[4]

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

if __name__ == '__main__':
    status = main()
    sys.exit(status)
