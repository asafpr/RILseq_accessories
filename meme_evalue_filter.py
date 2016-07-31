import os, sys
import numpy
import argparse

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
        '-d', '--meme_parent_dir', default = "/home/hosts/disk16/niv/RILseq-meme-runs-new-fix",
        help='The meme parent directory to use')
    args = parser.parse_args(argv)
    return args


def main(argv=None, statsTupleList=[]):
    args = process_command_line(argv)

    known_srnas = ["dsrA","chiX","mcaS","rprA","gcvB","ryjA","micL","omrB","omrA","rydC","oxyS","micA","sgrS"
        ,"fnrS","rybB","ryhB","micF","cyaR","sdsR","glmZ","spf","arcZ","micC"]
    knownEvalue = []
    sigEvalue = {}
    directory_list = next(os.walk(args.meme_parent_dir))[1]
    for directory in directory_list:
        meme_filename = args.meme_parent_dir+'/'+directory+'/meme.txt'
        if os.path.isfile(meme_filename) and os.stat(meme_filename).st_size>0:
            with open(meme_filename) as meme_file:
                lines = meme_file.readlines()
                for line in lines:
                    if "E-value = " in line:
                        eval = float(line.split("E-value = ")[1].strip())
                        if directory.split('_')[0] in known_srnas and eval < 10:
                            knownEvalue.append(eval)
                            print "EVAL:\t"+directory.split('_')[0]+"\t" + directory + "\t" + str(eval)
                        if (eval <= 10) and ("IGR" in directory or "TU" in directory or "AS" in directory):
                            sigEvalue[directory] = eval
                    if line.startswith("data"):
                        n = line.strip().split()[-1]
                        print "INPUT: " + directory + " " + str(n)


    print numpy.mean(knownEvalue)
    for i in sigEvalue:
        print i


if __name__ == '__main__':
    status = main()
    sys.exit(status)
