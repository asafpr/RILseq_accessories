"""
This script holds a class that parses a meme file holding motifs.
"""
import sys, os
from os import listdir
from collections import defaultdict
import numpy
from Bio import SeqIO
from copy import deepcopy
import warnings

class InputFileParser(object):
    """
    The Class parses a file to this example output:
    excpected_dict = {
    ('agrA', 3646086, 3646169, 1) : "TTAGATGATGGCTATCTCACTCCAGTCAGAGCCACCAACTCAGGGCTGGAAAGTAAAAAACCGACGCAAAGTCGGTTTTTTTAC",
    ('arcZ', 3348599, 3348719, 1) : "GTGCGGCCTGAAAAACAGTGCTGTGCCCTTGTAACTCATCATAATAATTTACGGCGCAGCCAAGATTTCCCTGGTGTTGGCGCAGTATTCGCGCACCCCGGTCTAGCCGGGGTCATTTTTT",
    ('arrS', 3656009, 3656077, -1) : "GTAATCCGATTTAAATATCGAGTCTCCTTGTTTCGACTTAAGCTGGCAATTGGATTGCCAGCTTTCTTT",
    ('chiX', 506428, 506511, 1) : "ACACCGTCGCTTAAAGTGACGGCATAATAATAAAAAAATGAAATTCCTCTTTGACGGGCCAATAGCGATATTGGCCATTTTTTT"
    }
    """
        
    def __init__(self, filename):   
        self.filename = filename
        self.elements = {}

    def ParseInputFile(self):
        """
        The function parses a fasta input file.
        """
        fasta_sequences = SeqIO.parse(open(self.filename),'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            name_list = name.strip().split('_')
            try:                
                #name_tuple = ( name_list[0], int(name_list[1])-50, int(name_list[2])+50, self.strand_processing(name_list[3]) )
                ##changed to take the correct coordinates from the input_file.fa
                name_tuple = ( name_list[0], int(name_list[1]), int(name_list[2]), self.strand_processing(name_list[3]) )
            except ValueError:
                name_tuple = ( name_list[0]+'_'+name_list[1], int(name_list[2]), int(name_list[3]), self.strand_processing(name_list[4]) )
                #warnings.warn("Problem parsing the inputfile.fa. File: "+ self.filename + " Name: " + name)
            self.elements[name_tuple] = sequence
        return self.elements    

    def strand_processing(self, strand):
        """
        The function gets a strand string or integer and returns
        the strand in a integer 1/-1 format
        """
        if not (strand=="+" or strand=="-" or strand==1 or strand==-1):
            #warnings.warn("Strand must be 1 or -1!")
            raise ValueError("strand processing")
        try:
            return int(strand+"1")
        except TypeError:
            return strand
        
class MemeFileParser(object):
    """
    The Class parses a meme file to this example output:
    """
        
    def __init__(self, filename):
        
        self.filename = filename
        self.elements = defaultdict(int)
        
    def ParseMemeFile(self):
        """
        Read the meme file and return a dictionary with sequence name and length,
        sequence name and position of motif and motif length.
        Works for motif 1 only
        Arguments:
        - `memefile`: meme.txt file
        """
        seq_len = {}
        mot_pos_list = []
        mot_pvals_list = []
        evals = []
        mlens = []
        conss = []
        memein = open(self.filename)
        line = memein.next()
        while not line.startswith('Sequence name'):
            line = memein.next()
        _ = memein.next()  # Read the ----- line
        for line in memein:
            if line.startswith('*'):
                break
            spl = line.strip().split()
            seq_len[spl[0]] = int(spl[2])
            try:
                seq_len[spl[3]] = int(spl[5])
            except IndexError:
                pass  # short line, only 1 sequence

        for line in memein:
            self.elements.clear()

            # Read until the MOTIF line
            if line.startswith("MOTIF"):
                evals.append(float("%.3g" % float(line.strip().split()[-1])))
                mlens.append(int(line.strip().split()[5]))
            if "Sequence name" in line and "Start" in line:
                _ = memein.next()
                table_line = memein.next()
                while not table_line.startswith('-'):
                    spl = table_line.strip().split()
                    target_key = spl[0]+'*'+spl[1]
                    if target_key not in self.elements:
                        try:
                            self.elements[target_key] = [int(spl[1]), float("%.3g" % (float(spl[2]))),
                                                     spl[3], spl[4], spl[5]]
                        except IndexError:
                            self.elements[target_key] = [int(spl[1]), float("%.3g" % (float(spl[2]))),
                                                     spl[3], spl[4], ""]
                    else:  #the same name already exists
                        warnings.warn("There are multiple targets with the same name: %s" % (spl[0]))
                        print "%s: There are multiple targets with the same name: %s" % (self.filename, spl[0])
                    table_line = memein.next()
                #  print self.elements
                mot_pos_list.append(deepcopy(self.elements))
            if "regular expression" in line:
                _ = memein.next()
                conss.append(memein.next().strip().replace('T', 'U'))
        
        return seq_len, mot_pos_list, mlens, evals, conss
