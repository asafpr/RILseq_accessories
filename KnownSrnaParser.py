"""
This script holds a class that parses a known srna db containing sRNAs and targets.
"""
import sys, os
from os import listdir
from collections import defaultdict
import numpy
from Bio import SeqIO

class KnownSrnaParser(object):
    """
    The Class parses a Known Srna file to this example output:
    excpected_dict = {
    ('arcZ', "GTGCGGCCTGAAAAACAGTGCTGTGCCCTTGTAACTCATCATAATAATTTACGGCGCAGCCAAGATTTCCCTGGTGTTGGCGCAGTATTCGCGCACCCCGGTCTAGCCGGGGTCATTTTTT") :
    [('eptB', 69, 74), ('flhD', 72, 79), ('rpoS', 64, 71), ('sdaC', 70, 74), ('tpx', 69,74)]
    }
    """
        
    def __init__(self, filename):
        
        self.filename = filename
        self.elements = defaultdict(list)
        self.coverage = {}
        
    def ParseKnownSrna(self):
        """
        The function parses a known srna file.
        """
        self.elements.clear()
        
        for line in open(self.filename, 'r'):
            splitted = line.split('\t')
            srna = splitted[1][0].lower() + splitted[1][1:]
            sequence = splitted[0].split(' ')[0]
            target = splitted[2]
            target_coords = splitted[3].strip().split(' ')
            for coords in target_coords:
                self.elements[(srna,sequence)].append((target, int(coords.split('-')[0]), int(coords.split('-')[1])))
                
        return self.elements
    
    def GetCoverageDict(self):
        """
        The function uses the parsed known srna dictionary and calculates the coverage of the binding sites.
        The function returns a dictionary in this format:
        {
            (<srna_name>, <srna_seq>) : [list_of_coverage_per_position]
            .
            .
        }
        """
        self.coverage.clear()
        for key in self.elements:
            self.coverage[key] = [0 for i in range(len(key[1]))] # creating a new zero list in the length of the srna.
            for target, start, end in self.elements[key]:
                for i in range(start-1, end):
                    self.coverage[key][i] += 1
        return self.coverage

