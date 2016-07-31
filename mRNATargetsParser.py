"""
This script holds a class that parses a mRNA target file that holds binding sites for sRNA's.
"""
import sys, os
from os import listdir
from collections import defaultdict
import numpy
from Bio import SeqIO
import warnings

class mRNATargetsParser(object):
    """
    The Class parses a mRNA targets file to this example output:
    excpected_dict = {'arcZ' :
        [('eptB',"CCAGGG", -1, 3708513, 3708508, -10, -15),
        ('flhD',"CGACATCA", -1, 1976275, 1976282, -61, -54),
        ('rpoS',"GGGAAATC", -1, 2865670, 2865677, -104, -97),
        ('sdaC',"CCAGG", 1, 2926238, 2926242, -13, -9),
        ('tpx',"CCAGGG", -1, 1386810, 1386815, 20, 25)]
    }

    """
        
    def __init__(self, filename):
        
        self.filename = filename
        self.elements = defaultdict(list)
        
    def ParsemRNATargets(self):
        """
        The function parses a mRNA targets file.
        """
        self.elements.clear()
        for line in open(self.filename, 'r'):
            splitted = line.split()
            sequence = splitted[0]
            mrna = splitted[2]
            if splitted[5] == "reverse":
                strand = -1
            elif splitted[5] == "forward":
                strand = 1
            else:
                warnings.warn("There is no strand in the mRNA targets file")
            
            binding = line.split("binding: ")[1]
            splitted_binding = binding.split()
            start_bind = splitted_binding[0]
            end_bind = splitted_binding[1]
            
            bind = line.split("BIND ")[1]
            splitted_bind = bind.split()
            srna = splitted_bind[0]
            srna = srna[0].lower() + srna[1:]

            cur = line.split("cur: ")[1]
            splitted_cur = cur.split()
            start_cur = splitted_cur[0]
            end_cur = splitted_cur[1]
            self.elements[srna].append( (mrna, sequence, int(strand), int(start_bind), int(end_bind), int(start_cur), int(end_cur) ) )
                
        return self.elements
