"""
This script holds a class that parses an hit list file created by mast.
"""
import sys, os
from os import listdir
from collections import defaultdict
import numpy
from Bio import SeqIO
from scipy.stats import rankdata

class MastHitListParser(object):
    """
    The Class parses a mast hitlist file to this example output:
    excpected_dict = {
        'micA' : [(-1, 4, 14, 1559.93, 2.02e-07)]
    }
    """
        
    def __init__(self, input_string, is_file = True):
        
        self.filename = None
        self.input_string = None
        if is_file:
            self.filename = input_string
        else:
            self.input_string = input_string
        self.elements = defaultdict(list)
        
    def ParseMastHitList(self):
        """
        The function parses a known srna file.
        """
        self.elements.clear()
        if self.filename:
            data_set = open(self.filename, 'r')
        else:
            data_set = self.input_string.strip().split('\n')
        for line in data_set:
            if not line.startswith('#'):
                splitted = line.split()
                srna = splitted[0]
                #TODO: revert changes.
                #if splitted[1][0] == '-':
                motif_index = splitted[1]
                #else:
                    #continue # We ignore all '+' strand!
                    #strand = 1
                #print line
                start = splitted[2]
                end = splitted[3]
                score1 = splitted[4]
                p_value = splitted[5]

                try:
                    self.elements[srna].append( ( int(motif_index), int(start), int(end), float(score1), float(p_value) ) )
                except ValueError:
                    print score1, p_value
        return self.elements
    
    def FilterByElement(self, directory, element, motif_dict):
        """
        The function returns the site/sites of an element by it's size.
        The function also returns the rank of the element, sorted by it's p_value.
        """
        elements_sites = []
        element_rank_list = []
        element_name_list = []
        for element_key in self.elements:
            for tup in self.elements[element_key]:
                element_rank_list.append(tup[4]) # list of p_value
                element_name_list.append(element_key+'_'+str(tup[0])+'_'+str(tup[1])+'_'+str(tup[2])) # list of element_key+from+to
        #print element_name_list
        ranked_list = rankdata(element_rank_list)
        
        for element_key in self.elements:
            if element_key.lower().startswith(element.lower()):
                #print element_rank_list
                #print ranked_list
                for x, tup in enumerate(self.elements[element_key]):
                    motif_dict[element].append( (directory, x+1, tup[1], tup[2], tup[4], ranked_list[element_name_list.index(element_key+'_'+str(tup[0])+'_'+str(tup[1])+'_'+str(tup[2]))]) )
        