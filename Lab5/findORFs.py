#!/usr/bin/env python3
# Name: Fabrice Kurmann (fkurmann)
# Group Members: James Arnold, Ronnie Voskoboynik, Tucker Morrison, Yahaira Uribe

########################################################################
# Design 
#
# Pseudocode: Main
# create command line instance
# create OFR finder instance based on args
# create fasta instance
# for sequence in fasta
#   findGenes(sequence)
#   format print results
#
# Pseudocode: ORF Finder
#
# for left to right, right to left
#   for entire sequence, starting at positions 0, 1, 2
#     while stop not encountered
#       fill start list
#     fill gene list, empty start list
#       
#     add all dangling (no stop codon for gene) to gene list
#
# sort gene list by length, gene start position
#   
# for gene in gene list
#   report coding frame, start position, stop position, length
#
#
########################################################################

from sequenceAnalysis import FastAreader
from sequenceAnalysis import OrfFinder

'''This program takes a fasta input file given by STDIN and uses the ORF Finder class to return an output file with
correctly formatted information about the genes in the sequences in the input file.'''

########################################################################
# CommandLine
########################################################################
class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=True, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main 
########################################################################
def main(inFile = None, options = None):
    '''Use FastAreader to parse an fa file of codons into an instance of the NucParams class.
    Sort the amino acids and codons in alphabetical order, then return formatted information
    about the sequence.'''

    # Create class instances
    thisCommandLine = CommandLine(options)
    
    # Extra credit command line options
    minLength = thisCommandLine.args.minGene
    biggestOnly = thisCommandLine.args.longestGene
    startList = thisCommandLine.args.start
    # If startList was changed, remove defaults
    if len(startList) != 1:
        startList = startList[1:]
    stopList = thisCommandLine.args.stop
    # If stopList was changed, remove defaults
    if len(stopList) != 3:
        stopList = stopList[3:]
    
    myReader = FastAreader(inFile) 
    myORF = OrfFinder(minLength, biggestOnly, startList, stopList)

    # Read in sequences
    for head, seq in myReader.readFasta():
        foundORFs = myORF.findGenes(seq)

        # Print header line, then all genes found
        print (head)
        for gene in foundORFs:
            print ('{:+d} {:>5d}..{:>5d} {:>5d}'.format(gene[1], gene[2], gene[3], gene[4])) 
    
if __name__ == "__main__":
    main()
