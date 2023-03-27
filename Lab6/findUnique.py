#!/usr/bin/env python3
# Name: Fabrice Kurmann (fkurmann)
# Group Members: James Arnold, Ronnie Voskoboynik, Tucker Morrison, Yahaira Uribe

########################################################################
# Design 
#
# Main:
# Combine head/sequences into list, sort by head
# 
# For sequence, read into UniqueFinder object, which stores all sequence data
# UniqueFinder class deals with process of finding unique and essential sets described below
# 
# For sorted head/sequences items
#   Print head and sequence, then unique sequence items in order of position
#     Determine position by searching for each subsequence, storing position together with sequence in tuple
#
# UniqueFinder:
#
# Initizlie with empty powerset, uniqe, and essential lists
# Each list will store one set corresponding to data for one sequence each
# Call powerset function, unique function, essential function
#
# Powerset function:
#   For each sequence in sequences
#     Traverse sequence with nested for loops finding each subsequence, add all to set, no duplicates
#
# Unique function:
#   For each set in powersets
#     Subtract the union of all other sets to end with subsequences only in this powerset
#
# Essential function:
#   For each set in unique
#     For each item in set
#       Compare item with all other items, if item inside another item, delete the other, it is non essential
#
# Now that all above functions have run, shrinking list of sets, return remaining set to main.
#
########################################################################

########################################################################
# FastA reader 
########################################################################
import sys
import re
class FastAreader:
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

########################################################################
# Unique finder 
########################################################################
class UniqueFinder:
  '''Class to find, unique subsequences of genes sequences that have been read in.'''

  def __init__ (self, seqs):
    '''Constructor for new object, usually one object created per sequence.'''
    self.sequences = seqs
    self.powersets = []
    self.uniques = []
    self.essentials = []

    # Call functions to create sets and then filter them
    self.powerset()
    self.unique()
    self.essential()

  def powerset (self):
    '''Parse sequence to get all subsequences and arrange in a set to form a powerset.'''
    
    # Add a set of subsequences for each sequence
    for sequence in self.sequences:
      currentPowerset = set()

      # Remove whitespace in input, capitalize
      sequence = sequence[1]
      sequence = sequence.upper()
      sequence = sequence.strip()

      # Loop over all subsets of the sequence, add subsets that haven't already been added
      for i in range(0, len(sequence)):
        for j in range(i, len(sequence) + 1):
          subsequence = sequence[i:j]
          if (len(subsequence) == 0):
             continue
          if not (subsequence in currentPowerset):
            currentPowerset.add(subsequence)
      
      self.powersets.append(currentPowerset)

  def unique(self):
    '''Filter the unique subsequences from the set of all subsequeces for every other gene sequence'''

    # For every sequence's set 
    for powerset in self.powersets:
      otherUnion = set()
      # Sum up all sets of different sequences
      for otherset in self.powersets:
        if (powerset != otherset):
          otherUnion = otherUnion.union(otherset)
      
      # Find the items that only belong in this set
      self.uniques.append(powerset - otherUnion)

  def essential(self):
    '''Filter the essential subsequences from the unique subsequences.'''

    for uniqueSet in self.uniques:
      nonEssential = set()
      # Compare items of the same sequences' unique set to see which items are minimal
      for setMember in uniqueSet:
        for otherMember in uniqueSet:
          # Remove non minimal items
          if ((setMember != otherMember) and (setMember in otherMember)):
            nonEssential.add(otherMember)

      self.essentials.append(uniqueSet - nonEssential)

########################################################################
# Main
########################################################################
def main(inCL=None):
  '''Main method that takes input from stdin, creates uniqueFinder objects, 
  and returns results to stdout.'''
  
  # List for head/sequence tuples
  hs = []

  myReader = FastAreader(inCL) 

  # Read in sequences, remove unwanted characters
  for head, seq in myReader.readFasta():
    head = head.replace(' ', '')
    seq = seq.replace('.', '')
    seq = seq.replace('_', '')
    seq = seq.replace('-', '')

    # Combine header and sequence lists so they are easier to sort
    hs.append((head, seq))

  # Sort header lines alphabetically and sequences in line with headers
  hs = sorted(hs, key=lambda x: x[0])

  # Make the uniqueFinder instance to complete the set operations
  unique = UniqueFinder(hs)
  unique = unique.essentials

  # Print the outputs as format suggests: head, sequence, essential subsets
  for index in range(0, len(hs)):
    print(hs[index][0])
    print(hs[index][1])
    uniquePrintList = []

    # Find all positions where the desired substring appears
    for item in unique[index]:
      positions = []
      start_index = 0
      # print(start_index)
      while (start_index < len(hs[index][1])):
        position = (hs[index][1]).find(item, start_index)
        if position == -1:
          break
        else:
          positions.append(position)
          start_index = position + 1

      # positions = [i.start() for i in re.finditer(re.escape(item), (hs[index][1]))]
      for position in positions:
        uniquePrintList.append((position, item))
      
    # Sort essential subsets by position before printing
    uniquePrintList = sorted(uniquePrintList, key=lambda x: x[0])
    for item in uniquePrintList:
      print(('.' * item[0]) + item[1])
 
if __name__ == "__main__":
  main() 