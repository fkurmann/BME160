from sequenceAnalysis import FastAreader
from sequenceAnalysis import NucParams
'''This program takes properly formatted genomic information from either an input file or stdin
and returns datapoints on the analyzed input.'''

def main (fileName=None):
    '''Use FastAreader to parse an fa file of codons into an instance of the NucParams class.
    Sort the amino acids and codons in alphabetical order, then return formatted information
    about the sequence.'''
    # Create class instances
    myReader = FastAreader(fileName) 
    myNuc = NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)
        
    # Sort codons in alphabetic order, by Amino Acid, excluding codons that code to nothing
    codons = []
    for codon, aminoAcid in myNuc.rnaCodonTable.items():
      codons.append((codon, aminoAcid))
    codonsSorted = sorted(codons, key=lambda x: (x[1], x[0]))

    # Output part 1: Sequence length
    megaBases = myNuc.nucCount() / 1000000
    print('{0:0.2f} Mb'.format(megaBases))
    print('')
    # Output part 2: GC content
    gcContent = (((myNuc.nucComp['G'] + myNuc.nucComp['C']) * 100) / (sum(myNuc.nucComp.values()) - myNuc.nucComp['N']))
    print('{0:0.1f}%'.format(gcContent))
    print('')

    # Output part 3: Relative codon usage per codon
    for tuple in codonsSorted:
      print ('{:s} : {:s} {:5.1f} ({:6d})'.format(tuple[0], tuple[1], (myNuc.codonComp[tuple[0]] / myNuc.aaComp[tuple[1]]) * 100, myNuc.codonComp[tuple[0]]))

if __name__ == "__main__":
    main()