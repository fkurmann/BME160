########################################################################
# ORF finder 
########################################################################
class OrfFinder:
    '''Class to find, sort, and return ORFs of length 100 or longer in a FASTA file'''

    def __init__ (self, minLength, biggestOnly, startList, stopList):
        '''Constructor for new object, usually one object created per input instance.'''

        # Dictionary for base pairings
        self.basePairings = {'A':'T',
                            'T':'A',
                            'G':'C',
                            'C':'G'}
        
        # Store command line specified settings
        self.startCodons = startList
        self.minLength = minLength
        self.biggestOnly = biggestOnly
        self.stopCodons = stopList

    def findGenes (self, sequence):
        '''Find all ORFs in a sequence, then index coding frames for all 3 possible positions,
        reverse and flip the sequence and repeat the process for the opposie direction'''
        
        # Remove whitespace in input, capitalize
        sequence = sequence.upper()
        sequence = sequence.strip()

        # Build opposing sequence by reversing sequence order and finding corresponding bases
        opposingSequence = ''
        for i in range(len(sequence) - 1, -1, -1):
            opposingSequence += self.basePairings[sequence[i]]
        # print('Sequence and Opposing sequence:')
        # print(opposingSequence)
        # print(sequence)
        sequences = [sequence, opposingSequence]

        # Properties of each verified ORF reside in tuples
        # (sequence, coding frame, start pos, stop pos, len)

        # Properties of each found "frame ORF" reside in tuples
        # (start pos, stop pos)
        
        ORFs = []
        opposingTrack = -1
        # Iterate over both directions and sides of the sequence
        for sequence in sequences:
            opposingTrack += 1
            # Iterate over all three coding frames
            for j in range(3):
                # List of all start positions for the currently read frame, 1 always included too
                startPositions = [1]
                # List of all ORFs in this particular frame
                frameORFs = []
                # Iterate over codons in sequence
                for i in range(0, len(sequence), 3):
                    codon = sequence[i: i + 3]

                    # Check for start codon
                    if (codon in self.startCodons):                        
                        if ((i + 1 + j) not in startPositions):
                            startPositions.append(i + 1 + j)

                    # Check for stop codon
                    if (codon in self.stopCodons):
                        for start in startPositions:
                            # Add all ORFs to the list by index, start to stop
                            frameORFs.append((start, i + 3 + j))
                            # If only getting largest sequence, break after earliest start
                            if (self.biggestOnly == True):
                                break
                        # Clear start
                        startPositions = []

                # Add any dangling genes
                for start in startPositions:
                    # Add all ORFs to the list by index, start to stop + 1
                    frameORFs.append((start, len(sequence)))
                    # If only getting largest sequence, break after earliest start
                    if (self.biggestOnly == True):
                        break

                # Consider only ORFs or a minimum length
                for ORF in frameORFs:
                    if (ORF[1] - ORF[0] >= self.minLength):
                        # Add ORFs in upper sequence
                        if (opposingTrack == 0):
                            ORFs.append((sequence[ORF[0] - 1:ORF[1]], j+1, ORF[0], ORF[1], (ORF[1] - ORF[0] + 1)))
                        if (opposingTrack == 1):
                            if ((ORF[1] - ORF[0] + 1) == len(sequence)):
                                ORFs.append((sequence[ORF[0] - 1:ORF[1]], (-1)*(j+1), ORF[0], ORF[1], (ORF[1] - ORF[0] + 1)))
                            else:
                                ORFs.append((sequence[ORF[0] - 1:ORF[1]], (-1)*(j+1), len(sequence) - ORF[1] + 1, len(sequence) - ORF[0] + 1, (ORF[1] - ORF[0] + 1)))
                
                # Move onto next frame
                sequence += sequence[0]
                sequence = sequence[1:]
        
        # Sorts the genes found in a sequence based on length(desc), then start position(asc)
        ORFs = sorted(ORFs, key=lambda x: ((-1* x[4]), (x[2])))
        
        return (ORFs)

########################################################################
# Sequence feature analysis
########################################################################
class ProteinParam:
    '''Class to get properties of sequences of amino acids'''
  # These tables are for calculating:
  #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
  #     absorbance at 280 nm (aa2abs280)
  #     pKa of positively charged Amino Acids (aa2chargePos)
  #     pKa of negatively charged Amino acids (aa2chargeNeg)
  #     and the constants aaNterm and aaCterm for pKa of the respective termini
  #  Feel free to move these to appropriate methods as you like

  # As written, these are accessed as class attributes, for example:
  # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        '''Capitalize input sequence, then find values of all 20 amino acids in this sequence, fill in composition table'''
        self.sequence = protein.upper()
        self.composition = {
        'A': self.sequence.count('A'),  'G': self.sequence.count('G'),  'M': self.sequence.count('M'), 'S': self.sequence.count('S'), 'C': self.sequence.count('C'),
        'H': self.sequence.count('H'), 'N': self.sequence.count('N'), 'T': self.sequence.count('T'), 'D': self.sequence.count('D'), 'I': self.sequence.count('I'),
        'P': self.sequence.count('P'), 'V': self.sequence.count('V'), 'E': self.sequence.count('E'), 'K': self.sequence.count('K'), 'Q': self.sequence.count('Q'),
        'W': self.sequence.count('W'),  'F': self.sequence.count('F'), 'L': self.sequence.count('L'), 'R': self.sequence.count('R'), 'Y': self.sequence.count('Y')
        }

    def aaCount (self):
        '''Traverse over composition dictionary and add up all values'''
        total = 0
        for value in self.composition.values():
            total += value
        return total

    def pI (self):
        '''Loop over all possible pH values and check their resulting charge, whichever pH results in the most neuteral is returned'''
        # Store most neuteral charge and the pH at which it appears
        mostNeuteral = 1000
        returnValue = 0
        for i in range (0, 1401, 1):
            i = i/100
            # Compare charges based on absolute value of charge
            if self._charge_(i) > 0 and self._charge_(i) < mostNeuteral:
                mostNeuteral = self._charge_(i)
                returnValue = i
            if self._charge_(i) < 0 and (-1*self._charge_(i)) < mostNeuteral:
                mostNeuteral = (-1*self._charge_(i))
                returnValue = i
        return returnValue

    def aaComposition (self) :
        return self.composition

    def _charge_ (self, pH):
        '''Add up all positive charges, then subtract all negative charges'''
        netCharge = 0
        # Add/subtract N and C terminus charges 
        netCharge += (10**self.aaNterm)/(10**self.aaNterm+10**pH)
        netCharge -= (10**pH)/(10**self.aaCterm+10**pH)
        # Add/subtract charges from charged amino acids
        for key, value in self.aa2chargePos.items():
            netCharge += self.composition[key] * ((10**value)/(10**value+10**pH))
        for key, value in self.aa2chargeNeg.items():
            netCharge -= self.composition[key] * ((10**pH)/(10**value+10**pH))
        return netCharge

    def molarExtinction (self):
        '''Add up extinction coefficient by multiplying coefficient table values by counts of respective aa's'''
        returnValue = 0
        for key, value in self.aa2abs280.items():
            returnValue += self.composition[key] * value
        return returnValue

    def massExtinction (self):
        '''Divide molar extinction by molecular weight, preventing division by 0'''
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''Sum weights of individual amino acids minus released water weight, one extra water weight unit to start'''
        sumWeight = self.mwH2O
        for key, value in self.composition.items():
            sumWeight += value * (self.aa2mw[key] - self.mwH2O)
        return sumWeight

########################################################################
# FastA reader 
########################################################################
import sys
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
# Nucleotide/amino acid counter 
########################################################################
class NucParams:
    ''' Class to count occurances of various amino acids, codons, and nucleotides in a genomic sequence.'''

    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        '''Constructor for new object, usually one object created per input instance.'''
        # Initialize amino acid dictionary
        self.aaComp = {
        'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
        'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
        'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
        'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0, '-': 0
        }
        # Initialize nucleic acid dictionary
        self.nucComp = {
        'A': 0, 'C': 0, 'T': 0,'G': 0, 'N': 0, 'U':0
        }
        # Initialize codon dictionary
        self.codonComp = {
        # U
        'UUU': 0, 'UCU': 0, 'UAU': 0, 'UGU': 0,  # UxU
        'UUC': 0, 'UCC': 0, 'UAC': 0, 'UGC': 0,  # UxC
        'UUA': 0, 'UCA': 0, 'UAA': 0, 'UGA': 0,  # UxA
        'UUG': 0, 'UCG': 0, 'UAG': 0, 'UGG': 0,  # UxG
        # C
        'CUU': 0, 'CCU': 0, 'CAU': 0, 'CGU': 0,  # CxU
        'CUC': 0, 'CCC': 0, 'CAC': 0, 'CGC': 0,  # CxC
        'CUA': 0, 'CCA': 0, 'CAA': 0, 'CGA': 0,  # CxA
        'CUG': 0, 'CCG': 0, 'CAG': 0, 'CGG': 0,  # CxG
        # A
        'AUU': 0, 'ACU': 0, 'AAU': 0, 'AGU': 0,  # AxU
        'AUC': 0, 'ACC': 0, 'AAC': 0, 'AGC': 0,  # AxC
        'AUA': 0, 'ACA': 0, 'AAA': 0, 'AGA': 0,  # AxA
        'AUG': 0, 'ACG': 0, 'AAG': 0, 'AGG': 0,  # AxG
        # G
        'GUU': 0, 'GCU': 0, 'GAU': 0, 'GGU': 0,  # GxU
        'GUC': 0, 'GCC': 0, 'GAC': 0, 'GGC': 0,  # GxC
        'GUA': 0, 'GCA': 0, 'GAA': 0, 'GGA': 0,  # GxA
        'GUG': 0, 'GCG': 0, 'GAG': 0, 'GGG': 0  # GxG
        }
        # Fill amino acid dictionary with original input string
        self.addSequence(inString)
        
    def addSequence (self, inSeq):
        '''Clean formatting of input data, then look up codons in either RNA or DNA table,
        add corresponding amino acids to composition dictionary,
        add every nucleotide, including N to the nucelotide composition dictionary,
        add every valid codon, in RNA formatting to the codon composition dictionary.'''
        # Remove whitespace in input, capitalize
        sequence = inSeq.upper()
        sequence = sequence.strip()

        # Parse through input in groups of three
        for i in range(0, len(sequence), 3):
            codon = sequence[i: i+3]
            # Fill amino acid dictionary one by one, excluding codons that code to nothing
            # In RNA table
            if(self.rnaCodonTable.get(codon) != None):
                self.aaComp[self.rnaCodonTable.get(codon)] += 1
            # In DNA table
            elif(self.dnaCodonTable.get(codon) != None):
                self.aaComp[self.dnaCodonTable.get(codon)] += 1
            # Invalid codon, don't translate or update aa count
            else:
                pass

            # Fill nucleotide dictionary
            for letter in codon:
                if(letter in 'ACGTUN'):
                    self.nucComp[letter] += 1

            # Fill codon dictionary
            # Filter out invalid codons
            if (codon.find('N') != -1):
                pass
            # Convert codon to RNA, add to codon table
            rnaCodon = codon.replace('T', 'U')
            if(self.codonComp.get(rnaCodon) != None):
                self.codonComp[rnaCodon] += 1
        
    def aaComposition(self):
        '''Return object's aa dictionary'''
        return self.aaComp
    def nucComposition(self):
        '''Return object's nuc dictionary'''
        return self.nucComp
    def codonComposition(self):
        '''Return object's codon dictionary'''
        return self.codonComp
    def nucCount(self):
        return sum(self.nucComp.values())
