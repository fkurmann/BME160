#!/usr/bin/env python3
# Name: Fabrice Kurmann (fkurmann)
# Group Members: James Arnold, Ronnie Voskoboynik, Tucker Morrison, Yahaira Uribe

''' This program returns key pieces of information regarding a sequence of amino acids the user enters through standard input.

Inspection Info Question Answers:

The data is saved in a dictionary mapping the number of each amino acid to their respective names.

I was able to quickly get the counts of individual amino acids as well as the total count of amino acids through simple key value calls to the dictionary or loops over the elements of the dictionary.

For charge, I looped over the dictionaries containing positively and negatively charged amino acids and multiplied these values by the terms specified in the equation given. Ph was just another parameter to the function.

Since range for loops do not allow for floating point numbers, I multiplied the entire range by 100 and indexed by 1. Then inside the loop, I divided the indexing variable by that same 100.

I called the molarExtinction function in the mass extinction function to prevent having to write the code multiple times.

I made use of all of the provided dictionaries such that the only one I had to build myself was one mapping amino acid prefixes to amino acid counts.
'''

class ProteinParam :
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

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()